from numpy import *
import sys, time
from interpola import interpola

set_printoptions(precision=4, linewidth=240, suppress=True)
def moving_average(a, n=7):
    ret = cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

def gompertz(deltat, n, r, k):
    '''
    Curva de crecimiento de Gomperz, caso especial de la GL == 0
    '''
    if (k > 1e-6):
        return deltat * r * n * log(k / n)
    else:
        return deltat * r * n * log(1e-6 / n)

def brody(deltat, n, r, k):
    '''
    Brody o Monomolecular
    '''
    return deltat * r * (k - n)
    
def rk4(modelo, dx, y, r, k):
    #integrador Runge-Kutta 4
    k1 = modelo(dx, y, r, k)
    k2 = modelo(dx, y + k1 / 2, r, k)
    k3 = modelo(dx, y + k2 / 2, r, k)
    k4 = modelo(dx, y + k3, r, k)
    return y + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0

#simula en el tiempo la población en función de un K y r variables
def simulan(n0 = 1e-4, k = ones((100,)), r = 1e-2, modelo = brody):
    nt = k.shape[0]
    n    = zeros((k.shape[0],))
    ni   = zeros(((n.shape[0] - 1) * 24,))
    n[0] = n0
    ii = 0
    for i in range(1, n.shape[0]):
        m = n[i-1]
        for j in range(24):
            m = rk4(modelo, 1/24., m, r, k[i])
            ni[ii] = m
            ii += 1
        n[i] = m
    return n

def leeclima(loi, lai):
    """
    Read climate data series from an numpy array containing daily data
    """
    y,x = (int((round(lai,1) + 42.0)*10), int((round(loi,1) + 73) * 10))
    climaf = load("datos_clima_eras5_patnor.npz")
    t0 = time.mktime((2000,1,1,0,0,0,0,0,0))
    tiempo = (climaf["tiempo"] - t0) / (86400 * 365.25) + 2000
    clima = empty((tiempo.shape[0],5))
    clima[:,0] = tiempo

    clima[:,1] = climaf["temperaturas"][:,y,x]
    clima[:,2] = climaf["lluvias"][:,y,x]
    #clima[:,5] = climaf["ETP"][:,y,x]
    clima[:,4] = climaf["radiacion"][:,y,x] / 10000000.0
    clima[:,3] = climaf["evaporacion"][:,y,x]
    #print(tiempo[-10:])
    return (clima)

from scipy.special import logit, expit
from scipy.optimize import fmin_bfgs, fmin_l_bfgs_b

def convolclima(params, clima):
    transf = fft.fft(clima)
    from scipy.stats import gamma
    curva = gamma.pdf(arange(clima.shape[0]), params[0], scale = params[1])
    curva = fft.fft(curva)
    cct = curva * transf
    clest = fft.ifft(cct).real
    return clest

def klineal_clima(params, modelo, climat, clima, nt, t0, ajuste = True):
    inn = where(abs(climat - t0) < 1/365)[0][0]
    ini = inn - 400
    fin = where(abs(climat - nt[-1]) < 1/365)[0][-1] + 30
    estima = zeros_like(nt)
    clima = clima[ini:fin]
    climat = climat[ini:fin]
    #print("nt0",nt[0], "t0", t0, "n0", n0, "cl0", climat[0], "cl400", climat[400])
    n0, oot, pet, oor, per, oop, pep, c, e, r = params
    kt = expit(oot + pet * clima[:,0])
    kr = expit(oor + per * clima[:,3])
    climac = convolclima([c, e], clima[:,1])

    kp = expit(oop + pep * climac)
    k = kt * kp * kr
        
    n = nan_to_num(simulan(n0, k = k[400:], r = r))
    j = 0
    for i in range(nt.shape[0]):
        try:
            while nt[i] > climat[400:][j]:
                j+=1
            estima[i] = n[j]
        except IndexError:
            print("Error de indice",i,j)
            pass
    return estima, (climat[400:], kt[400:], kr[400:], kp[400:], n)


def ajuste(params, climat, clima, nt, ndvi, modelo = 5, t0 = 2000):
    aj0 = 0
    estima, Kc = klineal_clima(params, modelo, climat, clima, nt, t0)
    ajj = sum((ndvi - estima) ** 2) + aj0
    return ajj, estima, Kc


def randomfit(nt, ndvi, lerror, climav, climat, ventana = 46, i = 46, modelo = 5, t0 = 2000, sdmax = 0.1, iteraciones = 10, actualiza = False, hijos = 1024, rmxx = 0, rmnn = -1):
    
    from multiprocessing import Pool, TimeoutError
    pool = Pool(processes=8)
    
    print("inicio",ventana, i, t0)
    params = lerror[0][2]
    optimos = params
        
    
    bounds = array([[0,1],[-20,20],[-2,2],[-20,20],[-2,2],[-20,20],[-100,1000],[0,200],[0,200],[0,100]])
    sdp    = array([ 0.1,        1,  0.1,       1,  0.1,       1,       1,      1,     50,     1])
    
    
    
        
    ntt = nt[(i-ventana):i]
    ive = ndvi[(i-ventana):i]
    scc = ((ive - ive.mean())**2).sum()

    salida = []
    for i in range(len(lerror)):
        params = lerror[i][2]
        params[:] = clip(params, bounds[:,0], bounds[:,1])
        salida.append([i, params, pool.apply_async(ajuste, args = (params, climat, climav, ntt, ive, modelo, t0))])
        
    for i, params, si in salida:
        lerror[i][0], lerror[i][3], lerror[i][6] = si.get()
        lerror[i][4] = 1-lerror[i][0] / scc
        lerror[i][5] = ventana * log(lerror[i][0] / ventana)
        
        
    lerror.sort()
    ermax = lerror[-1][0] * 1.01
    pool.close()
    pool = Pool(processes=8)
    
    if lerror[0][4] < 0:
        iteraciones *= 2
        hijos *= 4
        sdp[0] = 0.25
    elif actualiza:
        sdp[0] = 1e-2

    #sdp[0] *= max(abs(ndvi[0] - lerror[0][3][0])/0.2, 1)
        
    j = 0
    while (j < iteraciones or lerror[0][4] <= rmxx or lerror[-1][4] <= rmnn) and j < 100:
        likerel = clip(exp((lerror[0][5] - lerror[-1][5]) / 2), 0, 1)
        print("\rJ %2d SDP %6.4f %6.4f dN0 %6.4f r20 %9.6f r2m %9.6f likerel %9.6f params: r %8.4f "%(j, sdmax, sdp[0], abs(ndvi[0] - lerror[0][3][0]), lerror[0][4], lerror[-1][4],
                                                                                              likerel, lerror[0][2][-1]))
        le = len(lerror)
        vueltas = 0
        while(len(lerror) < hijos) and vueltas < 32:
            salida = []
            
            for i in range(le):
                rr = random.random()
                if rr < exp((lerror[0][5] - lerror[i][5]) / 2):
                    sdd = sdp.copy()
                    for ii in range(10):
                        sdd[0] *= max(abs(ndvi[3] - lerror[i][3][3])/(0.02 * (ii + 1)), 1)
                    params = random.normal(lerror[i][2], sdd * sdmax)
                    params = clip(params, bounds[:,0], bounds[:,1])
                    salida.append([i, params, pool.apply_async(ajuste, args = (params, climat, climav, ntt, ive, modelo, t0))])
            
            for i, params, si in salida:
                    errori, estima, Kc = si.get() #ajuste(params, climat, climav, ntt, ive, modelo = modelo, t0 = t0, n0 = n0)
                    r2 = 1 - errori / scc
                    llk = ventana * log(errori / ventana)
                    try:
                        rlk = exp((lerror[255][5] - llk) / 2)
                        print("\rLen %5d %4d %4d r20 %9.4f r2i %9.4f llk %8.2f rlk %8.4f E0 %6.4f E1 %6.4f E2 %6.4f E3 %6.4f P: r %8.4f         "%(
                            len(lerror), i, vueltas, lerror[0][4], clip(r2, -100, 1), clip(llk, -10000, 10000), clip(rlk, 0, 1),
                            abs(ndvi[0] - lerror[i][3][0]), abs(ndvi[1] - lerror[i][3][1]), abs(ndvi[2] - lerror[i][3][2]), abs(ndvi[3] - lerror[i][3][3]), lerror[i][2][-1]), 
                              end = " ");sys.stdout.flush()
                        if rr < rlk:
                            lerror.append([errori, random.random(), params.copy(), estima, r2, llk, Kc])
                    except IndexError:
                        lerror.append([errori, random.random(), params.copy(), estima, r2, llk, Kc])
                    #else: print(end = "Rechaza")
            vueltas += 1
        lerror.sort()
        lerror = lerror[:256]
        j += 1
        
    pool.close()
    return lerror
    
    
import warnings


def main(ve = 46):
    """
    Main function of the process
    parameters:
    ve: window width
    """
    #import modulse containing the database management for the output
    from SME_sqlite import creabase, guardadato
    #load data from the last parameter in the command line
    a1 = load(sys.argv[-3])

    
    try:
        creabase(sys.argv[-3])
    except: pass
    
    ndvi = a1["cubo"][0,0]/10000
    nt = a1["nt"][0,0]
    ndvi = interpola(ndvi, nt)
    coords = float(sys.argv[-2]), float(sys.argv[-1])
    clima = leeclima(coords[0], coords[1])

    optimos = array([0.6134,  4.8331,  0.0129,  3.2775,  0.1159, -2.1898,  8.7557,  6.7986, 45.6596,  7.1825])
    aparams = zeros((nt.shape[0], 256, optimos.shape[0])) - 1e9
    ar2 = zeros((nt.shape[0], 256)) - 1e9
    aest = zeros((nt.shape[0], 256 * ve)) - 1e9
    
    lerror = []
    for i in range(8):
        lerror.append([88, i, optimos, 0, 0, 0, [[],[],[],[],[]]])

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        lerror = randomfit(nt, ndvi, lerror, clima[:,1:], clima[:,0], ventana = 46, i = 46, modelo = 5, sdmax = 1, t0 = nt[0] - 16/365, iteraciones=1)

    for sdi in linspace(0,-1.2,7):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            lerror = randomfit(nt, ndvi, lerror, clima[:,1:], clima[:,0], ventana = 46, i = 46, modelo = 5, sdmax = power(10, sdi), t0 = nt[0] - 16/365, iteraciones = 10, rmnn = lerror[len(lerror)//2][4])
        #optimos = lerror[0][2]

    ID = 0
    for i in range(len(lerror)):
        guardadato(sys.argv[-3], lerror[i], nt, ve, ve, ID)
        ID += 1
    
    for i in range(ve + 1, nt.shape[0]):
        for l in range(len(lerror)):
            lerror[l][2][0] = lerror[l][3][1]
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            lerror = randomfit(nt, ndvi, lerror, clima[:,1:], clima[:,0], ventana = ve, i = i, modelo = 5, sdmax = 0.1, t0 = nt[i-ve-1], iteraciones = 1, actualiza = True, rmnn = lerror[-1][4] - 0.5)
        for i in range(len(lerror)):
            guardadato(sys.argv[-3], lerror[i], nt, i, ve, ID)
            ID += 1

if __name__ == "__main__": main()
