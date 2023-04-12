import numpy as nu
from scipy import signal

def spike_remove(serie):
    pp = 2.5 * serie.std()
    picos = signal.find_peaks(serie, prominence = pp, width = [1,3])[0]
    serie[picos] = -9900
    pp = 3.0 * serie.std()
    picos = signal.find_peaks(serie, height = serie.mean() + pp, width = 3)[0]
    serie[picos] = -9900
    return serie

def slft(serie, nt):
    x = nu.linspace(nt.min(),nt.max(),10000)
    dt = nt.max() - nt.min()
    pendiente, origen = nu.ma.polyfit(nt, serie, 1)
    media = (origen + pendiente * nt)
    seriec = serie - media
    interp = media + nu.zeros_like(nt)
    for i in range(1,13):
        ff = i / dt
        #coeficiente real
        a = 2 * (seriec * nu.cos(nt * 2 * nu.pi * ff)).mean()
        #coeficiente imaginario
        b = 2 * (seriec * nu.sin(nt * 2 * nu.pi * ff)).mean()
        interp += a * nu.cos(nt * 2 * nu.pi * ff) + b * nu.sin(nt * 2 * nu.pi * ff)
    donde = nu.where(serie.mask == True)
    serie[donde].mask = False
    serie[donde] = interp[donde]
    return serie
    
def localiza(serie, i):
    D = 8
    E = 8
    k = 0
    while sum(1 - serie[i-D:i+E+1].mask) < 24:
        if i-D < 0: D = i
        if not k%2: D +=1
        else: E += 1
        k += 1
    return D, E

def interpola_o(serie, nt, nne = False, spiker = True):
    if spiker and (serie <= 0).sum() > 50:
        serie = interpola_o(serie, nt, nne = nne, spiker = False)
    if spiker:
        serie = spike_remove(serie)
        serie = -spike_remove(-serie)
        
    serie = nu.ma.masked_array(serie, mask = (serie < 0) + (serie > 1))
    donde = nu.where(serie.mask)[0]
    for i in donde:
        D, E = localiza(serie, i)
        serie[i-D:i+E+1] = slft(serie[i-D:i+E+1], nt[i-D:i+E+1])
    return serie

    
def interpola(serie, nt, nne = False, pasadas = 1):
    for i in range(pasadas):
        serie = interpola_o(serie, nt, nne = nne)
        serie = 1 - interpola_o(1 - serie, nt, nne = nne)
    return serie

