import numpy as nu
import sqlite3 as sq
import os

def creabase(lugar):
    try:
        #Trata de conectar a la base de datos
        con = sq.connect("%s.db"%lugar)
    except IOError:
        print("Error, creando base por sistema")
        #de lo contrario la crea usando comandos del sistema
        f = os.popen("sqlite %s.db"%lugar)
        f.write(".tables\n.exit\n")
        con = sq.connect("%s.db"%lugar)
        
    cur = con.cursor()
    cur.execute('create table ajuste("ID" INT, "fecha" REAL,"r2" REAL, "like" REAL)')
    cur.execute('create table params("ID" INT, "fecha" REAL, "n0" REAL, "oot" REAL, "pet" REAL, "oor" REAL, "per" REAL, "oop" REAL, "pep" REAL, "c" REAL, "e" REAL, "r" REAL)')
    cur.execute('create table estimacion("ID" INT, "fecha" REAL, "nest", REAL)')
    cur.execute('create table K_clima("ID" INT, "fecha" REAL, "ktemp" REAL, "krad" REAL, "kprec" REAL, "nest" REAL)')
    con.commit()
    con.close()

def guardadato(lugar, ajparams, nt, ix, ve, ID):
    con = sq.connect("%s.db"%lugar)
    cur = con.cursor()

    #inidice del inicio de las simulaciones
    i0 = ix - ve
    
    #guardando datos de ajuste
    cur.execute("insert into ajuste values (%d, %f, %f, %f)"%(ID, nt[i0], ajparams[2], ajparam[5]))

    #guardando parametros
    cur.execute("insert into K_clima values (%d, %f, %f, %f, %f, %f, %f, %f, %f, %f)"%(ID, ajparams[2][0], ajparams[2][1],
                                                                                       ajparams[2][2], ajparams[2][3],
                                                                                       ajparams[2][4], ajparams[2][5],
                                                                                       ajparams[2][5], ajparams[2][7],
                                                                                       ajparams[2][8], ajparams[2][9]))
    
    #guardando estimacion de ndvi por imagen
    for i in range(ve):
        cur.execute("insert into estimacion values (%d, %f, %f, %f)"%(ID, nt[i0 + i], ajparams[3][i]))
    
    #guardando estimacion dia x dia
    for j in range(len(ajparams[6][0])):
        cur.execute("insert into K_clima values (%d, %f, %f, %f, %f, %f)"%(ID, ajparams[6][0][i], ajparams[6][1][i],
                                                                           ajparams[6][2][i], ajparams[6][3][i],
                                                                           ajparams[6][4][i]))
    con.commit()
    con.close()
    
