import mip
import math
import random
import numpy as np
from base import *


def PrimeSearch(num,mod,base):
    """
    Opis:
       PrimeSearch  vrne  razcep na prafaktorje iz baze ce
       tak razcep obstaja, v nasprotnem primeru vrne None.

     Definicija:
       PrimeSearch(num,mod,base)

     Vhodni podatki:
       num...stevilo, ki ga hocemo razcepiti na prafaktorje
       mod...modul s katerim delamo
       base...baza prastevil s katerimi delamo

     Izhodni  podatek:
       seznam_prafaktorjev ,kjer prvi element pomeni
       ali je spredaj - ali + (ce je element 0
       potem je +, ce je element 1 potem je -)
       ce tak razcep obstaja(npr. [1,2,0,3,1] pri bazi
       [2,3,5,7] pomeni da se to stevilo
       zapise kot -2**2*3**0*5**3*7), sicer None
    """

    factors = [0]*len(base)
    num1 = NumberMod(num,mod)
    Prime = num1.isPrime()
    #ce je stevilo prastevilo pogledamo
    #ce lahko razcepimo stevilo - mod
    if Prime and num1.num not in base:
        num1 = NumberMod(-num,mod)
        
    stevilo = num1.num
    st = 0
    for el in base:
        while stevilo % el == 0:
            factors[st] += 1
            stevilo  = stevilo // el
        st += 1

    if Prime and num1.num not in base and stevilo == 1:
        odg = [1] + factors
    else:
        odg = [0] + factors

    if (odg != [1] +[0]*len(base) and
        odg != [0]*(len(base)+1) and
        stevilo == 1):
        return odg
    else:
        return None
    
def FindSystemMod(mod,g,base):
    """
    Opis:
       FindSystemMod  vrne  sistem enacb, ki jih moramo
       resiti, da resimo problem diksretnega logaritma.

     Definicija:
       FindSystemMod(mod,g,base)

     Vhodni podatki:
       g...generator multiplikativne grupe $\Z_mod$.
       mod...modul s katerim delamo
       base...baza prastevil s katerimi delamo

     Izhodni  podatek:
       par (A,b)...A je matrika sistema,
                   b je desna stran sistema
    """
    #kratnik doloca koliko enacb je v sistemu
    #(n pomeni dolzina_baze*n enacb)
    kratnik = 5
    #pove ali dodamo enacbo v sistem ali ne
    dodamo = False 
    #pove ali ze imamo pokrite vse spremenljivke
    #z dosedanjimi enacbami
    VSE = False
    #porazdelitev spremenljivk v posameznih
    #enacbah ki smo jih dodali
    mamo = [[0]*len(base) for i in range(kratnik*len(base))]
    #katere spremenljivke smo ze pokrili
    mamoSPR = [0]*len(base)
    st = 0#stevilo dodanih enacb
    A = []
    b = []
    for i in range(1,mod-1):
        #razcepimo stevilo na faktorje iz baze
        odg = PrimeSearch(g**i,mod,base)
        if odg != None:
            #pogledamo ce smo dobili kaksno
            #novo spremenljivko
            if VSE == False:
                for j in range(1,len(base)+1):
                    #zamik zarad +- 1 spredi
                    if odg[j] != 0 and mamoSPR[j-1] == 0:
                        #ce se nimamo pokrite te
                        #spremenljivke to enacbo
                        #vzamemo
                        mamoSPR[j-1] = 1
                        dodamo = True
                        
            elif VSE and len(A) < kratnik*len(base):
                #nimamo se dovolj velikega sistema
                #preverimo ce tako enacbo ze mamo
                dodamo = True
                tmp = list(map(bool,odg[1:]))
                for k in range(st):
                    if mamo[k] == tmp:
                        dodamo = False
            if dodamo:
                mamo[st] = list(map(bool,odg[1:]))
                st += 1
                A.append(odg[1:])
                b.append(i-((mod-1)//2)*odg[0])
                dodamo = False
                if np.prod(mamoSPR) == 1:
                    VSE = True
            

        if (np.prod(mamoSPR) != 0 and
            len(A) >= kratnik*len(base)):
            #prekinemo iskanje ce imamo dovolj enacb
            break

    
    return(A,b)
        
        

def SolveSystemMod(A,b,mod):
    """
Opis:
   SolveSystemMod  vrne  resitev sistema linearnih enacb
   podanih z matrikama A in b (Ax = b) modulo stevilo mod.
   Metoda za resevanje sistem pretvori v novo obliko tako,
   da vsako vrstico pretvori v novo enacbo in potem
   celotni sistem resi s pomocjo celostevilskega linearnega
   programa iz paketa mip

   Zgled:
       A = [[0	,1	,0	,0	,0	,0	],
            [2	,0	,0	,1	,0	,1	],
            [0	,0	,3	,0	,0	,0	],
            [1	,0	,2	,0	,0	,0	],
            [0	,0	,1	,0	,1	,0	],
            [0	,0	,0	,0	,0	,1	]]
       b = [1,-584,25,-578,-554,87]

       mod = 1216

       posamezno vrstico pretvori tako, da namesto x2 = 1
       dobimo vrstico x2 + mod*X1(dodatna spremenljivka) = 1
       (za vsako vrstico uvedemo novo dodatno spremenljivko)

 Definicija:
   SolveSystemMod(A,b,mod)

 Vhodni podatki:
   A...matrika podana kot [[...],[...],[...],...,[...]]
   b...desna stran sistema podana s seznamom
   mod...modulo glede na katerega resujemo sistem

 Izhodni  podatek:
   x  seznam vrednosti resitve sistema
"""
    n =len(A)
    dol_vrstice = len(A[0])
    Kopija = A[:]
    #pretvorimo matriko A v primerno obliko
    for i in range(n):
        tmp = [0]*n
        tmp[i] = mod
        Kopija[i] = Kopija[i] + tmp

    #naredimo linearni program(zanimajo nas samo
    #spremenljivke, ki predstavljajo bazo,
    #katere omejimo med 0 in mod-1)
    m = mip.Model()
    x = [ m.add_var(var_type=mip.INTEGER,lb=1,ub = mod-1)
          for i in range(dol_vrstice) ]
    #spremenljivk tok kukr jih je v vrstici
    for i in range(n):
        x.append(m.add_var(var_type=mip.INTEGER,lb = -10,ub = 10))
        #dodatna spremenljivka za vsako vrstico posebaj
    m.objective = mip.minimize(
        mip.xsum(0*x[j] for j in range(dol_vrstice)))
    for i in range(n):
        m += mip.xsum(
            Kopija[i][j]*x[j] for j in range(dol_vrstice + n)) == b[i]
    m.optimize()

    #zaokrozimo za lepsi izpis(vcasih
    #lahko pride resitev 0.999999999 namesto 1, kar je isto)
    resitev = [int(round(x[i].x)) for i in range(dol_vrstice)]
    
    return(resitev)


def IndexCalculus(g,mod,a,base):
    """
    Opis:
       IndexCalculus  vrne  tak $k$, da velja
       $$g^k \equiv a (mod mod)$$, kjer je g generator
       multiplikativne grupe $\Z_mod$.

       Zgled:
           g = 3
           mod = 1217
           a = 37
           base = [2,3,5,7,11,13]

     Definicija:
       IndexCalculus(g,mod,a,base)

     Vhodni podatki:
       g...generator multiplikativne grupe $\Z_mod$.
       mod...modul s katerim delamo
       a...desna stran problema, ki ga resujemo
       base...baza prastevil s katerimi delamo

     Izhodni  podatek:
       $k$, da velja $$g^k \equiv a (mod mod)$$
    """
    Temp = FindSystemMod(mod,g,base)
    A = Temp[0]
    b = Temp[1]
    resitev = SolveSystemMod(A,b,mod-1)
    k = 0

    for i in range(1,mod-1):
        faktorji = PrimeSearch((g**i)*a,mod,base)
        if  faktorji != None:
            k = sum([a*b for a,b in zip(faktorji[1:],resitev)])
            k+=faktorji[0]*((mod-1)//2) -i 
            return NumberMod(k,mod-1).num


#primer
base = [2,3,5,7,11,13]
g = 27
mod = 1217
k = IndexCalculus(g,mod,37,base)
