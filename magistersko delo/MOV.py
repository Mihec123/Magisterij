import math
import random
import numpy as np
from base import *


    
def BabyGiant(P):
    """
    Opis:
       Funkcija BabyGiant je implementacija algoritma
       Baby step, Giant step za iskanje reda tocke P

     Definicija:
       BabyGiant(P)

     Vhodni podatki:
       P...razred Point, ki predstavljatocko na
           elipticni krivulji

     Izhodni  podatek:
       int k, ki predstavlja red tocke P (kP = $\infty$)
    """
    q = P.mod
    Q = (q+1)*P
    m = int(q**(1/4))+1
    seznam = []

    for j in range(m+1):
        seznam.append(j*P)


    tmp = (2*m)*P
    for k in range(-m,m+1):
        T = Q + k*tmp
        if T in seznam:
            j = seznam.index(T)
            break
        elif -1*T in seznam:
            j = seznam.index(-1*T)
            break
    st1 = q+1+2*m*k -j
    st2 = q+1+2*m*k +j
    if st1*P == INFPoint:
        M = st1
    else:
        M = st2

    odg = pomozna(M,P)

    return odg
    


def pomozna(M,P):
    """
    Opis:
       Funkcija pomozna predstavlja rekurzivni del funkcije
       BabyGiant

     Definicija:
       pomozna(M,P)

     Vhodni podatki:
       M...stevilo, ki ga testiramo ce predstavlja red
           tocke P
       P...razred Point, ki predstavljatocko na
           elipticni krivulji

     Izhodni  podatek:
       int k, ki predstavlja red tocke P (kP = $\infty$)
    """
    faktorji = Factor(M)
    for el in faktorji:
        if (M//el)*P == INFPoint:
            rez = pomozna(M//el,P)
            return rez
        else:
            return M


##def MOV(E,m,N):
##    E1 = ElipticCurve(E.a,E.b,E.mod**m)
##    while True:
##        T = E1.rand()
##        red = BabyGiant(T)
##        d = gcd(red,N)
##        T1 = (M//d)*T
    
    
test = BabyGiant(Point(-10,21,557,2,3))
print(test)
