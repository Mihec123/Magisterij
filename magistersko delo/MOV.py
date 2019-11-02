import math
import random
import numpy as np
from base import *


##def BabyGiant(P,q,a):
##    (x1,y1) = P
##    modul = x1.modul
##    unit = (NumberMod(0,modul),NumberMod(1,modul))
##    Q = CubicCurveMulti(P,q+1,a)
##    m = int(q**1/4)+1
##    seznam = [None]*(m+1)
##    seznam[0] = unit
##    for j in range(1,m+1):
##        seznam[j] = CubicCurveMulti(P,j,a)
##
##    print(seznam)
##
##    mP = CubicCurveMulti(P,2*m,a)
##    for k in range(-m,m+1):
##        tmp = CubicCurveMulti(mP,k,a)
##        tmp1 = CubicCurveSum(Q,tmp,a)
##        (x2,y2) = tmp1
##        tmp2 = (-x2,-y2)
##        if (tmp1 in seznam):
##            j = seznam.index(tmp1)
##            print(k)
##            break
##        elif (tmp2 in seznam):
##            j = seznam.index(tmp2)
##            print(k)
##            break
##    print(k)
##    stevilo1 = q+1+2*m*k-j
##    stevilo2 = q+1+2*m*k+j
##    P1 = CubicCurveMulti(P,stevilo1,a)
##    P2 = CubicCurveMulti(P,stevilo2,a)
##    if P1 == unit:
##        M = stevilo1
##    else:
##        M = stevilo2
##    print(P1)
##    print(P2)
##    print(M)
##
##    faktorji = [M] #ni ok
    
def BabyGiant(P):
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

    #ni se dokoncan

    odg = pomozna(M,P)

    return odg
    


def pomozna(M,P):
    faktorji = Factor(M)
    for el in faktorji:
        if (M//el)*P == INFPoint:
            rez = pomozna(M//el,P)
            return rez
        else:
            return M
        
    
test = BabyGiant(Point(-10,21,557,2,3))
print(test)
