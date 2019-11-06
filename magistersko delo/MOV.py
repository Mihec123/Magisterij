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
       P...razred Point, ki predstavlja tocko na
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


def BabyGiantCount(E):
    q = E.mod
    mejasp = q+1-2*(q**0.5)
    mejazg = q+1+2*(q**0.5)
    seznamst = [i for i in range(int(mejasp),int(mejazg+1))]
    redi = []
    
    #vsaj en red rabmo
    T = E.rand()
    red = BabyGiant(T)
    redi.append(red)
    lcm = red

    while True:
        T = E.rand()
        red = BabyGiant(T)
        redi.append(red)
        lcm = lcm*red/gcd(lcm,red)
        st = 0
        for el in seznamst:
            if el%lcm == 0 and st ==0:
                N = el
                st += 1
            elif el%lcm == 0 and st >0:
                st += 1
                break

        if st == 1:
            return N
    
    


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
       P...razred Point, ki predstavlja tocko na
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


def naklon(P,Q):
    """
    Opis:
       Funkcija naklon izracuna naklon premice med tockama
       P in Q, ki lezita na elipticni krivulji

     Definicija:
       naklon(P,Q)

     Vhodni podatki:
       P...razred Point, ki predstavlja tocko na
           elipticni krivulji
       Q...razred Point, ki predstavlja tocko na
           elipticni krivulji   

     Izhodni  podatek:
       stevilo po modolu P.mod, ki predstavlja naklon
    """
    mod = P.mod
    if P == Q:
        st = (3*(P.x**2)+P.a) % mod
        im = NumberMod(2*(P.y),mod).inverse().num
        rez = (st*im) % mod
        return rez
    else:
        st = (Q.y-P.y) % mod
        im = NumberMod(Q.x-P.x,P.mod).inverse().num
        rez = (st*im) % mod
        return rez



def g(P,Q,X):
    """
    Opis:
       Funkcija g del funkcije potrebne za izracun Weilovega
       parjenja, s pomocjo Millerjevega algoritma. Funkcija izracuna
       funkcijo $g_{P,Q}(X)$

     Definicija:
       g(P,Q,X)

     Vhodni podatki:
       P...razred Point, ki predstavlja tocko na
           elipticni krivulji
       Q...razred Point, ki predstavlja tocko na
           elipticni krivulji
       X...razred Point, ki predstavlja tocko na
           elipticni krivulji

     Izhodni  podatek:
       stevilo, ki predtavlja vrednost $g_{P,Q}(X)$
    """
    if P.x == Q.x and P.y != Q.y:
        rez = (X.x-P.x) % P.mod
        return rez
    else:
        #lambda je naklon premice
        lam = naklon(P,Q)
        st = (X.y-P.y-lam*(X.x-P.x)) % P.mod
        im = NumberMod(X.x+P.x+Q.x-lam**2,P.mod).inverse().num
        rez = (st*im) % P.mod
        return rez

def Miller(P,X,m):
    """
    Opis:
       Funkcija Miller je implementacija Millerjevega algoritma
       potrebnega za izracun Weilovega parjenja. Ce je
       $$e_m(P,Q) = (f_P(Q+S)/f_P(S)) / (f_Q(P-S)/f_Q(-S)) $$
       potem funkcija Miller predstavlja izracun vrednosti
       $f_P(X)$.

     Definicija:
       Miller(P,X,m)

     Vhodni podatki:
       P...razred Point, ki predstavlja tocko na
           elipticni krivulji
       X...razred Point, ki predstavlja tocko na
           elipticni krivulji
       m...red tocke P

     Izhodni  podatek:
       vrednost funkcije $f_P(X)$
    """
    binarno = bin(m)[2:] #vrne niz porezemo 0b na zacetku niza
    n = len(binarno)
    T = P
    f = 1
    for i in range(1,n):
        f = (f**2 * g(T,T,X)) % P.mod
        T = 2*T
        if int(binarno[i]) == 1:
            f = (f* g(T,P,X)) % P.mod
            T = T + P
    return f


def WeilPairing(P,Q,S,N):
    """
    Opis:
       Funkcija WeilPairing je implementacija Weilovega
       parjenja $e_N(P,Q)$, kjer je
       $$e_N(P,Q) = (f_P(Q+S)/f_P(S)) / (f_Q(P-S)/f_Q(-S)) $$

     Definicija:
       WeilPairing(P,Q,S,N)

     Vhodni podatki:
       P...razred Point, ki predstavlja tocko na
           elipticni krivulji
       Q...razred Point, ki predstavlja tocko na
           elipticni krivulji
       S...razred Point, ki predstavlja tocko, ki
           ni v podgrupi generirani z P,Q
       N...red tocke P

     Izhodni  podatek:
       int k, ki predstavlja red tocke P (kP = $\infty$)
    """
    
    fpQS = Miller(P,Q+S,N)
    fpS = Miller(P,S,N)
    fqPS = Miller(Q,P-S,N)
    fqS = Miller(Q,-S,N)

    fpS = NumberMod(fpS,P.mod).inverse().num
    eN1 = (fpQS * fpS) % P.mod
    
    fqS = NumberMod(fqS,P.mod).inverse().num
    eN2 = (fqPS * fqS) % P.mod
    eN2 = NumberMod(eN2,P.mod).inverse().num

    eN  = (eN1*eN2) % P.mod

    return eN
    

P = Point(30,34,631,617,5)
Q = Point(30,34,631,121,244)
S = Point(30,34,631,0,36)

eN = WeilPairing(P,Q,S,5)
    
    


def MOV(E,m,N,P,Q):
    E1 = ElipticCurve(E.a,E.b,E.mod**m)
    while True:
        print("tuki")
        T = E1.rand()
        red = BabyGiant(T)
        d = gcd(red,N)
        print(d)
        T1 = (red//d)*T
        print(T1)
        print(BabyGiant(T1))
        while True:
            S = E1.rand()
            red1 = BabyGiant(S)
            print(red1)
            if red1 > red:
                break
        zeta1 = WeilPairing(P,T1,S,N)
        zeta2 = WeilPairing(Q,T1,S,N)
        break

    return (zeta1,zeta2)


def projektivnaAdd(P,Q):
    x1 = P.x
    x2 = Q.x
    y1 = P.y
    y2 = Q.y
    z1 = 1
    z2 = 1
    a = P.a
    b = P.b
    mod = P.mod
    x31 = (x1*y2-x2*y1)*(y1*z2+y2*z1)+ (x1*z2-x2*z1)*y1*y2 \
            -a*(x1*z2 + x2*z1)*(x1*z2-x2*z1)-3*b*(x1*z2-x2*z1)*z1*z2
    y31 = -3*x1*x2*(x1*y2-x2*y1)-y1*y2*(y1*z2-y2*z1)-a*(x1*y2-x2*y1)*z1*z2 \
          + a*(x1*z2+x2*z1)*(y1*z2-y2*z1)+3*b*(y1*z2-y2*z1)*z1*z2
    z31 = 3*x1*x2*(x1*z2-x2*z1)-(y1*z2+y2*z1)*(y1*z2-y2*z1) \
          + a*(x1*z2-x2*z1)*z1*z2

    x31 = x31 % mod
    y31 = y31 % mod
    z31 = z31 % mod

    x32 = y1*y2*(x1*y2+x2*y1)-a*x1*x2*(y1*z2+y2*z1)-a*(x1*y2+x2*y1)*(x1*z2+x2*z1) \
          -3*b*(x1*y2+x2*y1)*z1*z2 - 3*b*(x1*z2+x2*z1)*(y1*z2+y2*z1)\
          + a*a*(y1*z2+y2*z1)*z1*z2
    y32 = y1*y1*y2*y2 + 3*a*x1*x1*x2*x2+9*b*x1*x2*(x1*z2+x2*z1)\
          -a*a*x1*z2*(x1*z2+2*x2*z1)-a*a*x2*z1*(2*x1*z2+x2*z1)\
          -3*a*b*z1*z2*(x1*z2+x2*z1)-(a*a*a+9*b*b)*z1*z1*z2*z2
    z32 = 3*x1*x2*(x1*y2+x2*y1)+y1*y2*(y1*z2+y2*z1)+a*(x1*y2+x2*y1)*z1*z2\
          +a*(x1*z2+x2*z1)*(y1*z2+y2*z1)+3*b*(y1*z2+y2*z1)*z1*z2

    x32 = x32 % mod
    y32 = y32 % mod
    z32 = z32 % mod

    x33 = (x1*y2+x2*y1)*(x1*y2-x2*y1) + a*x1*x2*(x1*z2-x2*z1)\
          +3*b*(x1*z2+x2*z1)*(x1*z2-x2*z1)-a*a*(x1*z2-x2*z1)*z1*z2
    y33 = (x1*y2-x2*y1)*y1*y2-3*a*x1*x2*(y1*z2-y2*z1)\
          +a*(x1*y2+x2*y1)*(x1*z2-x2*z1)+3*b*(x1*y2-x2*y1)*z1*z2\
          -3*b*(x1*z2+x2*z1)*(y1*z2-y2*z1)+a*a*(y1*z2-y2*z1*z1*z2)
    z33 = -(x1*y2+x2*y1)*(y1*z2-y2*z1)-(x1*z2-x2*z1)*y1*y2 \
          -a*(x1*z2+x2*z1)*(x1*z2-x2*z1)-3*b*(x1*z2-x2*z1)*z1*z2

    x33 = x33 % mod
    y33 = y33 % mod
    z33 = z33 % mod

    tocka = ((x31+x32+x33)%mod,(y31+y32+y33)%mod,(z31+z32+z33)%mod)
    print([[x31,y31,z31],[x32,y32,z32],[x33,y33,z33]])

    try:
        inv = NumberMod(tocka[2],mod).inverse().num
        tocka = (tocka[0]*inv %mod,tocka[1]*inv% mod,1)
    except:
        pass
 

    #return [[x31,y31,z31],[x32,y32,z32],[x33,y33,z33]]
    return tocka

##E = ElipticCurve(1,1,13)
##P = E.rand()
##Q = E.rand()

P = Point(-1,1,25,1,1)
Q = Point(-1,1,25,21,4)
        

E = ElipticCurve(0,1,23)
P = Point(0,1,23,12,2)
m = 2
N = 24
Q = Point(0,1,23,1,18)
l = 12
E1 = ElipticCurve(0,1,23**m)
##a = MOV(E,m,N,P,Q)
a = Point(0,1,529,253,1)
    
##test = BabyGiant(Point(-10,21,557,2,3))
##print(test)
