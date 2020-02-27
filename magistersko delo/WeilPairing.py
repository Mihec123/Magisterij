from base import *


def naklon(P,Q):
    """
    Opis:
       Funkcija naklon izracuna naklon premice med tockama
       P in Q, ki lezita na elipticni krivulji

     Definicija:
       naklon(P,Q)

     Vhodni podatki:
       P...razred Point(ElipticCurve), ki predstavlja
           tocko na elipticni krivulji in je definiran
           v dodatku base
       Q...razred Point(ElipticCurve), ki predstavlja
           tocko na elipticni krivulji in je definiran
           v dodatku base

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
       parjenja, s pomocjo Millerjevega algoritma.
       Funkcija izracuna funkcijo $g_{P,Q}(X)$

     Definicija:
       g(P,Q,X)

     Vhodni podatki:
       P...razred Point(ElipticCurve), ki predstavlja
           tocko na elipticni krivulji in je definiran
           v dodatku base
       Q...razred Point(ElipticCurve), ki predstavlja
           tocko na elipticni krivulji in je definiran
           v dodatku base
       X...razred Point(ElipticCurve), ki predstavlja
           tocko na elipticni krivulji in je definiran
           v dodatku base

     Izhodni  podatek:
       stevilo, ki predstavlja vrednost $g_{P,Q}(X)$
    """
    if P.x == Q.x and P.y != Q.y:
        rez = (X.x-P.x) % P.mod
        return rez
    else:
        #lambda je naklon premice
        lam = naklon(P,Q)
        st = (X.y-P.y-lam*(X.x-P.x)) % P.mod
        im = NumberMod(X.x+P.x+Q.x-lam**2,P.mod)
        im = im.inverse().num
        rez = (st*im) % P.mod
        return rez

def Miller(P,X,m):
    """
    Opis:
       Funkcija Miller je implementacija Millerjevega
       algoritma potrebnega za izracun Weilovega
       parjenja. Ce je
       $$e_m(P,Q)=(f_P(Q+S)/f_P(S))/(f_Q(P-S)/f_Q(-S))$$
       potem funkcija Miller predstavlja izracun
       vrednosti $f_P(X)$.

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
    binarno = bin(m)[2:]
    #vrne niz porezemo 0b na zacetku niza
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
       $$e_N(P,Q)=(f_P(Q+S)/f_P(S))/(f_Q(P-S)/f_Q(-S))$$

     Definicija:
       WeilPairing(P,Q,S,N)

     Vhodni podatki:
       P...razred Point(ElipticCurve), ki predstavlja
           tocko na elipticni krivulji in je definiran
           v dodatku base
       Q...razred Point(ElipticCurve), ki predstavlja
           tocko na elipticni krivulji in je definiran
           v dodatku base
       S...razred Point(ElipticCurve), ki predstavlja
           tocko, ki se ne nahaja v podgrupi
           generirani z P,Q
       N...veckratnik reda tocke P

     Izhodni  podatek:
       int eN, ki predstavlja vrednost Weilovega parjenja
    """

    #test ali je S v podgrupi generirani z P,Q
    
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

#Zgled uporabe na primeru 6.7
E = ElipticCurve(30,34,631)
P = Point(30,34,631,36,60)
Q = Point(30,34,631,121,387)
S = Point(30,34,631,0,36)
eN = WeilPairing(P,Q,S,5)
