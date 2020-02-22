from base import *
import random
from fractions import Fraction

def Lift(P,Q,meja = 20):
    """
    Opis:
       Funkcija Lift sprejme tocki P,Q in ju
       vlozi v $\Z$, tako da so koordinate
       modulo p še vedno enake.
       

     Definicija:
       Lift(P,Q,meja = 20)

     Vhodni podatki:
       P...razred Point(ElipticCurve), ki
           predstavlja tocko na elipticni
           krivulji in je definiran v dodatku base
       Q...razred Point(ElipticCurve), ki
           predstavlja tocko na elipticni
           krivulji in je definiran v dodatku base
       meja...predsatvlja zgornjo mejo iskanja
           nakljucnega stevila v algoritmu

     Izhodni  podatek:
       Tocki P,Q vlozeni v $\Z$
    """
    p = P.mod  
    x1 = P.x + random.randint(0,meja)*p
    x2 = Q.x + random.randint(0,meja)*p
    if P.x != Q.x:
        y1 = P.y + random.randint(0,meja)*p
        razl = abs(x2-x1)
        i=0
        while True:
            y2 = Q.y + i*p
            if pow(y2,2,razl) == pow(y1,2,razl):
                break
            i += 1
        A1 = int((y2**2-y1**2)/(x2-x1)/
                 -(x2**3-x1**3)/(x2-x1))
        B1 = y1**2-x1**3-A1*x1
    else:
        x2 = x1
        y1 = P.y + random.randint(0,meja)*p
        A1 = A + random.randint(0,meja)*p
        B1 = y1**2-x1**3-A1*x1
        if P.y == Q.y:
            y2 = y1
        else:
            y2 = -y1
    return(Point(A1,B1,p**2,x1,y1),Point(A1,B1,p**2,x2,y2))

def AnomalneAttack(P,Q):
    """
    Opis:
       Funkcija AnomalneAttack izvede resuje
       problem diskretnega logaritma
       $$kP = Q$$
       nad anomalnimi krivuljami
       

     Definicija:
       AnomalneAttack(P,Q)

     Vhodni podatki:
       P...razred Point(ElipticCurve), ki
           predstavlja tocko na elipticni
           krivulji in je definiran v dodatku base
       Q...razred Point(ElipticCurve), ki
           predstavlja tocko na elipticni
           krivulji in je definiran v dodatku base

     Izhodni  podatek:
        diskretni logaritem tocke $P$
    """
    (P1,Q1) = Lift(P,Q)
    p = P.mod
    P2 = (p-1)*P1
    Q2 = (p-1)*Q1
    m1 = p*Fraction(P2.y-P1.y,P2.x-P1.x)
    m2 = p*Fraction(Q2.y-Q1.y,Q2.x-Q1.x)
    k = m1/m2
    k1 = k.numerator
    k2 = k.denominator
    k2 = NumberMod(k2,p).inverse().num
    k = (k1*k2) % p
    return k
