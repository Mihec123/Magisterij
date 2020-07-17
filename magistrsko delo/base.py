import math
import random

def razEvk(a, b,inverz = True):
    """
    Opis:
       funkcija razEvk predstavlja razsirjen
       Evklidov algoritem

     Definicija:
       razEvk(a, b,inverz = True)

     Vhodni podatki:
       a...prvo stevilo
       b...drugo stevilo
       inverz...ce zelimo racunati samo inverz
           ne potrebujemo dodatne vrednosti
           ki nam jo da razsirjen Evklidov
           algoritem

     Izhodni  podatek:
       stevili $gcd(a,b),x,y$ tako, da
       $ax+by = gcd(a,b)$
       ali v primeru ko je inverz True
       stevili $gcd(a,b),a^-1$
       
    """
    if inverz == False:
        s0, s1, t0, t1 = 1, 0, 0, 1
        while b != 0:
            q, a, b = a // b, b, a % b
            s0, s1 = s1, s0 - q * s1
            t0, t1 = t1, t0 - q * t1
        return  a, s0, t0
    else:
        s0, s1 = 1, 0
        while b != 0:
            q, a, b = a // b, b, a % b
            s0, s1 = s1, s0 - q * s1
        return  a, s0

class NumberMod:
    """
    Opis:
       Razred za predstavitev stevil po nekem
       modulu p

     Definicija:
       NumberMod(a,mod)

     Vhodni podatki:
       a...stevilo, ki ga hocemo predstaviti
       mod...modul po katerem delamo

     Izhodni  podatek:
       Razred stevila po modolu
    """
    def __init__(self,a, n=None):
        self.original = a
        self.modul = n
        if self.modul == None:
            self.num = a
        else:
            self.num = (a%self.modul)
            self.numS = self.num
            self.small(self.modul)
            
    def __neg__(self):
        return NumberMod(-self.num,self.modul)
            
    def __add__(self,other):
        """sestevanje stevil po modolu"""
        if self.modul != other.modul:
            raise Exception('numbers have different modul')
        else:
            novonum = NumberMod((self.num + other.num)\
                                % self.modul,self.modul)
            return novonum
        
    def __sub__(self,other):
        """odstevanje stevil po modolu"""
        if self.modul != other.modul:
            raise Exception('numbers have different modul')
        elif self.modul == None:
            return NumberMod(self.num-other.num,self.modul)
        else:
            novonum = NumberMod((self.num - other.num)\
                                % self.modul,self.modul)
            return novonum
        
    def __mul__(self,other):
        """mnozenje stevil po modolu"""
        if self.modul != other.modul:
            raise Exception('numbers have different modul')
        else:
            novonum = NumberMod((self.num * other.num)\
                                % self.modul,self.modul)
            return novonum
        
    def __truediv__(self,other):
        """mnozenje z inverzom stevila po modolu"""
        if self.modul != other.modul:
            raise Exception('numbers have different modul')
        else:
            temp = other.inverse()
            novonum = self*temp
            return novonum

    def __lt__(self, other):
        """manjse kot"""
        if self.num > 0 and other.num > 0:
            return self.num < other.num
        elif self.num < 0 and other.num > 0:
            return (self.modul + self.num) < other.num
        elif self.num > 0 and other.num < 0:
            return self.num < (other.modul + other.num)
        else:
            return self.num < other.num

    def __le__(self, other):
        """manjse ali enako"""
        if self.num > 0 and other.num > 0:
            return self.num <= other.num
        elif self.num < 0 and other.num > 0:
            return (self.modul + self.num) <= other.num
        elif self.num > 0 and other.num < 0:
            return self.num <= (other.modul + other.num)
        else:
            return self.num <= other.num

    def __eq__(self, other):
        """enakost"""
        return self.num == other.num

    def __ne__(self, other):
        """ne enakost"""
        return self.num != other.num

    def __gt__(self, other):
        """vecje kot"""
        if self.num > 0 and other.num > 0:
            return self.num > other.num
        elif self.num < 0 and other.num > 0:
            return (self.modul + self.num) > other.num
        elif self.num > 0 and other.num < 0:
            return self.num > (other.modul + other.num)
        else:
            return self.num > other.num

    def __ge__(self, other):
        """vecje ali enako"""
        if self.num > 0 and other.num > 0:
            return self.num >= other.num
        elif self.num < 0 and other.num > 0:
            return (self.modul + self.num) >= other.num
        elif self.num > 0 and other.num < 0:
            return self.num >= (other.modul + other.num)
        else:
            return self.num >= other.num

    def __pow__(self,b):
        p = self.modul
        if not isinstance(b,int):
            raise Exception('power must be an integer')
        else:
            if b == 0:
                return NumberMod(1,p)
            elif b == 1:
                return NumberMod(self.num,p)
            else:
                if b%2 == 0:
                    st = NumberMod(self.num*self.num,p)
                    return st**(b//2)
                else:
                    st = NumberMod(self.num,p)**(b-1)
                    return NumberMod(self.num,p)*st

    
    def __str__(self):
        return str(self.num)
    
    def __repr__(self):
        return str(self.num)
    
    def inverse(self):
        """inverz stevila po modolu"""
        if self.modul == None:
            raise Exception('modul equals None')
        n = self.modul
        if self.num < 0:
            a = self.num + n
        else:
            a = self.num
        if gcd(a,n) != 1:
            raise Exception('inverse does not exist')
        else:
            return NumberMod(razEvk(a,n)[1],n)      
        
    def small(self,n):
        """ce je stevilo veliko ga proba funkcija zapisati
        s predznakom minus"""
        if n == None:
            pass
        else:
            if ((n-1)/2) < self.numS:
                self.numS = self.numS-n
            elif self.numS < 0 and -((n-1)/2) > self.numS:
                self.numS = self.numS+n
                
    def isPrime(self,eps = 1/100000):
        """Miller-Rabinov test za prastevilskost
        z napako eps"""
        #Napaka posameznega koraka je manjsa kot 1/4
        k = int(math.log(eps)/math.log(1/4))+1
        r = 0
        n = self.num
        if n == 3:
            return True
        if n < 0:
            n += self.modul
        if n % 2 == 0:
            return False
        temp = n-1
        d = temp
        while temp % 2 == 0:
            r += 1
            d = int(d/2)
            temp = d
        i = 0
        while i < k:
            a = random.randint(2,n-2)
            xprej = pow(a,d,n)
            j = 1
            while j < r+1:
                xnov = pow(xprej,2,n)
                if (xnov % n == 1 and
                (xprej != n-1 and xprej != 1)):
                    #xprej je Millerjeva prica
                    return False
                j += 1
                xprej = xnov
            if xnov % n != 1:
                #xnov Fermatova prica
                return False
            i += 1
        return True
    
    def gcd(self,other):
        p = self.modul
        q = other.modul
        if p == q:
            return NumberMod(gcd(self.num,other.num),p)
        else:
            raise Exception('numbers have different modul')


def gcd(a,b):
    """
    Opis:
       Funkcija gcd poisce najvecji
       skupni delitelj stevil $a,b$

     Definicija:
       gcd(a,b)

     Vhodni podatki:
       a...prvo stevilo
       b...drugo stevilo

     Izhodni  podatek:
       najvecji skupni delitelj stevil
       $a,b$
    """
    if b == 0:
        return a
    else:
        return gcd(b,a % b)


class ElipticCurve:
    """
    Opis:
       Razred elipticnih krivulj v Weierstrassovi obliki
       modulo p
       y^2 = x^3 + ax+b mod p

     Definicija:
       ElipticCurve(a,b,mod)

     Vhodni podatki:
       a...parameter a v enacbi
       b...parameter b v enacbi
       mod...modul p po katerem delamo

     Izhodni  podatek:
       Razred elipticne krivulje
    """

    def __init__(self, a, b, modulo):
        self.a = a
        self.b = b
        self.mod = modulo
        elipticna = self.isEliptic()
        if not elipticna:
            raise Exception('Curve is not Eliptic')
    def __str__(self):
        """Pri izpisu print nam vrne
        krivuljo v lepi obliki"""
        if self.a > 0 and self.b > 0:
            return "y^2 = x^3 + {0}x + {1} mod {2}".format(
                self.a,self.b,self.mod)
        elif self.a > 0 and self.b <0:
            return "y^2 = x^3 + {0}x - {1} mod {2}".format(
                self.a,abs(self.b),self.mod)
        elif self.a > 0 and self.b == 0:
            return "y^2 = x^3 + {0}x mod {1}".format(
                self.a,self.mod)
        elif self.a == 0 and self.b > 0:
            return "y^2 = x^3 + {0} mod {1}".format(
                self.b,self.mod)
        elif self.a == 0 and self.b == 0:
            return "y^2 = x^3 mod {0}".format(
                self.mod)
        elif self.a == 0 and self.b <0:
            return "y^2 = x^3 - {0} mod {1}".format(
                abs(self.b),self.mod)
        elif self.a < 0 and self.b >0:
            return "y^2 = x^3 - {0}x + {1} mod {2}".format(
                abs(self.a),self.b,self.mod)
        elif self.a < 0 and self.b == 0:
            return "y^2 = x^3 - {0}x mod {1}".format(
                abs(self.a),self.mod)
        elif self.a < 0 and self.b <0:
            return "y^2 = x^3 - {0}x - {1} mod {2}".format(
                abs(self.a),abs(self.b),self.mod)

    def __repr__(self):
        return str(self)



    def rand(self):
        """
        Opis:
           Funkcija rand generira nakljucno tocko na
           elipticni krivulji E
           y^2 = x^3 + ax+b mod p

         Definicija:
           E.rand()

         Vhodni podatki:
           ni vhodnih podatkov

         Izhodni  podatek:
           Razred Point, ki predstavlja tocko na
           elipticni krivulji
        """
        najdl = False
        while not najdl:
            x = random.randint(0,self.mod-1)
            y2 = NumberMod(x,self.mod)**3\
                 + NumberMod(self.a*x,self.mod)\
                 + NumberMod(self.b,self.mod)
            y2 = y2.num
            if pow(y2,(self.mod-1)//2,self.mod) == 1:
                y = TonelliShanks(y2,self.mod)
                najdl = True
        return Point(self.a,self.b,self.mod,x,y)

    def isOn(self,P):
        """
        Opis:
           Funkcija isOn preveri ali tocka P res
           lezi na elipticni krivulji E

         Definicija:
           E.isOn(P)

         Vhodni podatki:
           P...razred Point, ki predstavlja tocko
               na elipticni krivulji

         Izhodni  podatek:
           True/False
        """
        
        if (P.a != self.a or P.b != self.b
            or P.mod != self.mod):
            return False
        else:
            x = NumberMod(P.x,self.mod)**3\
                + NumberMod(self.a*P.x,self.mod)\
                + NumberMod(self.b,self.mod)
            y = (P.y**2) % self.mod
            return x.num==y
    def isEliptic(self):
        """Preverimo ali je dana krivulja res elipticna"""
        if not (self.a ==None and self.b== None
                and self.mod == None):
            #ce to velja je tocka v neskoncnosti
            #nas ne zanima
            temp = (4*self.a**3 + 27*self.b**2) % self.mod
            if temp != 0:
                return True
            else:
                return False
        else:
            return True

           
INF = "inf"

class Point(ElipticCurve):
    """
    Opis:
       Razred tock na elipticni krivulji, ki je
       podrazred razdreda ElipticCurve.

       Tocka $\infty$ je podana kot
       INFPoint = Point(None,None,None,INF,INF),
       kjer je spremenljivka INF = "inf"

     Definicija:
       Point(a,b,mod,x,y)

     Vhodni podatki:
       a...parameter a v enacbi
       b...parameter b v enacbi
       mod...modul p po katerem delamo
       x...x-koordinate tocke
       y...y-koordinata tocke

     Izhodni  podatek:
       Razred tocke na elipticni krivulji
    """

    def __init__(self, a,b, mod, x, y):
        super().__init__(a,b, mod)
        self.x = x
        self.y = y
        if ((self.x == INF or self.y == INF)
            and self.y != self.x):
            #preverimo ali je tocka v neskoncnosti,
            # v tem primeru zahtevamo, da sta
            #obe koordianti INF
            raise Exception('Point not infinity')

    def __str__(self):
        """Za lepsi izpis tocke po klicu print"""
        if self.x != INF:
            return "({0},{1}) mod {2}".format(
                self.x,self.y,self.mod)
        else:
            return u"\u221e"

    def __repr__(self):
        """za lepsi izpis tocke po klicu return"""
        if self.x != INF:
            return "({0},{1}) mod {2}".format(
                self.x,self.y,self.mod)
        else:
            return u"\u221e"

    def __add__(self,Q):
        """Algoritem za sestevanje tock nad isto elipticno
        krivuljo."""
        infty = False
        if self.x == INF or Q.x == INF:
            infty = True
        
        if (not(self.a == Q.a and self.b == Q.b
                and self.mod == Q.mod)
            and not infty):
            raise Exception("""Points don\'t lie on the
                same curve, or have different modulus""")

        if self.x == INF:
            return Q
        elif Q.x == INF:
            return self

        elif self.x != Q.x:
            m = NumberMod(Q.x-self.x,self.mod).inverse().num
            #print("stevec: ",Q.y-self.y)
            m = (m*(Q.y-self.y)) % self.mod
            x3 = (m**2-self.x-Q.x) % self.mod
            y3 = (m*(self.x-x3) - self.y) % self.mod
            return Point(self.a,self.b,self.mod,x3,y3)

        elif self.x == Q.x and self.y != Q.y:
            return Point(None,None,None,INF,INF)

        elif (self.x == Q.x and self.y == Q.y
              and self.y != 0):
            m = NumberMod(2*self.y,self.mod).inverse().num
            m = (m*(3*(self.x**2)+self.a)) % self.mod
            x3 = (m**2-2*self.x) % self.mod
            y3 = (m*(self.x-x3) - self.y) % self.mod
            return Point(self.a,self.b,self.mod,x3,y3)
        
        else:
            return Point(None,None,None,INF,INF)

    def __neg__(self):
        if self == INFPoint:
            return INFPoint
        else:
            return Point(self.a,self.b,self.mod,
                         self.x,-self.y % self.mod)


    def __sub__(self,Q):
        return self + (-Q)

    

    def __rmul__(self,num):
        """Funkcija predstavlja mnozenje stevila z tocko
        npr. 5P = P+P+P+P+P"""
        if not isinstance(num,int):
            raise Exception('Can only multiply with an int')

        if num < 0:
            return -(abs(num)*self)

        elif num == 0:
            return Point(None,None,None,INF,INF)
        elif num == 1:
            return self
        else:
            if num % 2 == 0:
                pol = int(num /2)
                return (pol*self + pol*self)
            else:
                return (self+(num-1)*self)

    def __eq__(self, Q):
        """Dve tocki sta enaki, ce lezita na
        isti krivulji in imata iste koordiante
        ali pa predstavljata tocko v neskoncnosti"""
        if self.x == INF and Q.x == INF:
            return True
        elif (self.a == Q.a and self.b == Q.b
              and self.x == Q.x and self.y == Q.y
              and self.mod == Q.mod):
            return True
        else:
            return False

INFPoint = Point(None,None,None,INF,INF)


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
    
    return M

def Factor(n):
    """
    Opis:
       Factor je enostavna neucinkovita funkcija
       za iskanje fakotrojev stevila
     Definicija:
       Factor(n)
     Vhodni podatki:
       n...stevilo, ki ga hocemo faktorizirati
     Izhodni  podatek:
       Seznam faktorjev stevila n
    """
    faktorji = []

    tmp = NumberMod(n,None).isPrime()
    if tmp:
        faktorji.append(n)
        return faktorji
    else:
        tmp = n
        while tmp %2 == 0:
            faktorji.append(2)
            tmp = tmp//2
        i = 3
        while i*i < tmp:

            while tmp % i == 0:
                faktorji.append(i)
                tmp = tmp//i

            i+=2

        if tmp != 0:
            faktorji.append(tmp)

    return faktorji

def TonelliShanks(n,p1):
    """
    Opis:
       Funkcija TonelliShanks izracuna kvadraticni
       ostanek stevila $n$ modulo $p1$

     Definicija:
       TonelliShanks(n,p1)

     Vhodni podatki:
       n...stevilo katerega kvadraticni
           ostanek nas zanima
       p1...modul po katerem delamo

     Izhodni  podatek:
       kvadraticni ostanek stevila $n$,
       ce ta obstaja
    """
    p = p1-1
    s = 0
    while True:
        if p %2 == 0:
            s+=1
            p = p//2
        else:
            break
    Q = p
    
    while True:
        z = random.randint(0,p1-1)
        if pow(z,(p1-1)//2,p1) == p1-1:
            break
    M = s
    c = pow(z,Q,p1)
    t = pow(n,Q,p1)
    R = pow(n,(Q+1)//2,p1)

    while True:
        if t == 0:
            return 0
        elif t == 1:
            return R
        else:
            i = 1
            while i < M:
                if pow(t,pow(2,i),p1) == 1:
                    break
                i += 1
            b = pow(c,pow(2,M-i-1),p1)
            M = i
            c = pow(b,2,p1)
            t = (t*pow(b,2,p1)) % p1
            R = (R*b) % p1

#Zgled definicij
#Definirajmo elipticno krivuljo E
E = ElipticCurve(30,34,631)
#ter tocki P,Q na krivulji E
P = Point(30,34,631,36,60)
Q = Point(30,34,631,121,387)
#nakljucno tocko dobimo z naslednjim klicem
Nakljucna = E.rand()
