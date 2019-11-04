import math
import random

def razEvk(a, b,inverz = True):
    """extended quclidian algorithm can compute inverse"""
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
    """NumberMod(a,n = None) represents the number a modulo n. If n= None
    then the number behaives as normal otherwise it is represented modulo n
    in a way such that it is writen either with + or - depending on which number
    is smaller. Example b = NumberMod(15,17) will be represented as -2 instead of
    15. To perform any calculations between numbers they must have the same
    modulo."""
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
        """addition returns num()"""
        if self.modul != other.modul:
            raise Exception('numbers have different modul')
        else:
            novonum = NumberMod((self.num + other.num) % self.modul,self.modul)
            return novonum
        
    def __sub__(self,other):
        """substracting returns num()"""
        if self.modul != other.modul:
            raise Exception('numbers have different modul')
        elif self.modul == None:
            return NumberMod(self.num-other.num,self.modul)
        else:
            novonum = NumberMod((self.num - other.num) % self.modul,self.modul)
            return novonum
        
    def __mul__(self,other):
        """multiplaying returns num()"""
        if self.modul != other.modul:
            raise Exception('numbers have different modul')
        else:
            novonum = NumberMod((self.num * other.num) % self.modul,self.modul)
            return novonum
        
    def __truediv__(self,other):
        """multiplying by inverse of other returns num()"""
        if self.modul != other.modul:
            raise Exception('numbers have different modul')
        else:
            temp = other.inverse()
            novonum = self*temp
            return novonum

    def __lt__(self, other):
        """less then returns True or False in modulo arithemtic"""
        if self.num > 0 and other.num > 0:
            return self.num < other.num
        elif self.num < 0 and other.num > 0:
            return (self.modul + self.num) < other.num
        elif self.num > 0 and other.num < 0:
            return self.num < (other.modul + other.num)
        else:
            return self.num < other.num

    def __le__(self, other):
        """less or equal returns True or False in modulo arithemtic"""
        if self.num > 0 and other.num > 0:
            return self.num <= other.num
        elif self.num < 0 and other.num > 0:
            return (self.modul + self.num) <= other.num
        elif self.num > 0 and other.num < 0:
            return self.num <= (other.modul + other.num)
        else:
            return self.num <= other.num

    def __eq__(self, other):
        """equal returns True or False in modulo arithemtic"""
        return self.num == other.num

    def __ne__(self, other):
        """not equal returns True or False in modulo arithemtic"""
        return self.num != other.num

    def __gt__(self, other):
        """greater then returns True or False in modulo arithemtic"""
        if self.num > 0 and other.num > 0:
            return self.num > other.num
        elif self.num < 0 and other.num > 0:
            return (self.modul + self.num) > other.num
        elif self.num > 0 and other.num < 0:
            return self.num > (other.modul + other.num)
        else:
            return self.num > other.num

    def __ge__(self, other):
        """greater or equal then returns True or False in modulo arithemtic"""
        if self.num > 0 and other.num > 0:
            return self.num >= other.num
        elif self.num < 0 and other.num > 0:
            return (self.modul + self.num) >= other.num
        elif self.num > 0 and other.num < 0:
            return self.num >= (other.modul + other.num)
        else:
            return self.num >= other.num

    def __pow__(self,b):
        if not isinstance(b,int):
            raise Exception('power must be an integer')
        else:
            if b == 0:
                return NumberMod(1,self.modul)
            elif b == 1:
                return NumberMod(self.num,self.modul)
            else:
                if b%2 == 0:
                    return NumberMod(self.num*self.num,self.modul)**(b//2)
                else:
                    return NumberMod(self.num,self.modul)*(NumberMod(self.num,self.modul)**(b-1))

    
    def __str__(self):
        return str(self.num)
    
    def __repr__(self):
        return str(self.num)
    
    def inverse(self):
        """inverse of self.num by modulo self.modul returns num()"""
        if self.modul == None:
            raise Exception('modul equals None')
        n = self.modul
        if self.num < 0:
            a = self.num + n
        else:
            a = self.num
        if n % a == 0:
            raise Exception('inverse does not exist')
        else:
            return NumberMod(razEvk(a,n)[1],n)      
        
    def small(self,n):
        """correts self.num by changing it into the smallest value by modulo n
        i.e. if n= 31 and self.num = 30 it changes self.num into
        self.num = -1. It does this if (n-1)/2 is lees then self.num"""
        if n == None:
            pass
        else:
            if ((n-1)/2) < self.numS:
                self.numS = self.numS-n
            elif self.numS < 0 and -((n-1)/2) > self.numS:
                self.numS = self.numS+n
                
    def isPrime(self,eps = 1/100000):
        """Miller-Rabin test for primes with error eps"""
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
                if xnov % n == 1 and (xprej != n-1 and xprej != 1):
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
        if self.modul == other.modul:
            return NumberMod(gcd(self.num,other.num),self.modul)
        else:
            raise Exception('numbers have different modul')

    def pow(self,exp1):
        """computes NumberMod(a,b)^exp1 by squaring"""
        exp = NumberMod(exp1,self.modul)
        base = self
        ans = NumberMod(1,self.modul)
        while exp > NumberMod(0,self.modul):
            if exp.num % 2 == 0:
                exp = NumberMod(exp.num /2,self.modul)
                base = base*base
            else:
                exp = NumberMod(exp.num -1,self.modul)
                ans = ans*base

                exp = NumberMod(exp.num /2,self.modul)
                base = base*base
                                 
        return ans



def testPrimitiveRoot(g,m,razcepPhi,Phi):
    """
    Opis:
       testPrimitiveRoot  testira ali je element g generator grupe $\Z_{m}$. Pri tem moramo podati razcep $\fi(m-1)$ ter $\fi(m)$


     Zgled:
         Phi = 1216
         razcepPhi =[2,19]
         m = 1217
         g = 23
         odg = testPrimitiveRoot(g,m,razcepPhi,Phi)

     Definicija:
       testPrimitiveRoot(g,m,razcepPhi,Phi)

     Vhodni podatki:
       g...stevilo ki ga testiramo
       m...modul s katerim delamo
       razcepPhi...razcep stevila podan s seznamom dobljenega z Eulerjevo phi funkcijo (npr. [2,3] bi bil razcep stevila, ki se zapise kot 2**x*3**y)
       Phi...vrednost $\phi(m)$

     Izhodni  podatek:
       Boolean True/False
    """
    skupni = razEvk(g,m,False)
    if skupni[0] != 1:
        return False
    for el in razcepPhi:
        if NumberMod(g**(Phi//el),m) == NumberMod(1,m):
            return False
    return True

def gcd(a,b):
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
    def __str__(self):
        """Pri izpisu print nam vrne krivuljo v lepi obliki"""
        if self.a > 0 and self.b > 0:
            return "y^2 = x^3 + {0}x + {1} mod {2}".format(self.a,self.b,self.mod)
        elif self.a > 0 and self.b <0:
            return "y^2 = x^3 + {0}x - {1} mod {2}".format(self.a,abs(self.b),self.mod)
        elif self.a > 0 and self.b == 0:
            return "y^2 = x^3 + {0}x mod {1}".format(self.a,self.mod)
        elif self.a == 0 and self.b > 0:
            return "y^2 = x^3 + {0} mod {1}".format(self.b,self.mod)
        elif self.a == 0 and self.b == 0:
            return "y^2 = x^3 mod {0}".format(self.mod)
        elif self.a == 0 and self.b <0:
            return "y^2 = x^3 - {0} mod {1}".format(abs(self.b),self.mod)
        elif self.a < 0 and self.b >0:
            return "y^2 = x^3 - {0}x + {1} mod {2}".format(abs(self.a),self.b,self.mod)
        elif self.a < 0 and self.b == 0:
            return "y^2 = x^3 - {0}x mod {1}".format(abs(self.a),self.mod)
        elif self.a < 0 and self.b <0:
            return "y^2 = x^3 - {0}x - {1} mod {2}".format(abs(self.a),abs(self.b),self.mod)



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
            x = random.randint(0,self.mod)
            y2 = NumberMod(x,self.mod)**3 + NumberMod(self.a*x,self.mod) + NumberMod(self.b,self.mod)
            y2 = y2.num
            for i in range(self.mod+1):         
                y = (i**2) %self.mod
                if y == y2:
                    najdl = True
                    break
        return Point(self.a,self.b,self.mod,x,i)

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
        
        if P.a != self.a or P.b != self.b or P.mod != self.mod:
            return False
        else:
            x = NumberMod(P.x,self.mod)**3 + NumberMod(self.a*P.x,self.mod) + NumberMod(self.b,self.mod)
            y = (P.y**2) % self.mod
            return x.num==y

           
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
        if (self.x == INF or self.y == INF) and self.y != self.x:
            #preverimo ali je tocka v neskoncnosti, v tem primeru zahtevamo,
            #da sta obe koordianti INF
            raise Exception('Point is not a proper infinity')

    def __str__(self):
        """Za lepsi izpis tocke po klicu print"""
        if self.x != INF:
            return "({0},{1}) mod {2}".format(self.x,self.y,self.mod)
        else:
            return u"\u221e"

    def __repr__(self):
        """za lepsi izpis tocke po klicu return"""
        if self.x != INF:
            return "({0},{1}) mod {2}".format(self.x,self.y,self.mod)
        else:
            return u"\u221e"

    def __add__(self,Q):
        """Algoritem za sestevanje tock nad isto elipticno
        krivuljo."""
        infty = False
        if self.x == INF or Q.x == INF:
            infty = True
        
        if not(self.a == Q.a and self.b == Q.b and self.mod == Q.mod) and not infty:
            raise Exception('Points don\'t lie on the same curve, or have different modulus')

        if self.x == INF:
            return Q
        elif Q.x == INF:
            return self

        elif self.x != Q.x:
            m = NumberMod(Q.x-self.x,self.mod).inverse().num
            m = (m*(Q.y-self.y)) % self.mod
            x3 = (m**2-self.x-Q.x) % self.mod
            y3 = (m*(self.x-x3) - self.y) % self.mod
            return Point(self.a,self.b,self.mod,x3,y3)

        elif self.x == Q.x and self.y != Q.y:
            return Point(None,None,None,INF,INF)

        elif self.x == Q.x and self.y == Q.y and self.y != 0:
            m = NumberMod(2*self.y,self.mod).inverse().num
            m = (m*(3*(self.x**2)+self.a)) % self.mod
            x3 = (m**2-2*self.x) % self.mod
            y3 = (m*(self.x-x3) - self.y) % self.mod
            return Point(self.a,self.b,self.mod,x3,y3)
        
        else:
            return Point(None,None,None,INF,INF)

    def __neg__(self):
        return Point(self.a,self.b,self.mod,self.x,-self.y % self.mod)


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
        elif self.a == Q.a and self.b == Q.b and self.x == Q.x and self.y == Q.y and self.mod == Q.mod:
            return True
        else:
            return False

INFPoint = Point(None,None,None,INF,INF)




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
