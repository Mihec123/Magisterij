import mip
import math
import random
import numpy as np

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
        #Napaka posameznega koraka je manjša kot 1/4
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
                    #xprej je Millerjeva priča
                    return False
                j += 1
                xprej = xnov
            if xnov % n != 1:
                #xnov Fermatova priča
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


def PrimeSearch(num,mod,base):
    """
    Opis:
       PrimeSearch  vrne  razcep na prafaktorje iz baze ce tak razcep obstaja, v nasprotnem primeru vrne None.

     Definicija:
       PrimeSearch(num,mod,base)

     Vhodni podatki:
       num...stevilo, ki ga hocemo razcepiti na prafaktorje
       mod...modul s katerim delamo
       base...baza praštevil s katerimi delamo

     Izhodni  podatek:
       seznam_prafaktorjev ,kjer prvi element pomeni ali je spredaj - ali + (ce je element 0 potem je +, ce je element 1 potem je -)
       ce tak razcep obstaja(npr. [1,2,0,3,1] pri bazi [2,3,5,7] pomeni da se to stevilo zapise kot -2**2*3**0*5**3*7), sicer None
    """

    factors = [0]*len(base)
    num1 = NumberMod(num,mod)
    Prime = num1.isPrime() #ce je stevilo prastevilo pogledamo ce lahko razcepimo stevilo - mod
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

    if odg != [1] +[0]*len(base) and odg != [0]*(len(base)+1) and stevilo == 1:
        return odg
    else:
        return None
    
def FindSystemMod(mod,g,base):
    """
    Opis:
       FindSystemMod  vrne  sistem enacb, ki jih moramo resiti, da resimo problem diksretnega logaritma.

     Definicija:
       FindSystemMod(mod,g,base)

     Vhodni podatki:
       g...generator multiplikativne grupe $\Z_mod$.
       mod...modul s katerim delamo
       base...baza praštevil s katerimi delamo

     Izhodni  podatek:
       par (A,b)...A je matrika sistema, b je desna stran sistema
    """
    kratnik = 5 #dolocamo koliko enacb je v sistemu (n pomeni dolzina_baze*n enacb)
    dodamo = False #pove ali dodamo enacbo v sistem ali ne
    VSE = False #pove ali ze imamo pokrite vse spremenljivke z dosedanjimi enacbami
    mamo = [[0]*len(base) for i in range(kratnik*len(base))]#porazdelitev spremenljivk v posameznih enacbah ki smo jih dodali
    mamoSPR = [0]*len(base)#katere spremenljivke smo ze pokrili
    st = 0#stevilo dodanih enacb
    A = []
    b = []
    for i in range(1,mod-1):
        odg = PrimeSearch(g**i,mod,base)#razcepimo stevilo na faktorje iz baze
        if odg != None:
            #pogledamo ce smo dobili kaksno novo spremenljivko
            if VSE == False:
                for j in range(1,len(base)+1):#zamik zarad +- 1 spredi
                    if odg[j] != 0 and mamoSPR[j-1] == 0:#ce se nimamo pokrite te spremenljivke to enacbo vzamemo
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
            

        if np.prod(mamoSPR) != 0 and len(A) >= kratnik*len(base):#prekinemo iskanje ce imamo dovolj enacb
            break

    
    return(A,b)
        
        

def SolveSystemMod(A,b,mod):
    """
Opis:
   SolveSystemMod  vrne  rešitev sistema linearnih enacb
   podanih z matrikama A in b (Ax = b) modulo število mod.
   Metoda za reševanje sistem pretvori v novo obliko tako,
   da vsako vrstico pretvori v novo enacbo in potem celotni sistem
   resi s pomocjo celostevilskega linearnega programa iz paketa mip

   Zgled:
       A = [[0	,1	,0	,0	,0	,0	],
            [2	,0	,0	,1	,0	,1	],
            [0	,0	,3	,0	,0	,0	],
            [1	,0	,2	,0	,0	,0	],
            [0	,0	,1	,0	,1	,0	],
            [0	,0	,0	,0	,0	,1	]]
       b = [1,-584,25,-578,-554,87]

       mod = 1216

       posamezno vrstico pretvori tako, da namesto x2 = 1  dobimo vrstico
       x2 + mod*X1(dodatna spremenljivka) = 1 (za vsako vrstico uvedemo novo dodatno spremenljivko)

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

    #naredimo linearni program(zanima nas samo prvih 6 spremenljivk katere omejimo med 0 in mod-1)
    m = mip.Model()
    x = [ m.add_var(var_type=mip.INTEGER,lb=1,ub = mod-1) for i in range(dol_vrstice) ] #spremenljivk tok kukr jih je v vrstici
    for i in range(n):
        x.append(m.add_var(var_type=mip.INTEGER,lb = -10,ub = 10))#dodatna spremenljivka za vsako vrstico posebaj
    m.objective = mip.minimize(mip.xsum(0*x[j] for j in range(dol_vrstice)))
    for i in range(n):
        m += mip.xsum(Kopija[i][j]*x[j] for j in range(dol_vrstice + n)) == b[i]
    m.optimize()

    #zaokrozimo za lepsi izpis(vcasih lahko pride resitev 0.999999999 namesto 1, kar je isto)
    resitev = [int(round(x[i].x)) for i in range(dol_vrstice)]
    
    return(resitev)

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

def IndexCalculus(g,mod,a,base):
    """
    Opis:
       IndexCalculus  vrne  tak $k$, da velja $$g^k \equiv a (mod mod)$$, kjer je g generator
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
       base...baza praštevil s katerimi delamo

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
            k = sum([a*b for a,b in zip(faktorji[1:],resitev)]) + faktorji[0]*((mod-1)//2) -i 
            return NumberMod(k,mod-1).num

    
#primer1 iz magisterija
A = [[0	,1	,0	,0	,0	,0	],
            [2	,0	,0	,1	,0	,1	],
            [0	,0	,3	,0	,0	,0	],
            [1	,0	,2	,0	,0	,0	],
            [0	,0	,1	,0	,1	,0	],
            [0	,0	,0	,0	,0	,1	]]
b = [1,-584,25,-578,-554,87]

base = [2,3,5,7,11,13]
g = 27
mod = 1217
k = IndexCalculus(g,mod,37,base)

Phi1217 = 1216
base1216=[2,19]
odg = testPrimitiveRoot(23,mod,base1216,Phi1217)
