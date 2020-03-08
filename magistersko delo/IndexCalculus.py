def IndexCalculus(g,h,q,baza):
    """
    Opis:
       IndexCalculus  vrne  tak $k$, da velja
       $$g^k \equiv h (mod q)$$, kjer je g generator
       multiplikativne grupe $\Z_q$.

       Zgled:
           g = 3
           q = 1217
           h = 37
           baza = [-1,2,3,5,7,11,13]

     Definicija:
       IndexCalculus(g,h,q,baza)

     Vhodni podatki:
       g...generator multiplikativne grupe $\Z_q$.
       q...modul s katerim delamo
       h...desna stran problema, ki ga resujemo
       baza...baza prastevil s katerimi delamo more
               vsebovati tudi -1

     Izhodni  podatek:
       $k$, da velja $$g^k \equiv h (mod q)$$
    """
    
    seznam_relacij = []
    V = VectorSpace(GF(q),len(baza)+1)
    Z = GF(q)
    k = 1
    while True:
        e0 = 0
        st = g^k % q
        if is_prime(st) and st not in baza:
            st = abs(st-q)
            e0 = 1
        tmp = ecm.factor(st)
        if set(tmp).issubset(set(baza)):
            relacija = [0]*(len(baza)+1)
            relacija[0] = e0
            for i in range(1,len(baza)):
                relacija[i] = tmp.count(baza[i])
            relacija[-1] = k

            seznam_relacij.append(relacija)
            M = Matrix(Z,seznam_relacij)
            if len(V.linear_dependence(M)) >0:
                seznam_relacij = seznam_relacij[:-1]

            if len(seznam_relacij) > len(baza):
                break
        k+=1
            
    tmp = seznam_relacij[:]
    b = [0]*(len(baza)+1)
    for i in range(len(baza)+1):
        if seznam_relacij[i][0] == 1:
            b[i] = seznam_relacij[i][-1] - ((q-1)//2)
        else:
            b[i] = seznam_relacij[i][-1]
        tmp[i] = tmp[i][1:-1]
    
    
    M = Matrix(IntegerModRing(q-1),tmp)
    B = vector(IntegerModRing(q-1),b)
    X = M.solve_right(B)
    
    j = 0
    while True:
        rez = (power_mod(g,j,q)*h) % q
        tmp = ecm.factor(rez)
        if set(tmp).issubset(set(baza)):
            enacba = [0]*(len(baza)-1)
            for i in range(1,len(baza)):
                enacba[i-1] = tmp.count(baza[i])
            enacba = vector(IntegerModRing(q-1),enacba)
            
            return( enacba*X - j)
        j+=1


#Uporaba funkcije na primeru 7.5
baza = [-1,2,3,5,7,11,13]
g = 3
q = 1217
h = 37
Lh = IndexCalculus(g,h,q,baza)
