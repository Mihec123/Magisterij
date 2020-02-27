def MOV(E,q,m,P,Q,st_korakov = 10):
    """
    Opis:
       Funkcija MOV resuje problem diskretnega logaritma
       na krivulji E, podanega s tockama P,Q.
       Iscemo torej $k$, tako da velja $kP = Q$.

     Definicija:
       MOV(E,q,m,P,Q,st_korakov = 10)

     Vhodni podatki:
       E...Elipticna krivulja iz SAGE
       q...modul krivulje
       m...stevilo ki pove razsiritev prvotnega polja
           iz GF(q) gremo v GF(q^m) (GF = Galvajevo polje)
       P...tocka na elipticni krivulji iz SAGE
       Q...tocka na elipticni krivulji iz SAGE
       st_korakov...v algoritmu dobivamo rezultate po nekem
                    modulu in moramo na koncu s pomocjo
                    kitajskega izreka dobit koncni rezultat.
                    Lahko pa se zgodi da vedno dobivamo po
                    istih modulih in se program ne bi
                    koncal. st_korakov omeji kolikokrat
                    se ta del zanke izvede.

     Izhodni  podatek:
       stevilo k
    """
    N = P.order()
    k = GF(q^m,'x')
    E1 = E.base_extend(k)#razsirimo krivuljo na vecje polje
    P1 = E1(P)#razsirit moramo tudi tocki
    Q1 = E1(Q)
    
    odgovori = [[],[]]
    koraki = 0
    
    while True:
        #vsak ta del nam da odgovor po modolu d
        T = E1.random_element()
        M = T.order()
        d = gcd(N,M)
        T1 = int(M/d)*T
        w1 = P1.weil_pairing(T1,N)
        print(w1)
        w2 = Q1.weil_pairing(T1,N)
        print(w2)
        if w1 != w2:
            #k = Naivno(w1,w2)
            k = discrete_log(w2, w1)
            odgovori[0].append(k)
            odgovori[1].append(d)
        if lcm(odgovori[1]) >= q:
            break
        elif koraki > st_korakov:
            break
        koraki+=1
    #uporabimo kitajski izrek o ostankih
    k = crt(odgovori[0],odgovori[1])
    return k
