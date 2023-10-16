#### Otrzymany wynik to 464 ####
# liczy się 19 minut
# idea: dla ustalonego k tworzymy widmo k-merów każdej sekwencji
# sprawdzamy czy w danym widmie jest jakiś k-mer, którego nie ma w żadnym innym widmie. Jeśli nie ma takiego elementu, to wybrane k jest za małe
import math
import numpy as np
from Bio import SeqIO
import datetime

#funkcja zwracająca widmo dla danej sekwencji i danego k
def widmo(s,k):
    wynik=list()
    for i in range(0,len(s)-k+1):
        wynik.append(s.seq[i:i+k])
    return set(wynik)


print(datetime.datetime.now())
sekwencje=[] #lista sekwencji z pliku yeast.fa
najdluzsza=0 #tu będzie zapisana długość najdłuższej sekwencji
najkrotsza=100000000 #tu będzie zapisana długość najkrótszej sekwencji
for record in SeqIO.parse("yeast.fa", "fasta"):
    sekwencje.append(record)
    if najdluzsza<len(record.seq):
        najdluzsza=len(record.seq)
    if najkrotsza>len(record.seq):
        najkrotsza=len(record.seq)
print('sekwencji jest', len(sekwencje), 'najdłuższa ma', najdluzsza, 'najkrotsza ma', najkrotsza)


sprawdzone=[0]*len(sekwencje) #w tej liście jako 0 oznaczamy sekwencje, dla których trzeba jeszcze znaleźć najmniejsze takie k, które działa
#jeśli dla jakiejś sekwencji istnieje unikalny fragment o długości x, to tym bardziej istnieje unikalny dłuższy fragment

def usun_przeciecia(lista, k):
    macierz_odwiedzin=np.diag([1]*len(lista)) #w tej macierzy oznaczamy pary sekwencji, których przecięcia już odjęliśmy
    for index1 in range(len(lista)):
        if sprawdzone[index1]==0:
            for index2 in range(len(lista)):
                if macierz_odwiedzin[index1][index2] == 0:
                    cap=lista[index1].intersection(lista[index2])
                    if len(cap)>0:
                        lista[index1] = lista[index1]-cap
                        lista[index2] = lista[index2]-cap
                    macierz_odwiedzin[index1][index2]=1
                    macierz_odwiedzin[index2][index1]=1

            if(len(lista[index1])==0):
                return index1   #index1 znaczy że to k nie działa dla sekwencji o tym indeksie, dla wcześniejszych działa
    return -1   #-1 znaczy że to k działa


#szukane k musi spełniać 4^k >= liczba sekwencji
k_down=math.ceil(math.log(len(sekwencje), 4))
k_up=najkrotsza

#szukane k jest pomiędzy k_down, a k_up (włącznie)
#szukamy takiego k, które działa, ale k-1 nie działa
while k_up-k_down>1:
    #ustal k
    k=math.ceil((k_up+k_down)/2)
    print('sprawdzam dla k =', k)
    # zrob widmo dla tego k
    widma=[]
    for sekwencja in sekwencje:
        widma.append(widmo(sekwencja, k))
    # sprawdź czy to k jest wystarczające, czy nie, na tej podstawie ogranicz przedział możliwych wartości k
    wynik = usun_przeciecia(widma,k)
    if wynik==-1:
        k_up=k
    else:
        k_down=k
        sprawdzone[0:wynik] = [k]*wynik

print('!!!szukane k to', k_up)

print(datetime.datetime.now())
