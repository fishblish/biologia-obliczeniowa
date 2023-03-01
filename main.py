import Bio.Seq
from Bio import SeqIO
import matplotlib.pyplot as plt
import networkx as nx

#record = SeqIO.read("test_fasta.txt", "fasta")
for record in SeqIO.parse("test_fasta.txt", "fasta"):
    print("%s %i" % (record.id, len(record)))
    print(record.seq)
    print(record.seq.reverse_complement())



def kmers(s,k):
    wynik=list()
    for i in range(0,len(s)-k+1):
        wynik.append(s.seq[i:i+k])
        wynik.append(s.seq.reverse_complement()[i:i+k])
    return set(wynik)

def widmo(s,k):
    wynik=list()
    for i in range(0,len(s)-k+1):
        wynik.append(s.seq[i:i+k])
    return set(wynik)

def zrob_graf(widmo):
    dict={}
    for wierzcholek in widmo:
        lista_zwiazanych=[]
        dict[wierzcholek]=lista_zwiazanych
        for x in ['A', 'T', 'G', 'C']:
            if (wierzcholek[1:]+x in widmo):
                lista_zwiazanych.append(wierzcholek[1:]+x)
                dict[wierzcholek] = lista_zwiazanych
    return dict

wynik1=widmo(record, 4)
print('sekwencja to', record.seq)
print('widmo to',wynik1)
graf1=zrob_graf(wynik1)
print('graf to',graf1)

#rysowanie grafu!

G = nx.DiGraph(graf1)
nx.draw(G, with_labels=True, font_weight='bold')
plt.show()

def hamilton(graph, start_v):
  size = len(graph)
  # if None we are -unvisiting- comming back and pop v
  to_visit = [None, start_v]
  path = []
  while(to_visit):
    v = to_visit.pop()
    if v :
      path.append(v)
      if len(path) == size:
        break
      for x in set(graph[v])-set(path):
        to_visit.append(None) # out
        to_visit.append(x) # in
    else: # if None we are comming back and pop v
      path.pop()
  return path

sciezka=[]
for wierzcholek in graf1:
    print(wierzcholek)
    sciezka=hamilton(graf1, wierzcholek)
    if sciezka != []:
        print(sciezka)
        break
if sciezka==[]:
    print('nie znaleziono sciezki hamiltona')