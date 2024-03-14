#script to simulate evolution of philogenetic tree for data in fasta file

from Bio import SeqIO
from Bio.Seq import Seq
import random
from Bio import Phylo
import numpy as np

# Wczytanie sekwencji
protein_seq = next(SeqIO.parse("PAH.fa", "fasta")).seq

# Przypisanie losowego kodonu dla każdego aminokwasu
codon_table = {
    "A": ["GCT", "GCC", "GCA", "GCG"],
    "C": ["TGT", "TGC"],
    "D": ["GAT", "GAC"],
    "E": ["GAA", "GAG"],
    "F": ["TTT", "TTC"],
    "G": ["GGT", "GGC", "GGA", "GGG"],
    "H": ["CAT", "CAC"],
    "I": ["ATT", "ATC", "ATA"],
    "K": ["AAA", "AAG"],
    "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
    "M": ["ATG"],
    "N": ["AAT", "AAC"],
    "P": ["CCT", "CCC", "CCA", "CCG"],
    "Q": ["CAA", "CAG"],
    "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
    "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
    "T": ["ACT", "ACC", "ACA", "ACG"],
    "V": ["GTT", "GTC", "GTA", "GTG"],
    "W": ["TGG"],
    "Y": ["TAT", "TAC"],
    "*": ["TAA", "TAG", "TGA"]
}

dna_seq = Seq("")
for aa in protein_seq:
    codons = codon_table[aa]
    dna_codon = random.choice(codons)
    dna_seq += Seq(dna_codon)

print("Sekwencja nukleotydowa", dna_seq)
seq_len=len(dna_seq)

# Wczytanie drzewa
tree = Phylo.read("tree", "newick")

print(tree)
Phylo.draw_ascii(tree)

bases = ["A", "C", "G", "T"]
q=1/3*1/10000
p=9999/10000
P = np.array([[p, q, q, q],
              [q, p, q, q],
              [q, q, p, q],
              [q, q, q, p]])

def ile_roznic(s1, s2):
    count = 0
    for i in range(min(len(s1), len(s2))):
        if s1[i] != s2[i]:
            count += 1
    count += abs(len(s1) - len(s2))
    return count


seq_first_change=''

#Symulacja mutacji
for node in tree.get_nonterminals():
    kroki=int(node.branch_length/100*10000)

    if(len(tree.get_path(node))==0):
        seq=dna_seq
    elif(len(tree.get_path(node))==1):
        seq=seq_first_change
    else:
        seq=tree.get_path(node)[-2].name
    seq_krok=seq
    for step in range(kroki):
        seq_new=''
        for i in seq_krok:
            seq_new += np.random.choice(bases, p=P[bases.index(i)])
        seq_krok=seq_new
    node.name=seq_new
    print('roznice', ile_roznic(dna_seq, node.name)/seq_len)
    if(len(tree.get_path(node))==0):
        seq_first_change=seq_new


for node in tree.get_terminals():
    kroki=int(node.branch_length/100*10000)

    if(len(tree.get_path(node))==0):
        seq=dna_seq
    elif(len(tree.get_path(node))==1):
        seq=seq_first_change
    else:
        seq=tree.get_path(node)[-2].name
    seq_krok=seq
    for step in range(kroki):
        seq_new=''
        for i in seq_krok:
            seq_new += np.random.choice(bases, p=P[bases.index(i)])
        seq_krok=seq_new
    lisc=node.name
    node.name=seq_new
    print('roznice', ile_roznic(dna_seq, node.name)/seq_len)
    print(lisc)


from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo
import matplotlib.pyplot as plt

#Przetłumaczenie na sekwencje białek i zapisanie do pliku fasta
protein_results=[]
names=['A', 'B', 'C', 'D', 'E']
for i, lisc in enumerate(tree.get_terminals()):
    protein_results.append(SeqRecord(Seq(lisc.name).translate()))
    protein_results[i].id=names[i]
print(protein_results)

with open('bialka.fa', 'w') as f:
    SeqIO.write(protein_results, f, 'fasta')

in_file = "bialka.fa"
clustalw_cline = ClustalwCommandline("/home/julia/Downloads/uni/clustalw-2.1-linux-x86_64-libcppstatic/clustalw2",
                                     infile=in_file, matrix="BLOSUM", outfile="bialka_aln.fa", output="FASTA")

stdout, stderr = clustalw_cline()
print(stdout)
print(stderr)

alignment = AlignIO.read('bialka_aln.fa', "fasta")

# Obliczenie macierzy odległości
calculator = DistanceCalculator('blosum62')
dm = calculator.get_distance(alignment)

# Wygenerowanie drzewa UPGMA
constructor = DistanceTreeConstructor(calculator, 'upgma')
tree_upgma = constructor.build_tree(alignment)

# Wygenerowanie drzewa NJ
constructor = DistanceTreeConstructor(calculator, 'nj')
tree_nj = constructor.build_tree(alignment)

print(tree_upgma)
print(tree_nj)

#Rysunki i zapisanie do plików
Phylo.draw(tree_upgma, do_show=False)
plt.savefig('upgma_tree.jpg')
Phylo.draw(tree_nj, do_show=False)
plt.savefig('nj_tree.jpg')
