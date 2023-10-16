from Bio import SeqIO
from Bio.Blast import NCBIWWW
import os
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
import json
import numpy as np
import pandas as pd
import math
from Bio import motifs

#część pierwsza

#stworzenie bazy sekwencji
os.system('makeblastdb -in data/genes_e_coli_new.fa -dbtype nucl -out base_db')

#wywołanie metody BLAST
os.system('tblastn -out output_blast.xml -outfmt 5 -query data/protein_fragments.fa -db base_db')


result_handle = open("output_blast.xml")
blast_records = NCBIXML.parse(result_handle)

blast_records = list(blast_records)

base={}
for record in SeqIO.parse('data/genes_e_coli_new.fa', 'fasta'):
    base[record.id]=record

#zapisanie wyniku xml do pliku fasta
result_sequences = []

for record in blast_records:
    alignment_name=record.alignments[0].title.split()[-3]
    result_sequences.append(SeqIO.SeqRecord(base[alignment_name].translate().seq, id=alignment_name+' '+record.query.split()[0]))


SeqIO.write(result_sequences, "blast_results.fa", "fasta")


#część druga

proteins=[]
for record in SeqIO.parse('data/protein_fragments.fa', 'fasta'):
    proteins.append(record)


os.system('mkdir hmmer_results')

#wywołanie metody HMMER dla każdej sekwencji białkowej
for query in proteins:
    com="curl -L -H 'Expect:' -H 'Accept:application/json' -F hmmdb=pfam -F seq='"+str(query.seq)+"' https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan -o 'hmmer_results/hmmer_output_"+query.id+".json'"
    os.system(com)


#wczystanie pliku wynikowego json i stworzenie dataframe'u zero-jedynkowego
hits_dict={}
input_proteins=[]
output_proteins=[]

for query in proteins:
    f = open('hmmer_results/hmmer_output_'+query.id+'.json')
    json_f = json.load(f)
    hits=[]
    for hit in json_f['results']['hits']:
        hits.append(hit['taxid'])
        output_proteins.append(hit['taxid'])
    hits_dict[query.id]=hits
    input_proteins.append(query.id)
    f.close()



tab = pd.DataFrame(0, index=input_proteins, columns=output_proteins)

for key in hits_dict.keys():
    for value in hits_dict[key]:
        tab.loc[key, value]=1

tab.to_csv('hmmer_result_table.csv', sep=';')



#część trzecia
#implementacja metody consensus
def consensus(seqs, k, bg={'A':1/4, 'C':1/4, 'G':1/4, 'T':1/4}):
    wynikowe_motywy=[]

    for i in range(len(seqs[0])-k+1):
        flaga=np.zeros(len(seqs))
        flaga[0]=1
        seq_set=[seqs[0][i:(i+k)]]

        while 0 in flaga:
            mf=motifs.create(seq_set)
            indeks_wykorzystany, wybrana_seq=wybierz_nastepna(mf, seqs, flaga, k, bg)
            flaga[indeks_wykorzystany]=1
            seq_set.append(wybrana_seq)

        wynikowe_motywy.append(motifs.create(seq_set))

    wynikowe_motywy=sorted(wynikowe_motywy, key=lambda i: -i.pssm.mean())
    for motyw in wynikowe_motywy[:5]:
        print('information content',motyw.pssm.mean())
        print('motyw',motyw)





def wybierz_nastepna(motif, seqs, flaga, k, bg):
    mf_lo=motif.counts.normalize(pseudocounts=1/len(seqs)).log_odds(background=bg)
    max_result=-1000
    for ind, seq in enumerate(np.array(seqs)):

        if(flaga[ind]==0):
            szanse=mf_lo.calculate(seq)
            max_arg=np.argmax(szanse)

            if szanse[max_arg]>max_result:
                max_result=szanse[max_arg]
                wybrana_seq=seq[max_arg:(max_arg+k)]
                wykorzystane=ind

    return(wykorzystane, wybrana_seq)




blast_results=[]
for record in SeqIO.parse('blast_results.fa', 'fasta'):
    blast_results.append(record)


proms={}
for record in SeqIO.parse('data/proms_e_coli_fixed.fa', 'fasta'):
    proms[record.id]=record

#podział na grupy A i B
promsA=[]
promsB=[]

#wyznaczenie rozkładu tła dla obu grup
A_background=[0,0,0,0]

B_background=[0,0,0,0]

for blast_record in blast_results:
    if 'groupA' in blast_record.description:
        promsA.append(proms[blast_record.id])
        sekwencja=proms[blast_record.id].seq
        A_background=np.add(A_background,[sekwencja.count('A'), sekwencja.count('C'), sekwencja.count('G'), sekwencja.count('T')])

    if 'groupB' in blast_record.description:
        promsB.append(proms[blast_record.id])
        sekwencja=proms[blast_record.id].seq
        B_background=np.add(B_background,[sekwencja.count('A'), sekwencja.count('C'), sekwencja.count('G'), sekwencja.count('T')])

A_background=A_background/sum(A_background)
B_background=B_background/sum(B_background)

A_bg={'A': A_background[0], 'C': A_background[1], 'G': A_background[2], 'T': A_background[3]}
B_bg={'A': B_background[0], 'C': B_background[1], 'G': B_background[2], 'T': B_background[3]}

#wywołanie metody consensus
print('Znalezione motywy dla grupy A:')
consensus([str(i.seq) for i in promsA], 10, A_bg)

print('Znalezione motywy dla grupy B:')
consensus([str(i.seq) for i in promsB], 10, B_bg)

