import Bio.Seq
from Bio import SeqIO

#record = SeqIO.read("test_fasta.txt", "fasta")
for record in SeqIO.parse("test_fasta.txt", "fasta"):
    print("%s %i" % (record.id, len(record)))
    print(record.seq)
    print(record.seq.reverse_complement())


def kmers(s,k):
    wynik=list()
    for i in range(0,len(s)-k):
        wynik.append(s.seq[i:i+k])
        wynik.append(s.seq.reverse_complement()[i:i+k])
    return set(wynik)
wynik1=kmers(record, 5)
print(wynik1)

def euler(kmers):
    slownik=()
    for kmer in kmers:
        slownik(kmer)=