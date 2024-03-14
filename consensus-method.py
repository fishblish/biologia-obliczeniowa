# script with implementation of consensus method for motif finding

import numpy as np
import math
import pandas as pd
from Bio import motifs

def consensus(seqs, k, pseudocount=0.0001, base=[1/4, 1/4, 1/4, 1/4]):
    wynikowe_motywy=[]

    for i in range(len(seqs[0])-k+1):
        flaga=np.zeros(len(seqs))
        flaga[0]=1
        seq_set=[seqs[0][i:(i+k)]]

        while 0 in flaga:
            mf=motifs.create(seq_set)
            indeks_wykorzystany, wybrana_seq=wybierz_nastepna(mf, seqs, flaga, k)
            flaga[indeks_wykorzystany]=1
            seq_set.append(wybrana_seq)

        wynikowe_motywy.append(motifs.create(seq_set))

    wynikowe_motywy=sorted(wynikowe_motywy, key=lambda i: -i.pssm.mean())
    for motyw in wynikowe_motywy[:5]:
        print('information content',motyw.pssm.mean())
        print('motyw',motyw)


def wybierz_nastepna(motif, seqs, flaga, k):
    mf_lo=motif.counts.normalize(pseudocounts=1/len(seqs)).log_odds()
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

#przykład
seqs=['ATCGT', 'CTGAT', 'CCCCC', 'AGCCT']

consensus(seqs, 3)
