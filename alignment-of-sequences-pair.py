# script to designate alignment for pair of sequences

from Bio.Seq import Seq
from Bio.Align import substitution_matrices
from Bio import Align


def get_split_positions(str, i_start):
    list=[]
    for i, letter in enumerate(str):
        if((i==0) & (letter!='*')):
            list.append(i+i_start)
        if((i>0) & (letter!='*') & (str[i-1]=='*')):
            list.append(i*3+i_start)
    return list

def uliniowienie(seq1, seq2, matrix = substitution_matrices.load("BLOSUM62")):

    #stworzenie trzech wersji tłumaczeń dla różnych ramek odczytu
    protein1=[]
    for i in range(0,3):
         rest=len(seq1[i:])%3
         if(rest!=0):
            protein1.append({'seq': seq1[i:-rest].translate(), 'index': i})
         else:
            protein1.append({'seq': seq1[i:].translate(), 'index': i})

    protein2 = []
    for i in range(0, 3):
        rest = len(seq2[i:]) % 3
        if (rest != 0):
            protein2.append({'seq': seq2[i:-rest].translate(), 'index': i})
        else:
            protein2.append({'seq': seq2[i:].translate(), 'index': i})

    #podział tak, aby kodon stopu nie był sekwencji do uliniowienia
    for i, el in enumerate(protein1):
        protein1[i]['index']=get_split_positions(el['seq'], el['index'])
        protein1[i]['seq']=el['seq'].split(sep='*')

    for i, el in enumerate(protein2):
        protein2[i]['index']=get_split_positions(el['seq'], el['index'])
        protein2[i]['seq']=el['seq'].split(sep='*')

    #usunięcie pustych sekwencji stworzonych przez split
    for i in range(3):
        while Seq('') in protein1[i]['seq']:
            protein1[i]['seq'].remove(Seq(''))
        while Seq('') in protein2[i]['seq']:
            protein2[i]['seq'].remove(Seq(''))

    #stworzenie alignera
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix=matrix
    aligner.extend_gap_score=-1/2
    aligner.open_gap_score=-1
    aligner.mode='local'
    score=-100000

    #szukanie uliniowienia
    for trans1 in range(3):
        for trans2 in range(3):
            for i1, el1 in enumerate(protein1[trans1]['seq']):
                for i2, el2 in enumerate(protein2[trans2]['seq']):
                    if(score<=aligner.score(el1, el2)):
                        score=aligner.score(el1, el2)
                        alignment=aligner.align(el1, el2)
                        frag1=el1
                        index_start1=protein1[trans1]['index'][i1]
                        index_start2=protein2[trans2]['index'][i2]
                        frag2=el2
    for i in alignment:
        print(i)
    nukl1=seq1[index_start1:(index_start1+alignment[0].shape[1]*3)]
    nukl2=seq2[index_start2:(index_start2+alignment[0].shape[1]*3)]
    print('WYNIK:', score, '\nuliniowienie:', nukl1, nukl2)
    return [nukl1, nukl2]


uliniowienie(Seq("ATTGGCGAGCACA"), Seq("GATAGATAGCTAGCTAGTA"))
