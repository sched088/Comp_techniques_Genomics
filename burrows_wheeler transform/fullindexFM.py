"""
1) the array storing BWT
2) the array storing the reference positions to the beginning of the sequence
3) the array storing the rank of the nucleotides.

Report how many mapped reads and their starting position on chrY to text file 'mapping_fullindexFM.txt'.
"""

import numpy as np

file_in = 'bwtSample1.fa'
# p = 'TAGTAGGGGCTTACGTGCTTAGTAGGTTTTAATGTGCTTACTAGGAG'

# Import BWT text from file

def import_genes(file_in):
    print("here")
    with open(file_in) as f:
        header = f.readline().rstrip()
        # print("Header: " + header)
        seq = ''
        for l, i in enumerate(f):
            seq += i.rstrip()
        # data_1 = {header1: data}
        # print(len(data))123   `
        print(len(seq))
        # print(seq)
    return seq, header

# Reference positions array

def reference_pos(t):
    sa = sorted([(t[i:], i) for i in range(len(t))])
    # print(sa)
    return list(map(lambda x: x[1], sa))

# Nucleotide rank array
def rank_array(seq):
    seq_list_sort = sorted(seq)
    seq_sort = "".join(seq_list_sort)
    r_a_f = []

    # Rank for sorted (F)
    rank = 0
    i = 0
    print(seq_sort)
    print(seq)
    while i < len(seq_sort):
        # if i == 0:
        #     rank = 0
        #     r_a_f.append(rank)
        if seq_sort[i] == seq_sort[i-1]:
            rank += 1
            r_a_f.append(rank)
        else:
            rank = 0
            r_a_f.append(rank)
        i += 1
    # return r_a_f

    # Rank for seq (L)
    dict = {}
    r_a_l = []
    for _ in seq:
        if _ not in dict:
            dict[_] = 0
        r_a_l.append(dict[_])
        dict[_] += 1
        totals = dict
    r_a = np.vstack((r_a_f, r_a_l)).T
    return r_a, seq_sort, totals

def fifm(F, L, r_a, P, totals):


    l_list = []
    value = int

    # for index, nucleotide in enumerate(reversed(P)):
    #     for n in F:
    #         if nucleotide == n:
    #             f_list.append(n.index)
    #
    #     print(f_list)

    print(F)
    f_list = []
    for v, n in enumerate(reversed(P)):
        # print("N+1: " + str(P[v]) + str(P[v+1]))
        print(v,n)
        # if column == 'F':
        for index, nucleotide in enumerate(F):
            print("nf in F: " + str(nucleotide))
            print(f_list)
            if n == nucleotide:
                f_list.append(index)
                f_list.append(index + (totals[n]-1))
                break
        if not f_list: break
        print("f_list: " + str(f_list))

        # for i in f_list:
        #     print("L[i]: " + str(L[i]))
        #     print("= " + str(P[:(-1*v)-1]))
        #     if L[i] == P[:(-1*v)-1]:
        #         value = L[i]
        #         l_list.append(r_a[i][1]) #rank of the characters in l
        # print("l_list: " + str(l_list))
        # print(value)



        # if column == 'L':
        #     for i in f_list:
        #         if L[i] == n:
        #             l_list.append(r_a[1][i])
        #
        #     check list of row numbers in L to see if matches P+1



seq, header = import_genes(file_in)
# L = 'gc$aaac'
L = 'ACTTGGCCCCCCTGTTGATGGAATTTCTGTTTTATGTAACGTAGTATTTA$GAG'
P = 'CAT'
r_a, F, totals = rank_array(L)
fifm(F, L, r_a, P, totals)




# def reverseBwt(bwtSeq):
#     table = [""] * len(bwtSeq)
#     for i in range(len(bwtSeq)):
#         table = sorted(bwtSeq[i] + table[i] for i in range(len(bwtSeq)))
#     s = [row for row in table if row.endswith('$')]
#     ss = "".join(s)
#     s = ss.rstrip('$')
#     print(table)
#     return s