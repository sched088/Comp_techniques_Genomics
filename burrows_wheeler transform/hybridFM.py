"""
1) the array storing BWT
2) the array storing the reference positions to the beginning of the sequence
3) the array storing the rank of the nucleotides.

Report how many mapped reads and their starting position on chrY to text file 'mapping_fullindexFM.txt'.
"""

import numpy as np

file_in = 'bwtChrY(25M-26M).fa'
read1 = 'SRR089545_1.fq'
read2 = 'SRR089545_2.fq'

# import short reads
def import_read1(read1):
    with open(read1) as f:
        # header = f.readline().rstrip()
        # print("Header: " + header)
        reads1 = []
        for l, i in enumerate(f):
            if (l-1) % 4 == 0:

                reads1.append(i.rstrip())

    return reads1 # , header

# def import_read(read2):

# Import BWT text from file

def import_genes(file_in):
    print("here")
    with open(file_in) as f:
        header = f.readline().rstrip()
        seq = ''
        for l, i in enumerate(f):
            seq += i.rstrip()

    return seq, header

# Get T from BWT

def reverseBwt(bwtSeq):
    print("??")
    table = [""] * len(bwtSeq)
    for i in range(len(bwtSeq)):
        table = sorted(bwtSeq[i] + table[i] for i in range(len(bwtSeq)))
    ss = [row for row in table if row.endswith('$')]
    print("reverseBWT" + str(ss))
    s = "".join(ss)
    # s = ss.rstrip('$')
    print(s)
    return s

# Reference positions array / suffix array

def reference_pos(t):
    print("T" + str(t))
    sarray = sorted([(t[i:], i) for i in range(len(t))])
    sarray = list(map(lambda x: x[1], sarray))
    dsa = {}
    for i, s in enumerate(sarray):
        dsa[i] = s
    sa = dsa
    print("SA: " + str(sa))

    """ Hybrid implementation"""

    hybrid_sa = {}
    for v, c in enumerate(sa):
        if c % 3 == 0:
            hybrid_sa[v] = c
    sa = hybrid_sa

    return sa


# Nucleotide rank array
def rank_array(seq):
    print("Seq " + str(seq))
    seq_list_sort = sorted(seq)
    seq_sort = "".join(seq_list_sort)
    r_a_f = []

    # Totals of each character and starting pos. of each
    totals = {}
    for char in seq:
        totals[char] = totals.get(char, 0) + 1
    first_loc = {}
    total_c = 0
    for char, count in sorted(totals.items()):
        first_loc[char] = total_c
        total_c += count

    # Rank for sorted (F)
    rank = 0
    i = 0

    while i < len(seq_sort):
        if seq_sort[i] == seq_sort[i-1]:
            rank += 1
            r_a_f.append(rank)
        else:
            rank = 0
            r_a_f.append(rank)
        i += 1

    # Rank for seq (L)
    dict = {}  # this dict is duplicated oops.
    r_a_l = []
    for _ in seq:
        if _ not in dict:
            dict[_] = 0
        r_a_l.append(dict[_])
        dict[_] += 1
        totals = dict
    # r_a = np.vstack((r_a_f, r_a_l)).T
    r_a = 'placeholder'

    print("SS: " + str(seq_sort))
    return r_a, seq_sort, totals, first_loc

def ranks(L, char, r):
    if r < 0: return 0
    ncc = 0
    while r >= 0:
        if L[r] == char:
            ncc += 1
        r -= 1
    return ncc

def fifm(F, L, r_a, P, totals, first_loc, sa):

    # for index, nucleotide in enumerate(reversed(P)):
    #     for n in F:
    #         if nucleotide == n:
    #             f_list.append(n.index)
    #
    def gotonext(r):
        return ranks(L, L[r], r-1) + first_loc[L[r]]

    start = 0
    stop = len(F) - 1
    for v, n in enumerate(reversed(P)):
        # print("P[v]" + str(n))
        start = ranks(L, n, start-1) + (first_loc[n])
        stop = ranks(L, n, stop) + (first_loc[n]-1)

        if stop < 1: break

    stop += 1
    alignments = []
    for r in range(start, stop):
        steps = 0
        print("r: " + str(r))
        while r not in hybrid_sa:
            r = gotonext(r)
            steps += 1
        alignments.append(sa[r] + steps)
    return alignments

print("1")
L, header = import_genes(file_in)
print("2")


print("3")

T = reverseBwt(L)
print("4")
r_a, F, totals, first_loc = rank_array(L)
print("5")
sa = reference_pos(T)
print("6")
for P in read1:
    print("P")
    alignments = fifm(F, L, r_a, P, totals, first_loc, sa)






# def reverseBwt(bwtSeq):
#     table = [""] * len(bwtSeq)
#     for i in range(len(bwtSeq)):
#         table = sorted(bwtSeq[i] + table[i] for i in range(len(bwtSeq)))
#     s = [row for row in table if row.endswith('$')]
#     ss = "".join(s)
#     s = ss.rstrip('$')
#     print(table)
#     return s