"""
1) the array storing BWT
2) the array storing the reference positions to the beginning of the sequence
3) the array storing the rank of the nucleotides.

Report how many mapped reads and their starting position on chrY to text file 'mapping_fullindexFM.txt'.
"""

import numpy as np

# Step 1: Import BWT and create array storing BWT

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


# Step 2: Take BWT and transform into T

def reverseBwt(bwtSeq):
    table = [""] * len(bwtSeq)
    for i in range(len(bwtSeq)):
        table = sorted(bwtSeq[i] + table[i] for i in range(len(bwtSeq)))
    t = [row for row in table if row.endswith('$')]
    st = "".join(s)
    t = st.rstrip('$')
    print(table)
    return t


# Step 3: Array Storing Reference Positions (aka Suffix Array)
# Only need offset positions

def reference_pos(t):
    sa = sorted([(t[i:], i) for i in range(len(t))])
    # print(sa)
    return list(map(lambda x: x[1], sa))

# Step 4: Array Storing Rank of Nucleotides

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
    r_a = np.vstack((r_a_f, r_a_l)).T
    print(dict)
    return r_a, seq_sort, dict

# Step 5:fullindexFA

def fifa(seq, dict, p, sa): # Seq = BWT
    f = 0
    l = len(seq) - 1

    for _ in range(len(p) - 1, -1, -1):  # or reversed(p)?
        if f < 0: ccount = 0
        while f >= 0:
            if seq[f] == p[_]:
                ccount += 1
            f -= 1
        ccount = 0
        while l >= 0:
            if seq[l] == p[_]:
                ccount += 1
            l -= 1
        l = ccount + dict[p[_]] - 1
        f = ccount + dict[p[_]]
        if l < 1: break
        print(dict[p[_]])
        print("f, l: " + str(f) + " " + str(l))

    for r in range(f, l):
        nsteps = 0
        while r not in sa:
            n = seq[r]
            i = r - 1
            if i < 0: return 0
            ccount = 0
            while i >= 0:
                if seq[i] == n:
                    ccount += 1
                i -= 1
            r = ccount + dict(n)
            nsteps += 1
        print("Here: " + str(sa[r] + nsteps))
        return sa[r] + nsteps

t = "TTGTGTGCATGTTGTTTCATCATTTAGAGATACATTGCGCTGCATCATGGTCA"
seq = "ACTTGGCCCCCCTGTTGATGGAATTTCTGTTTTATGTAACGTAGTATTTA$GAG"
p = "CAT"

sa = reference_pos(seq)
r_a, seq_sort, dicty = rank_array(seq)
fifa(seq, dicty, p, sa)

# print(sa)
# for i in sa:
#     print("%2d: %s" % (i, t[i:]))