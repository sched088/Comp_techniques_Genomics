"""
2. (20 points) Implement function reverseBwt that reverses the Burrows-Wheeler Transform:
    - Syntax: seq = reverseBwt(bwtSeq)
    - Input: bwtSeq that is the BWT of a sequence seq.
    - Output: seq that is the original nucleotide sequence such that BWT(seq) = bwtSeq.
Use file 'bwtHomoSapiens.fa  Download bwtHomoSapiens.fa' provided in section Datasets to get the input sequence,
and report the reversed BWT of that sequence to file 'rBwtSample2.fa'.
"""

import numpy as np

file_in = 'bwtHomoSapiens.fa'


# Read in Sequence
def import_genes(file_in):
    print("here")
    with open(file_in) as f:
        header = f.readline().rstrip()
        header += 'revBwtSeq'
        # print("Header: " + header)
        seq = ''
        for l, i in enumerate(f):
            seq += i.rstrip()
        # data_1 = {header1: data}
        # print(len(data))
        print(len(seq))
        # print(seq)
    return seq, header

"""
This first attempt is too computationally heavy and utilizes inverse transformation.
"""
# def reverseBwt(bwtSeq):
#     seq = []
#
#     for _ in bwtSeq:
#         seq.append(_)
#     # print(seq[:6])
#     seq = sorted(seq)
#     # print(seq[:6])
#
#     i = 0
#     while i < len(bwtSeq):
#         j = 0
#         while j < len(bwtSeq):
#             seq[j] = bwtSeq[j] + seq[j]  # Add
#             j += 1
#             # print(i)
#         # print(seq[:6])
#         seq = sorted(seq)  # Sort
#         # print(seq[:6])
#         i += 1
#         if i % 1000 == 0: print(str(i/329518) + " percent complete")
#
#     for _ in seq:
#         if _[-1] == '$':
#             rbwt = _
#             print(rbwt)
#             return rbwt

"""
Second implementation 45 min -- too slow. Two while loops.
"""
# def reverseBwt(bwtSeq):
#     # print(bwtSeq)
#     bwts = sorted(bwtSeq)
#     bwtss = "".join(bwts)
#     # print(bwtss)
#     seq = ''
#     seq += bwtSeq[0]  # Set starting char to avoid recurring if statement in while loop
#     find = seq[0]
#     loc = 0
#
#     i = 1
#     while i < len(bwtSeq):
#         # print("find: " + find)
#         num = bwtSeq[:loc+1].count(find)
#         # print(num)
#         bwt_index = find_nth(bwtss, find, num)
#         seq = bwtSeq[bwt_index] + seq
#         # print(seq)
#         find = seq[0]
#         loc = bwt_index
#         i += 1
#         if i % 1000 == 0: print(str(i/329518) + " percent complete")
#     seq = seq[1:]
#     print(seq)
#
# def find_nth(haystack, needle, n):
#     next = haystack.find(needle)
#     while next >= 0 and n > 1:
#         next = haystack.find(needle, next + len(needle))
#         n -= 1
#     return next

# Wiki implementation for speed.

def reverseBwt(bwtSeq):
    table = [""] * len(bwtSeq)
    for i in range(len(bwtSeq)):
        if i % 1000 == 0: print(str(i/329518) + " percent complete")
        table = sorted(bwtSeq[i] + table[i] for i in range(len(bwtSeq)))
    s = [row for row in table if row.endswith('$')]
    ss = "".join(s)
    s = ss.rstrip('$')
    print(table)
    return s


def bwtfasta(revbwtSeq, headerSeq):

    revbwtstr = ''
    for _ in revbwtSeq:
        revbwtstr += _

    with open("rBwtSample2.fa", "w") as file_out:
        file_out.write(headerSeq + "\n" + revbwtstr)

bwtSeq, header = import_genes(file_in)
# bwtSeq = 'gc$aaac'
revseq = reverseBwt(bwtSeq)
print(revseq)
