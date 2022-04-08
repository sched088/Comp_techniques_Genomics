"""
1. (20 points) Implement function bwt that generates the reversible transformation of a sequence:
    - Syntax: bwtSeq = bwt(seq)
    - Input: seq that is a single letter-code representation of a nucleotide sequence.
    - Output: bwtSeq that is BWT of seq.
Use file 'sample1.fa  Download sample1.fa' provided in section Datasets to get the input sequence,
and report the BWT of that sequence to file 'bwtSample1.fa'.
"""

"""
function BWT (string s)
    create a table, where the rows are all possible rotations of s
    sort rows alphabetically
    return (last column of the table)
"""

import numpy as np

file_in = "sample1.fa"


# Read in Sequence
def import_genes(file_in):
    with open(file_in) as f:
        header = f.readline().rstrip()
        # print("Header: " + header)
        seq = ''
        for l, i in enumerate(f):
            seq += i.rstrip()
        # data_1 = {header1: data}
        # print(len(data))
        # print(len(seq))
        # print(seq)
    return seq, header


def bwt(seq, header):
    seq += '$'
    # print(seq)
    # bwm = np.zeros((len(seq), 1), dtype='U1500')  # Burrows-Wheeler Matrix
    bwl = [] #Burrows-Wheeler list (quicker for sorting)

    i = 0
    while i < len(seq):
        bwl.append(seq[i:] + seq[0:i])
        i += 1
    bwls = sorted(bwl) #  Burrows-Wheeler list sorted
    bwtSeq = [] #  Burrows-Wheeler last letters only
    for _ in bwls:
        bwtSeq.append(_[-1])

    header += 'Seq'
    # print("".join(bwtSeq))
    return bwtSeq, header

def bwtfasta(bwtSeq, headerSeq):

    bwtstr = ''
    for _ in bwtSeq:
        bwtstr += _

    with open("bwtSample1.fa", "w") as file_out:
        file_out.write(headerSeq + "\n" + bwtstr)


seq, header = import_genes(file_in)
# seq = "TTGTGTGCATGTTGTTTCATCATTTAGAGATACATTGCGCTGCATCATGGTCA"
bwtSeq, headerSeq = bwt(seq, header)
bwtfasta(bwtSeq, headerSeq)
