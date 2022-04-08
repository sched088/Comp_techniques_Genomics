"""
Utilize Viterbi algorithm

Initialize start probabilities a01 + ... + a0k = 1
(equal probabilities of starting with any emission)

"""

import numpy as np
import random
import sys

args = sys.argv

e = np.array([
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 1],
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 1]])

p = 0.98

q = 0.999

# hardcoded in from CSV document
a = np.array([
    [0.180 * p, 0.274 * p, 0.426 * p, 0.120 * p, (1 - p) / 4, (1 - p) / 4, (1 - p) / 4, (1 - p) / 4],
    [0.171 * p, 0.368 * p, 0.274 * p, 0.188 * p, (1 - p) / 4, (1 - p) / 4, (1 - p) / 4, (1 - p) / 4],
    [0.161 * p, 0.339 * p, 0.375 * p, 0.125 * p, (1 - p) / 4, (1 - p) / 4, (1 - p) / 4, (1 - p) / 4],
    [0.079 * p, 0.355 * p, 0.384 * p, 0.182 * p, (1 - p) / 4, (1 - p) / 4, (1 - p) / 4, (1 - p) / 4],
    [(1 - q) / 4, (1 - q) / 4, (1 - q) / 4, (1 - q) / 4, 0.300 * q, 0.205 * q, 0.285 * q, 0.210 * q],
    [(1 - q) / 4, (1 - q) / 4, (1 - q) / 4, (1 - q) / 4, 0.322 * q, 0.298 * q, 0.078 * q, 0.302 * q],
    [(1 - q) / 4, (1 - q) / 4, (1 - q) / 4, (1 - q) / 4, 0.248 * q, 0.246 * q, 0.298 * q, 0.208 * q],
    [(1 - q) / 4, (1 - q) / 4, (1 - q) / 4, (1 - q) / 4, 0.177 * q, 0.239 * q, 0.292 * q, 0.292 * q]])

islands = 0
coding_regions = 0
seq_start = 43507093  # Assignment gives in 1-indexing, so subtract 1 for python indexing.
output_text = []
island_ends = []


def HMM(a, e, seq):
    print("Running HMM")
    start_b = random.randint(0, 3)  # equal distribution start
    emission_prob = e[:, start_b] * a[start_b, :]
    pi = []
    regions = []
    # emissions = (['A+', 'C+', 'G+', 'T+', 'A-', 'C-', 'G-', 'T-'])
    # only 4 letters in alphabet, ATCG -> then possible emissions are 2 for each one of those (A+ or A-)
    # V0 = 1
    for x in seq[1:]:
        if x == 'N': # for real fasta file. Would this throw off nucleotide positioning?
            continue
        val = ['A', 'C', 'G', 'T']

        Vi = a * np.row_stack(emission_prob)
        Vi = abs(np.ma.log(Vi))  # log to prevent approaching 0
        max_Vi = np.argmin(Vi, axis=0) # min because log to get max from 0
        emission_prob = e[:, val.index(x)] * Vi[max_Vi, np.linspace(0, a.shape[0] - 1, num=a.shape[0]).astype(int)]
        # linspace to create evenly spaced matrix values
        pi.append(max_Vi)
    ptrs = [np.argmax(emission_prob)]
    if np.argmax(emission_prob) >= 4:
        regions = ['-']
    else:
        regions =['+']
    while pi:
        max_Vi = pi.pop()
        ptrs.append(max_Vi[ptrs[-1]])
        if max_Vi[ptrs[-1]] >= 4:
            regions.append('-')
        else:
            regions.append('+')
    ptrs.reverse() # for scoring values not used in output
    regions.reverse() # for CpG y/n
    return ptrs, regions


def import_sequence(args):
    fasta = args
    with open(fasta) as f:
    # with open('seq1.fa') as f:
        header1 = f.readline().rstrip()
        data = ''
        for l, i in enumerate(f):
            data += i.rstrip()
        print(len(data))
    return data


def import_genes():
    with open("Chr21.txt") as f:
        # header1 = f.readline().rstrip()
        # print("Header: " + header)
        gene_data = []
        for l, i in enumerate(f):
            clean_text = i.rstrip()
            clean_text = clean_text.split('\t')
            gene_data.append(clean_text)

        for gene in gene_data:  # convert for search later
            gene[0] = int(gene[0])
            gene[1] = int(gene[1])
    return gene_data


def GeneIdentifier(seq, genes, islands, coding_regions):
    print("Running GeneIdentifier")
    island_starts = []
    island_lengths = []
    island_start = 0

    for i in range(len(seq)):
        if seq[i] == '+' and i == 0:
            islands += 1
            in_island = True
            island_start = i + seq_start
            island_starts.append(island_start)
        if seq[i] == '+' and seq[i - 1] == '-':
            islands += 1
            in_island = True
            island_start = i + seq_start
            island_starts.append(island_start)

        if seq[i] == '-' and seq[i - 1] == '+':
            island_end = i + seq_start
            island_ends.append(island_end)
            island_length = island_end - island_start
            island_lengths.append(island_length)
            in_island = False
            if island_length > 200:
                gene_ids = []
                for gene in genes:
                    if 500 >= (gene[0] - island_end) >= 0:  # check forward strand
                        gene_ids.append([gene[3], gene[2]])
                    if 500 >= (island_start - gene[1]) >= 0:  # check backward strand
                        gene_ids.append([gene[3], gene[2]])
                if len(gene_ids) == 0:
                    gene_ids.append("no_gene")
                else:
                    coding_regions += 1
                island_details = ["CpG Island " + str(islands) + ": " + str(island_length) + "bp (" + str(island_start)
                                  + "-" + str(island_end) + ") " + str(gene_ids)]
                output_text.append(island_details)

        if i == len(seq) - 1 and seq[i] == '+':  # condition ending sequence
            island_end = i + seq_start
            island_ends.append(island_end)
            island_length = island_end - island_start
            in_island = False
            if island_length > 200:
                gene_ids = []
                for gene in genes:
                    if 500 >= (gene[0] - island_end) >= 0:  # check forward strand
                        gene_ids.append([gene[3], gene[2]])
                    if 500 >= (island_start - gene[1]) >= 0:  # check backward strand
                        gene_ids.append([gene[3], gene[2]])
                if len(gene_ids) == 0:
                    gene_ids.append("no_gene")
                else:
                    coding_regions += 1
                island_details = ["CpG Island " + str(islands) + ": " + str(island_length) + "bp (" + str(island_start)
                                  + "-" + str(island_end) + ") " + str(gene_ids)]
                output_text.append(island_details)
    return islands, coding_regions, output_text


def save_file(islands, coding_regions, output_text):
    print(islands)
    print(coding_regions)
    print(output_text)
    print("Exporting...")
    file = open("Homework 2 Export", "w")
    file.write(f'Total CpG islands found: {islands}; {coding_regions} out of {islands} islands are followed by a '
               f'coding region \n')
    for output in output_text:
        file.write(str(output) + '\n')
    print("Done!")


if __name__ == '__main__':
    args = sys.argv
    match = False
    print("Your input arguments are: " + str(args))
    if len(args) != 2:
        args.append("no match")
    if args[1][-3:] != '.py':
        print('First argument should set "program_name.py"')
    if args[2][-6:] != '.fasta':
        print('Second argument should set to first FASTA file ".fasta"')
        match = True
    print(args[2])
    seq = import_sequence(args[2])
    genes = import_genes()
    scores, regions = HMM(a, e, seq)
    islands, coding_regions, output_text = GeneIdentifier(regions, genes, islands, coding_regions)
    save_file(islands, coding_regions, output_text)
    print("Total CpG islands found: " + str(islands) + "; " + str(coding_regions) + " out of " + str(islands) +
          " islands are followed by a coding region")
    for output in output_text:
        print(output)

# print(output_text)
# print(island_count)
# print(coding_regions)
