import pandas as pd
import gzip
import matplotlib.pyplot as plt
import re
import math
import sys


def exercise1(args):
    _, FASTA, annot, chr, x, y = args
    x = int(x)
    y = int(y)
    data = ''
    print(args)
    # open FASTA file
    f = open(FASTA, 'r')
    print("File found")
    #print(f.readline(5))
    for l, i in enumerate(f):
        print("enumerating...")
        if l != 0:
            data += i
    f.close()
    print("file opened")

    # genome = data.replace('N', '')
    whole_genome = data.replace('\n', '')
    whole_genome = whole_genome.lower()
    genome = whole_genome[x:y]

    print('*********************************************')
    print('Calculate the proportions of each nucleotides')
    print('*********************************************')
    print("rate of a: " + str((genome.count('a') / len(genome)) * 100) + '%')
    print("rate of t: " + str((genome.count('t') / len(genome)) * 100) + '%')
    print("rate of g: " + str((genome.count('g') / len(genome)) * 100) + '%')
    print("rate of c: " + str((genome.count('c') / len(genome)) * 100) + '%')

    summary = []
    p = math.ceil(len(genome) / 100)  # window size
    if p < 5:
        p = 5

    for i in range(p, len(genome), p):
        substring = genome[i - p:i - 1]
        prop_a = genome.count('a') / len(genome)
        prop_t = genome.count('t') / len(genome)
        prop_g = genome.count('g') / len(genome)
        prop_c = genome.count('c') / len(genome)
        sum_ = []

        if substring.count('n') / len(substring) > 0.8:
            # If the ratio of N is over 80 %,
            # just return the proportions of whole region.
            sum_.append(prop_a)
            sum_.append(prop_t)
            sum_.append(prop_g)
            sum_.append(prop_c)
        else:
            sum_ = []
            # i-1 - (i-p) + 1 = p
            sum_.append(substring.count('a') / p)
            sum_.append(substring.count('t') / p)
            sum_.append(substring.count('g') / p)
            sum_.append(substring.count('c') / p)

        summary.append(sum_)
    pd.DataFrame(summary, columns=['a', 't', 'g', 'c']).plot(figsize=[10, 10])
    plt.show()

    t = pd.read_table(
        'hg19_annotation.txt',
        header=None,
        names=[
            'gene names', 'transcript names', 'start',
            'end', 'chr', '5 to 3 or 3 to 5'
        ]
    )
    # print(t)

    t = t[t.chr == "'" + str(chr) + "'"].reset_index(drop=True)
    t = t.drop_duplicates('gene names').reset_index(drop=True)

    pd.DataFrame(summary, columns=['a', 't', 'g', 'c']).plot(figsize=[10, 10], zorder=1)

    for i in t.index:
        if t['start'][i] > x:
            if t['end'][i] < x + len(genome):
                print(t['gene names'][i], t['start'][i], t['end'][i])
                plt.hlines(0.5, (t['start'][i] - x) / p, (t['end'][i] - x) / p, colors='green', lw=4, zorder=10)
    plt.show()


if __name__ == '__main__':
    args = sys.argv
    if len(args) == 6:
        if (args[1][-3:] == '.fa'):
            if args[2][-4:] == '.txt':
                if (args[3][:3] == 'chr') & (args[3][3:].isnumeric()):
                    if (args[4].isnumeric()) & (args[4].isnumeric()):
                        if (int(args[4]) > 0) & (int(args[5]) > 0):
                            exercise1(args)
                        else:
                            print('Fourth and Fifth arguments should set positive number')
                    else:
                        print('Fourth and Fifth arguments should set integer')
                else:
                    print('Third argument should set "chr___')
            else:
                print('Second argument should set ".txt"')
        else:
            print('First argument should set ".fa"')
    else:
        print('The number of arguments is invalid. This file required 5 arguments.')
