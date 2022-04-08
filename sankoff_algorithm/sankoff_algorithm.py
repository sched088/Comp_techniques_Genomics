"""
Implement the Sankoff algorithm and test your algorithm with a few toy examples \
to demonstrate that your implementation is correct.

* As a toy example, if you got
                ['AUUCGUGAUU', 'AUUGAA-AUU', 'GUCCUCGGUU', 'GA-CACGAUC'],
the model should be returned
                ['AUUCGUGAUU', 'AUUGAA-AUU', 'GUCCUCGGUU', 'GA-CACGAUC', 'AUUCAUGAUU', 'GUUCACGAUU', 'AUUCAUGAUU'].
"""

"""
1) Iterate through each letter of each leaf sequence
    a) each letter is compared to neighboring branch's letter
    b) score is assigned by looking at score_matrix
    c) at each node, minimum score is taken and that letter is assigned. 

"""

# import numpy as np
import math

inf = math.inf

leafs = ['AUUCGUGAUU', 'AUUGAA-AUU', 'GUCCUCGGUU', 'GA-CACGAUC']
tree = []

score_matrix = [[0, 3, 4, 9, 8],
                [3, 0, 2, 4, 8],
                [4, 2, 0, 4, 8],
                [9, 4, 4, 0, 8],
                [8, 8, 8, 8, 0]]  # should corner 8 be 0?

score_init = [[0, inf, inf, inf, inf],
              [inf, 0, inf, inf, inf],
              [inf, inf, 0, inf, inf],
              [inf, inf, inf, 0, inf],
              [inf, inf, inf, inf, inf]]

parsimony_score = 0

rna = ["A", "U", "G", "C", "-"]


# iterate through pairs of leafs
def import_genes(file_in):
    with open(file_in) as f:
        seq = ''
        leafs = []
        num_seq = 0
        # header1 = f.readline().rstrip()
        # print("Header: " + header1)
        gene_data = ''
        for l, i in enumerate(f):
            if l == 0:
                num_seq += 1
                seq = str(i.rstrip())
                leafs.append(seq)
                seq = ''
                continue
            if i[0] == ">":
                num_seq += 1
                leafs.append(seq)
                seq = str(i.rstrip())
                leafs.append(seq)
                seq = ''
                continue
            seq = seq + str(i.rstrip())
            # gene_data += i.rstrip()
            # data_1 = {header1: data}
            # print(len(data))
            # print(len(gene_data))
        leafs.append(seq)
    print("There are " + str(num_seq) + " sequences being compared")
    print(leafs)
    return num_seq, leafs
    # return gene_data


def merge_layers(leafs, score_matrix, tree, parsimony_score, rna):
    tree.append(leafs)
    print("leafs_l----------------------------" + str(len(leafs)))
    new_leafs = []
    # print(leaf_length)
    for leaf in range(0, len(leafs), 2):
        leaf_length = len(leafs[leaf])
        # print(leafs[leaf + 1])
        node_scores = []
        # iterate through each letter in pairs of leafs
        for letter in range(leaf_length):
            # print(leaf_length)
            # print("letter: " + str(letter))
            # print("leaf: " + str(leaf))
            node_score = [0, 0, 0, 0, 0]
            pair_1 = leafs[leaf][letter]
            pair_2 = leafs[leaf + 1][letter]
            print("leaf: " + str(leaf) + "letter: " + str(letter) + " " + str(pair_2) + str(pair_1))
            if pair_1 not in rna:
                pair_1 = "-"
            if pair_2 not in rna:
                pair_2 = "-"
            # print(str(pair_1) + " " + str(pair_2))

            # iterate through RNA options for each letter pair
            for n in rna:
                # print(n)
                # print(node_score[rna.index(n)])
                # print(score_matrix[rna.index(n)][rna.index(pair_1)])
                # print(score_matrix[rna.index(n)][rna.index(pair_1)])
                node_score[rna.index(n)] = score_matrix[rna.index(n)][rna.index(pair_1)] + \
                                           score_matrix[rna.index(n)][rna.index(pair_2)]
            # print("Node score: " + str(node_score))
            node_scores.append(node_score)
            # print(node_scores)

        seq = ''
        for score in node_scores:
            # print("score" + str(score))
            # print(score.index(min(score)))
            # print(rna[score.index(min(score))])
            rna_min = rna[score.index(min(score))]
            if rna_min == 0:
                parsimony_score += 1
            seq = seq + rna_min
            # print(seq)
        new_leafs.append(seq)
    # print("new_leafs: " + str(new_leafs))

    if len(new_leafs) >= 2:
        print("length new leaves" + str(new_leafs))
        merge_layers(new_leafs, score_matrix, tree, parsimony_score, rna)
    else:
        tree.append(new_leafs)
        print(tree)
        # Build minimum score sequence
    # print(node_scores)
    # for score in node_scores:
    #     print(score)


num_seq, leafs = import_genes(file_in="aligned_seqs_hw3.txt")
# i = 0
# for leaf in leafs:
#     print(i)
#     print(leaf)
#     print(len(leaf))
#     i += 1
# merge_layers(leafs, score_matrix, tree, parsimony_score, rna)
#
# m1 = [>silva|U90316|1|1541| Salmonella typhimurium, >silva|EF643609|1|1542| Shigella flexneri]
# m2 = [m1 , >silva|AB045730|1|1534| Escherichia coli]
# m3 = [>silva|D78004|1|1483| Photorhabdus luminescens, >silva|AJ232231|1|1477| Yersinia pestis]
# m4 = [m2, m3]
# m5 = [>silva|BA000021|136582|138132| Wigglesworthia glossinidia endosymbiont of Glossina, >silva|M27039|1|1547| Buchnera aphidicola]
# m6 = [m4, m5]
# m7 = [>silva|AY513500|1|1491| Vibrio cholerae, >silva|CR378665|83603|85179| Photobacterium profundum]
# m8 = [m6, m7]
# m9 = [m8, >silva|AY613472|1|1538| Haemophilus influenzae]
# m10 = [>silva|EU918692|1|1499| Pasteurella multocida, >silva|AE014299|2973049|2974577| Shewanella oneidensis MR-1]
# m11 = [m10, m9]
# m12 = [>silva|EU620069|1|1570| Pseudomonas putida, >silva|AY574913|1|1531| Pseudomonas syringae]
# m13 = [m12, >silva|X06684|1|1537| Pseudomonas aeruginosa]
# m14 = [m13, >silva|EU071537|1|1509| uncultured Acinetobacter sp.]
# m15 = [m14, m11]
# m16 = [m15, >silva|AY342037|1|1457| Coxiella burnetii]
