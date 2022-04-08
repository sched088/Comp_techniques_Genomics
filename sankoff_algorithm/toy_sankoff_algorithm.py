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

# score_init = [[0, inf, inf, inf, inf],
#               [inf, 0, inf, inf, inf],
#               [inf, inf, 0, inf, inf],
#               [inf, inf, inf, 0, inf],
#               [inf, inf, inf, inf, inf]]

parsimony_score = 0

rna = ["A", "U", "G", "C", "-"]


# iterate through pairs of leafs

def merge_layers(leafs, score_matrix, tree, parsimony_score, rna):
    tree.append(leafs)
    # print("leafs_l----------------------------" + str(len(leafs)))
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
            # print("leaf: " + str(leaf) + "letter: " + str(letter) + " " + str(pair_2) + str(pair_1))
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
        # print("length new leaves" + str(new_leafs))
        merge_layers(new_leafs, score_matrix, tree, parsimony_score, rna)
    else:
        tree.append(new_leafs)
        print(tree)
        # Build minimum score sequence
    # print(node_scores)
    # for score in node_scores:
    #     print(score)


# num_seq, leafs = import_genes(file_in="aligned_seqs_hw3.txt")
# i = 0
# for leaf in leafs:
#     print(i)
#     print(leaf)
#     print(len(leaf))
#     i += 1
merge_layers(leafs, score_matrix, tree, parsimony_score, rna)
