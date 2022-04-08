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
merge = 0


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
    # print(leafs)
    return num_seq, leafs
    # return gene_data


def merge_layers(leafs, score_matrix, tree, parsimony_score, rna, merge):
    merge += 1
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
            parsimony_score = parsimony_score + min(score)
            seq = seq + rna_min
        # print("seq " + str(seq))
        new_leafs.append(seq)
    # print("new_leafs: " + str(new_leafs))

    print("Parsimony score at M" + str(merge) + " = " + str(parsimony_score))
    print("Most likely sequence at M" + str(merge) + " = " + str(new_leafs))
    print("")

    return parsimony_score, new_leafs, merge
    # if len(new_leafs) >= 2:
    #     print("length new leaves" + str(new_leafs))
    #     merge_layers(new_leafs, score_matrix, tree, parsimony_score, rna)
    # else:
    #     tree.append(new_leafs)
    #     print("Tree " + str(tree))
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

m1 = [leafs[leafs.index('>silva|U90316|1|1541| Salmonella typhimurium') + 1],
      leafs[leafs.index('>silva|EF643609|1|1542| Shigella flexneri') + 1]]
parsimony_score, m1, merge = merge_layers(m1, score_matrix, tree, parsimony_score, rna, merge)

m2 = [m1[0], leafs[leafs.index('>silva|AB045730|1|1534| Escherichia coli') + 1]]
parsimony_score, m2, merge = merge_layers(m2, score_matrix, tree, parsimony_score, rna, merge)

m3 = [leafs[leafs.index('>silva|D78004|1|1483| Photorhabdus luminescens') + 1],
      leafs[leafs.index('>silva|AJ232231|1|1477| Yersinia pestis') + 1]]
parsimony_score, m3, merge = merge_layers(m3, score_matrix, tree, parsimony_score, rna, merge)

m4 = [m2[0], m3[0]]
parsimony_score, m4, merge = merge_layers(m4, score_matrix, tree, parsimony_score, rna, merge)

m5 = [leafs[leafs.index('>silva|BA000021|136582|138132| Wigglesworthia glossinidia endosymbiont of '
                        'Glossina brevipalpis') + 1],
      leafs[leafs.index('>silva|M27039|1|1547| Buchnera aphidicola') + 1]]
parsimony_score, m5, merge = merge_layers(m5, score_matrix, tree, parsimony_score, rna, merge)

m6 = [m4[0], m5[0]]
parsimony_score, m6, merge = merge_layers(m6, score_matrix, tree, parsimony_score, rna, merge)

m7 = [leafs[leafs.index('>silva|AY513500|1|1491| Vibrio cholerae') + 1],
      leafs[leafs.index('>silva|CR378665|83603|85179| Photobacterium profundum') + 1]]
parsimony_score, m7, merge = merge_layers(m7, score_matrix, tree, parsimony_score, rna, merge)

m8 = [m6[0], m7[0]]
parsimony_score, m8, merge = merge_layers(m8, score_matrix, tree, parsimony_score, rna, merge)

m9 = [m8[0], leafs[leafs.index('>silva|AY613472|1|1538| Haemophilus influenzae') + 1]]
parsimony_score, m9, merge = merge_layers(m9, score_matrix, tree, parsimony_score, rna, merge)

m10 = [leafs[leafs.index('>silva|EU918692|1|1499| Pasteurella multocida') + 1],
       leafs[leafs.index('>silva|AE014299|2973049|2974577| Shewanella oneidensis MR-1') + 1]]
parsimony_score, m10, merge = merge_layers(m10, score_matrix, tree, parsimony_score, rna, merge)

m11 = [m10[0], m9[0]]
parsimony_score, m11, merge = merge_layers(m11, score_matrix, tree, parsimony_score, rna, merge)

m12 = [leafs[leafs.index('>silva|EU620069|1|1570| Pseudomonas putida') + 1],
       leafs[leafs.index('>silva|AY574913|1|1531| Pseudomonas syringae') + 1]]
parsimony_score, m12, merge = merge_layers(m12, score_matrix, tree, parsimony_score, rna, merge)

m13 = [m12[0], leafs[leafs.index('>silva|X06684|1|1537| Pseudomonas aeruginosa') + 1]]
parsimony_score, m13, merge = merge_layers(m13, score_matrix, tree, parsimony_score, rna, merge)

m14 = [m13[0], leafs[leafs.index('>silva|EU071537|1|1509| uncultured Acinetobacter sp.') + 1]]
parsimony_score, m14, merge = merge_layers(m14, score_matrix, tree, parsimony_score, rna, merge)

m15 = [m14[0], m11[0]]
parsimony_score, m15, merge = merge_layers(m15, score_matrix, tree, parsimony_score, rna, merge)

m16 = [m15[0], leafs[leafs.index('>silva|AY342037|1|1457| Coxiella burnetii') + 1]]
parsimony_score, m16, merge = merge_layers(m16, score_matrix, tree, parsimony_score, rna, merge)

m17 = [leafs[leafs.index('>silva|BX897699|1583137|1584628| Bartonella henselae') + 1],
       leafs[leafs.index('>silva|AJ250247|1|1455| Bartonella quintana') + 1]]
parsimony_score, m17, merge = merge_layers(m17, score_matrix, tree, parsimony_score, rna, merge)

m18 = [m17[0], leafs[leafs.index('>silva|AY513527|1|1412| Brucella melitensis') + 1]]
parsimony_score, m18, merge = merge_layers(m18, score_matrix, tree, parsimony_score, rna, merge)

m19 = [m18[0], leafs[leafs.index('>silva|AF208508|1|1472| Bradyrhizobium japonicum') + 1]]
parsimony_score, m19, merge = merge_layers(m19, score_matrix, tree, parsimony_score, rna, merge)

m20 = [m19[0], leafs[leafs.index('>silva|AE005673|2840172|2841610| Caulobacter crescentus CB15') + 1]]
parsimony_score, m20, merge = merge_layers(m20, score_matrix, tree, parsimony_score, rna, merge)

m21 = [leafs[leafs.index('>silva|AF084850|1|1436| Bdellovibrio bacteriovorus') + 1],
       leafs[leafs.index('>silva|AJ235272|177323|178829| Rickettsia prowazekii') + 1]]
parsimony_score, m21, merge = merge_layers(m21, score_matrix, tree, parsimony_score, rna, merge)

m22 = [m21[0], leafs[leafs.index('>silva|AY026912|1|1458| Wolbachia pipientis') + 1]]
parsimony_score, m22, merge = merge_layers(m22, score_matrix, tree, parsimony_score, rna, merge)

m23 = [m22[0], m20[0]]
parsimony_score, m23, merge = merge_layers(m23, score_matrix, tree, parsimony_score, rna, merge)

m24 = [leafs[leafs.index('>silva|U07573|1|1470| Helicobacter hepaticus') + 1],
       leafs[leafs.index('>silva|AY062899|1|1492| Helicobacter pylori') + 1]]
parsimony_score, m24, merge = merge_layers(m24, score_matrix, tree, parsimony_score, rna, merge)

m25 = [m24[0], leafs[leafs.index('>silva|AANK01000003|466|1969| Campylobacter jejuni subsp. jejuni 260.94') + 1]]
parsimony_score, m25, merge = merge_layers(m25, score_matrix, tree, parsimony_score, rna, merge)

m26 = [leafs[leafs.index('>silva|AY362360|1|1510| Desulfovibrio vulgaris') + 1],
       leafs[leafs.index('>silva|AE017180|1224149|1225696| Geobacter sulfurreducens PCA') + 1]]
parsimony_score, m26, merge = merge_layers(m26, score_matrix, tree, parsimony_score, rna, merge)

m27 = [m26[0], leafs[leafs.index('>silva|U12460|1|1436| Rickettsia conorii') + 1]]
parsimony_score, m27, merge = merge_layers(m27, score_matrix, tree, parsimony_score, rna, merge)

m28 = [m25[0], m27[0]]
parsimony_score, m28, merge = merge_layers(m28, score_matrix, tree, parsimony_score, rna, merge)

m35 = [m28[0], m23[0]]
parsimony_score, m35, merge = merge_layers(m35, score_matrix, tree, parsimony_score, rna, merge)

m29 = [leafs[leafs.index('>silva|BX640428|97283|98814| Bordetella parapertussis') + 1],
       leafs[leafs.index('>silva|BX640417|66117|67648| Bordetella pertussis') + 1]]
parsimony_score, m29, merge = merge_layers(m29, score_matrix, tree, parsimony_score, rna, merge)

m30 = [m29[0], leafs[leafs.index('>silva|AL646052|3277704|3279237| Ralstonia solanacearum') + 1]]
parsimony_score, m30, merge = merge_layers(m30, score_matrix, tree, parsimony_score, rna, merge)

m31 = [leafs[leafs.index('>silva|AF310340|1|1471| Neisseria meningitidis') + 1],
       leafs[leafs.index('>silva|EU693450|1|1505| Chromobacterium violaceum') + 1]]
parsimony_score, m31, merge = merge_layers(m31, score_matrix, tree, parsimony_score, rna, merge)

m32 = [m30[0], m31[0]]
parsimony_score, m32, merge = merge_layers(m32, score_matrix, tree, parsimony_score, rna, merge)

m33 = [leafs[leafs.index('>silva|AF442743|1|2165| Xanthomonas axonopodis pv. citri') + 1],
       leafs[leafs.index('>silva|AF310340|1|1471| Neisseria meningitidis') + 1]]
parsimony_score, m33, merge = merge_layers(m33, score_matrix, tree, parsimony_score, rna, merge)

m34 = [m33[0], m32[0]]
parsimony_score, m43, merge = merge_layers(m34, score_matrix, tree, parsimony_score, rna, merge)

m36 = [m35[0], m34[0]]
parsimony_score, m36, merge = merge_layers(m36, score_matrix, tree, parsimony_score, rna, merge)

m37 = [m36[0], m16[0]]
parsimony_score, m37, merge = merge_layers(m37, score_matrix, tree, parsimony_score, rna, merge)
