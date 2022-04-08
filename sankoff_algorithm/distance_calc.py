import Bio
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import *
import numpy

score_matrix = [[0, 3, 4, 9, 8],
                [3, 0, 2, 4, 8],
                [4, 2, 0, 4, 8],
                [9, 4, 4, 0, 8],
                [8, 8, 8, 8, 0]]

numpy.array(score_matrix)

alignment = Bio.AlignIO.read(open("aligned_seqs_hw3.fasta"), "fasta")
calc = DistanceCalculator('identity')
distance_matrix = calc.get_distance(alignment)

constructor = DistanceTreeConstructor(calc, 'nj')
tree = constructor.build_tree(alignment)
p_score = ParsimonyScorer(matrix=score_matrix)
print(p_score)
# print(distance_matrix)
print(tree)
with open('HW3 distance_matrix.txt', 'w') as f:
    f.write(str(distance_matrix))