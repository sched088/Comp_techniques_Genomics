
blosum = {}
lines = []
with open('BLOSUM62.txt', 'r') as f:
    lines = f.readlines()

count = 0
for line in lines:
    if line.startswith('#'):  # strips out lines starting with '#'
        continue
    line_list = line.split()
    if count == 0:  # pull out horizonal list of amino acids for making discitonary.
        reference_aa = line_list
    else:
        score_pairs = {}
        score_count = 0
        for score in line_list:
            if score.isalpha() or score == '*':
                continue
            score_pairs[reference_aa[score_count]] = score

            # print("score_pairs " + str(score_pairs))
            # print(reference_aa[score_count])
            # print("score: " + str(score))
            score_count += 1
            # print("Count: " + str(count))
            # print("l_l: " + str(line_list))
            # print("l_l count: " + str(line_list[count]))
            blosum[reference_aa[count-1]] = score_pairs
    # print(line_list)
    # print(blosum)
    # for letter in line:
    #     print(letter)
    # print(f'line {count}: {line}')
    count += 1

# print(blosum['A']['A'])
# print(blosum['E']['E'])
# print(blosum['R']['R'])
# print(blosum['*']['*'])
