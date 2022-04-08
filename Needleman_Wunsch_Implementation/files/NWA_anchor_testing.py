"""
This is the first homework for class 5481 Computational Genomics
Max Scheder

For anchored we can
1) Align up to point of first anchor
2) Align anchor region in case not perfect match
-repeat for anchors
3) Align until end

human_start, human_end : fly_start, fly_end
"""

import sys
import numpy

np = numpy
args = sys.argv

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
            score_count += 1
            blosum[reference_aa[count - 1]] = score_pairs
    count += 1


def _import(args):  # this code is confusing to me
    _, fasta1, fasta2, matches = args
    data_1 = ''
    data_2 = ''
    data_match = ''
    data = ''
    header = ''


''' Read in file data as dictionary '''
# Need to replace specific file name references with fasta1, fasta2 from args
# Need to strip '>' from headers

with open("Human_HOX.fa") as f:
    header1 = f.readline().rstrip()
    # print("Header: " + header)
    data = ''
    for l, i in enumerate(f):
        data += i.rstrip()
    data_1 = {header1: data}
# print(data_1['>fly_HOX'])

with open("fly_HOX.fa") as f:
    header2 = f.readline().rstrip()
    data = ''
    for l, i in enumerate(f):
        data += i.rstrip()
    data_2 = {header2: data}
# print(data_2['>human_HOX'])
# print(data_2['>human_HOX'][-1])
# print(str(data_1) + '\n' + str(data_2))

f_matches = []
with open("Match_HOX.txt") as f:
    matches = f.readlines()
    for match in matches:
        match_set = match.split()
        f_matches.append(match_set)
matches = f_matches

"""TEST DATA"""

# data_1 = "GCATATT"
# data_2 = "XGCATTA"
# print(data_1)
# print(data_2)

# ld1 = int(len(data_1['>fly_HOX']))
# ld2 = int(len(data_2['>human_HOX']))
ld1 = len(data_1[header1])
ld2 = len(data_2[header2])
print(ld1, ld2)


# for line in data:
#     matches.append(data)
# print("Matches" + str(matches))

# def blosum_lookup():
#     """We want to look up the column and row associated with i-1 str and j-1 str """
#     blosum = np.loadtxt('BLOSUM62.txt', dtype='str', comments='#')
#     print(blosum)
#     # return score
#
#
# blosum_lookup()


def _init_(data_1, data_2, d=-5):
    # print(ld1, ld2)
    """ Initialization:
        F(0,0) = 0
        F(0,j) = -j x d
        F(i,0) = -i x d
    """
    F = np.zeros((ld1 + 1, ld2 + 1))  # create array/matrix of sequences being compared
    F[:, 0] = np.linspace(0, -ld1, ld1 + 1)  # alternately could do loop through rows/columns and add penalty * i
    F[0, :] = np.linspace(0, -ld2, ld2 + 1)  # alternately could do loop through rows/columns and add penalty * i
    F *= 5
    # print(F)

    """ Pointers:
                    { Diag = (1,1)  if [case 1]
        ptr(i,j) =  { Left = (1,0)  if [case 2]
                    { Up = (0,1)    if [case 3]
    """
    ptr = F.astype(str)
    ptr[:, 0] = 'H'
    ptr[0, :] = 'V'
    ptr[0, 0] = 0
    # print(ptr)

    nwa_anchor(F, ptr, d)
    """ Induction:
                   { F(i-1, j-1) + s(xi, yj)    [case 1]
     F(i,j) = max  { F(i-1, j) - d              [case 2]
                   { F(i, j-1) - d              [case 3]

        This section populates two matrices, one with the values, the other with the pointer direction.
        Pointer direction not needed for assignment, but helped me visualize.

    For the implementation of anchoring we are taking the following strategy:
    1) Align up to point of first anchor
    2) Align anchor region in case not perfect match
    -repeat for anchors
    3) Align from end of last anchor until end

    """


def nwa_anchor(F, ptr, d):
    print("Enter Anchor")
    data_1_offset = 0
    data_2_offset = 0
    c_data_1 = []
    c_data_2 = []
    for _ in range(len(matches)):
        # print("Anchor #: " + str(_))
        # print("Matches1 Start: " + str(matches[_][0]))
        # print("Matches1 End: " + str(matches[_][1]))
        # print("Data_1 Offset: " + str(data_1_offset))
        # print("Data_2 Offset: " + str(data_2_offset))
        # pre = nwa63(F, ptr, data_1[data_1_offset:(int(matches[_][0])-1)],
        #            data_2[data_2_offset:(int(matches[_][2])-1)])
        # anchor = nwa63(F, ptr, data_1[int(matches[_][0]):int(matches[_][1])],
        #               data_2[int(matches[_][0]):(int(matches[_][1]))])
        v_seg = "pre"
        pre = nwa63(F, ptr, d, data_1_offset, int(matches[_][0]) - 1, data_2_offset, int(matches[_][2]) - 1, v_seg)
        # print("TEST ANCHOR: " + str(int(matches[_][0])))
        # print("END ANCHOR: ")
        # print(int(matches[_][1]), int(matches[_][3]))
        v_seg = "anchor"
        anchor = nwa63(F, ptr, d, int(matches[_][0])-1, int(matches[_][1]),
                       int(matches[_][2])-1, int(matches[_][3]), v_seg)
        # print("TEST POST: " + str(int(matches[len(matches) - 1][1])))
        # print("LD1 , LDS2: ")
        # print(ld1, ld2)
        # c_data_1 += anchor[0]
        # c_data_2 += anchor[1]
        c_data_1 += pre[0] + anchor[0]
        c_data_2 += pre[1] + anchor[1]
        data_1_offset = int(matches[_][1])
        data_2_offset = int(matches[_][3])
        print(_)
        print("pre: " + str(pre))
        print("anchor: " + str(anchor))
        print("Compiled Data_1: " + str(c_data_1))
        print("Compiled Data_2: " + str(c_data_2))
        print("Data_1 Offset: " + str(data_1_offset))
        print("Data_2 Offset: " + str(data_2_offset))
    #
    # print(int(matches[len(matches) - 1][1]))
    v_seg = "post"
    post = nwa63(F, ptr, d, int(matches[len(matches) - 1][1]), ld1, int(matches[len(matches) - 1][3]), ld2, v_seg)
    c_data_1 += post[0]
    c_data_2 += post[1]
    # print("post: " + str(post))

    print(*c_data_1, sep='')
    print(*c_data_2, sep='')
    return c_data_1, c_data_2


def nwa63(F, ptr, d, start_1, stop_1, start_2, stop_2, v_seg):
    print("NWA63")
    print(start_1 + 1, stop_1 + 1)
    print(start_2 + 1, stop_2 + 1)
    for i in range(start_1 + 1, stop_1 + 1):
        # print("populate i: " + str(i))
        for j in range(start_2 + 1, stop_2 + 1):
            # print("populate j: " + str(j))
            # print(stop_2)
            # print('hit')
            # print("i" + str(i))
            # print("j" + str(j))
            # print(data_1['>fly_HOX'][i-1])
            # print(data_2['>human_HOX'][j-1])
            # if data_1[i-1] == data_2[j-1]:
            # if data_1['>fly_HOX'][i - 1] == data_2['>human_HOX'][j - 1]:
            # print("1")
            diag = F[i - 1][j - 1] + int(blosum[data_1[header1][i - 1]][data_2[header2][j - 1]])
            # print("2")
            # print(F[i-1][j-1])
            # print("DATA 1: " + str(data_1[header1][i-1]))
            # print("DATA 2: " + str(data_2[header2][j-1]))
            # print(F)
            # print(F[i - 1][j])
            # print(F[i][j - 1])
            # print(d)
            # print(int(blosum[data_1[i-1]][data_2[j-1]]))
            # else:
            # diag = F[i - 1][j - 1] + 1
            hori = F[i - 1][j] + d
            vert = F[i][j - 1] + d
            # print(diag, vert, hori)
            F[i][j] = max(diag, vert, hori)
            # print("3")
            if diag == max(diag, vert, hori):
                ptr[i][j] = 'D'
                # print("DIAG")
            elif vert == max(diag, vert, hori):
                ptr[i][j] = 'V'
                # print("VERT")
            else:
                ptr[i][j] = 'H'
                # print("HORI")

    # print(F)
    # print(ptr)
    # print("i, j")
    # print(i, j)
    """ backwards trace through optimal alignment"""
    print("BACKTRACE")
    n = stop_1
    m = stop_2
    # if v_seg == "pre" or v_seg == "post":
    #     n = stop_1
    #     m = stop_2
    # else:
    #     n = stop_1 - 1
    #     m = stop_2 - 1
    print("n/stop_1, m/stop_2")
    print(n, m)
    print("start_1, start_2")
    print(start_1, start_2)
    r_data_1 = []
    r_data_2 = []

    # with open("Output_F.txt", "w") as outfile:
    #     outfile.write("\n".join(str(item) for item in F))
    # with open("Output_ptr.txt", "w") as outfile:
    #     outfile.write("\n".join(str(item) for item in ptr))
    # sys.exit()
    # print(F[n][m])
    # # print(F[n+1][m+1])
    # print(F[n][m+1])
    # print(F[n+1][m])
    # print(ptr[n][m])
    # # print(ptr[n+1][m+1])
    # print(ptr[n][m+1])
    # print(ptr[n+1][m])
    # print(data_1)
    # print(data_2)
    while n > start_1 or m > start_2:
        # print(start_1, start_2)
        print("Cyle: " + str(n) + " " + str(m))
        print("PTR " + str(ptr[n, m]))
        # print(ptr[n, m])
        # print(data_2[header2][m])
        # print(data_2[header2][m - 1])
        if ptr[n, m] == 'D':
            r_data_1.append(data_1[header1][n - 1])
            r_data_2.append(data_2[header2][m - 1])
            # r_data_1.append(data_1[n-1])
            # r_data_2.append(data_2[m-1])
            n -= 1
            m -= 1
            # print("D: " + str(n))
            # print("D: " + str(m))
            # print(F[n][m])
            # print(F[n + 1][m + 1])
            # print(F[n][m + 1])
            # print(F[n + 1][m])
            # print(ptr[n][m])
            # print(ptr[n + 1][m + 1])
            # print(ptr[n][m + 1])
            # print(ptr[n + 1][m])
        elif ptr[n, m] == 'H':
            r_data_1.append(data_1[header1][n - 1])
            # r_data_1.append(data_1[n-1])
            r_data_2.append('-')
            n -= 1
            # print("H: " + str(m))
            # print("H: " + str(n))
            # print(F[n][m])
            # print(F[n + 1][m + 1])
            # print(F[n][m + 1])
            # print(F[n + 1][m])
            # print(ptr[n][m])
            # print(ptr[n + 1][m + 1])
            # print(ptr[n][m + 1])
            # print(ptr[n + 1][m])
        elif ptr[n, m] == 'V':
            r_data_1.append('-')
            r_data_2.append(data_2[header2][m - 1])
            # r_data_2.append(data_2[m-1])
            m -= 1
            # print("V: " + str(m))
            # print("V: " + str(n))
            # print(F[n][m])
            # print(F[n + 1][m + 1])
            # print(F[n][m + 1])
            # print(F[n + 1][m])
            # print(ptr[n][m])
            # print(ptr[n + 1][m + 1])
            # print(ptr[n][m + 1])
            # print(ptr[n + 1][m])
        else:
            break
        print("PTR " + str(ptr[n, m]))
        # sys.exit()
    # print(r_data_1)
    r_data_1.reverse()
    r_data_2.reverse()

    # print(*r_data_1, sep='')
    # print(*r_data_2, sep='')
    return r_data_1, r_data_2


_init_(data_1, data_2)

# f = open("Fly_HOX.fa", 'r')
# for l, i in enumerate(f):
#     if l != 0:
#         data += i
#         data_1 = {"Fly_HOX.fa": data}
# f.close()
# print(data_1)

#
# f = open(fasta2, 'r')
# for l, i in enumerate(f)
#     if l != 0:
#         data_2 += i
# f.close()


# if __name__ == '__main__':
#     args = sys.argv
#     if len(args) == 6:
#         if (args[1][-3:] == '.py'):
#             if args[2][-6:] == '.fasta':
#                 if args[3][-6:] == '.fasta':
#                     if args[4][-4:] == '.txt':
#                         NWA(args)
#                     else:
#                         print('Fourth argument is optional and should be matches file ".txt"')
#                 else:
#                     print('Second argument should set to first FASTA file ".fasta"')
#             else:
#                 print('Second argument should set to first FASTA file ".fasta"')
#         else:
#             print('First argument should set "program_name.py"')
#     else:
#         print('The command for calling the program must be of form: program_name seq1.fasta seq2.fasta [matches.txt]')

