"""
This is the first homework for class 5481 Computational Genomics
Max Scheder
"""

import sys
import numpy

np = numpy
args = sys.argv


def _import(args):  # this code is confusing to me
    _, file_name, fasta1, fasta2, matches = args
    data_1 = []
    data_2 = []
    data_match = ''
    data = ''
    header = ''

    with open(fasta1) as f:
        header1 = f.readline().rstrip()
        data = ''
        for l, i in enumerate(f):
            data += i.rstrip()
        # data_1 = {header1: data}
    data_1[:0] = data
    
    with open(fasta2) as f:
        header2 = f.readline().rstrip()
        data = ''
        for l, i in enumerate(f):
            data += i.rstrip()
        # data_2 = {header2: data}
    data_2[:0] = data
    # print(data_2)
    # print(data_2[-1])
    # print(str(data_1) + '\n' + str(data_2))
    nwa(data_1, data_2)
    return data_1, data_2


def nwa(data_1, data_2, s = -3, m = 1, d = -2):
    ld1 = len(data_1)
    ld2 = len(data_2)
    """ Initialization:
        F(0,0) = 0
        F(0,j) = -j x d
        F(i,0) = -i x d
    """
    F = np.zeros((ld1 + 1, ld2 + 1)) #create array/matrix of sequences being compared
    F[:, 0] = np.linspace(0, -ld1*2, ld1 + 1) #alternately could do loop through rows/columns and add penalty * i
    F[0, :] = np.linspace(0, -ld2*2, ld2 + 1) #alternately could do loop through rows/columns and add penalty * i
    # print(F)
    """ Pointers:
                    { Diag = (1,1)  if [case 1]
        ptr(i,j) =  { Left = (1,0)  if [case 2]
                    { Up = (0,1)    if [case 3]
    """
    ptr = F.astype(str)
    ptr[:, 0] = 'V'
    ptr[0, :] = 'H'
    ptr[0, 0] = 0
    # print(ptr)

    """ Induction:
                   { F(i-1, j-1) + s(xi, yj)    [case 1]
     F(i,j) = max  { F(i-1, j) - d              [case 2]
                   { F(i, j-1) - d              [case 3]
                   
        This section populates two matrices, one with the values, the other with the pointer direction.
        Pointer direction not needed for assignment, but helped me visualize.
    """
    for i in range(1, ld1 + 1):
        # print("populate: " + str(i))
        for j in range(1, ld2 + 1):
            # print("populate: " + str(j))
            # print('hit')
            # if data_1[i-1] == data_2[j-1]:
            if data_1[i-1] == data_2[j-1]:
                diag = F[i-1][j-1] + m
            else:
                diag = F[i-1][j-1] + s
            hori = F[i-1][j] + d
            vert = F[i][j-1] + d
            F[i][j] = max(diag, vert, hori)
            if diag == max(diag, vert, hori):
                ptr[i][j] = 'D'
            elif vert == max(diag, vert, hori):
                ptr[i][j] = 'V'
            else:
                ptr[i][j] = 'H'
    # print(F)
    # print(ptr)
    """ backwards trace through optimal alignment"""
    n = ld1
    m = ld2
    r_data_1 = []
    r_data_2 = []
    # print(data_1)
    # print(data_2)
    while n > 0 or m > 0:
        # print("Cyle: " + str(n) + " " + str(m))
        # print("PTR " + str(ptr[n,m]))
        # print(ptr[n, m])
        # print(data_2[m])
        # print(data_2[m-1])
        if ptr[n, m] == 'D':
            r_data_1.append(data_1[n - 1])
            r_data_2.append(data_2[m - 1])

            n -= 1
            # print("D: " + str(n))
            m -= 1
            # print("D: " + str(m))
        elif ptr[n, m] == 'H':
            r_data_1.append(data_1[n - 1])
            r_data_2.append('-')
            n -= 1

        elif ptr[n, m] == 'V':
            r_data_1.append('-')
            r_data_2.append(data_2[m - 1])
            m -= 1

    r_data_1.reverse()
    r_data_2.reverse()

    alignment_score = int(F[ld1, ld2])

    print("Sequence 1: ")
    print(*r_data_1, sep='')
    print("")
    print("Sequence 2: ")
    print(*r_data_2, sep='')
    print("")
    print("Alignment Score: " + str(alignment_score))

    return r_data_1, r_data_2, alignment_score


# nwa(data_1, data_2)

if __name__ == '__main__':
    args = sys.argv
    match = False
    print("Your input arguments are: " + str(args))
    if len(args) == 4:
        args.append("no match")
    if args[1][-3:] != '.py':
        print('First argument should set "program_name.py"')
    if args[2][-3:] != '.fa':
        print('Second argument should set to first FASTA file ".fa"')
    if args[3][-3:] != '.fa':
        print('Third argument should set to second FASTA file ".fa"')
    if len(args) == 5:
        if args[4][-4:] == '.txt':
            match = True
    _import(args)
