"""
Write a program to find all the start-codon and stop-codon
and report all the ORFs longer than k (try multiple ks for the best result).

For example, you can try k = 60, 100, 200.

Compare your predictions with the ground-truth to make sure your program is working properly.
"""

"""
1) Read in file
2) Search Reading Frames for start codon
3) Count length of ORFs until stop codon
"""

# Start Codons: 'ATG'
# Stop Codons: 'TAA', 'TGA', 'TAG'


def import_sequence():
    with open("GCF_000005845.2_ASM584v2_genomic.fna") as f:
        header1 = f.readline().rstrip()
        # print("Header: " + header)
        data = ''
        for l, i in enumerate(f):
            data += i.rstrip()
        # data_1 = {header1: data}
    return data


def search_orf(data, k):
    """
    Frames starting at 0, 1, 2 forward and then len(data), len(data)-1, len(data)-2
    """

    orf_lst = []

    start_c = ['ATG', 'GTG', 'TTG']
    current_c_lst = []
    start_loc = 0
    start = False

    stop_c = ['TAA', 'TGA', 'TAG']
    # stop_c_lst = []  # not needed due to just adding the stop at the end of the start list.
    # stop_loc = 0  # not needed due to stop_loc = i+1 which is tracked in loop


#
# Forward Read First:
#
    rf = 3
    orf_c = 0

    rd = "-"
    while rf > 0:
        print(rd)
        print(rf)
        for i in range(len(data)-rf-1, rf, -3):
            codon = data[i]+data[i-1]+data[i-2]
            if start is False:
                if codon in start_c:
                    current_c_lst.extend([rf, rd, codon, i+1])
                    start_loc = i+1
                    start = True
            else:
                if codon in stop_c:
                    if (i+1 - start_loc) > k:
                        current_c_lst.insert(len(current_c_lst)-3, codon)
                        current_c_lst.extend([i+1])
                        current_c_lst.extend([(i+1 - start_loc)])
                        orf_lst.extend([current_c_lst])
                        current_c_lst = []
                        # print("ORF: " + str(orf_lst))
                        orf_c += 1
                        start = False
        rf -= 1

    rd = "+"
    while rf < 3:
        print(rd)
        print(rf)
        for i in range(rf, len(data)-rf, 3):
            codon = data[i:i+3]
            if start is False:
                if codon in start_c:
                    current_c_lst.extend([rf, rd, codon, i+1])
                    start_loc = i+1         
                    start = True
            else:
                if codon in stop_c:
                    if (i+1 - start_loc) > k:
                        current_c_lst.insert(len(current_c_lst)-3, codon)
                        current_c_lst.extend([i+1])
                        current_c_lst.extend([(i+1 - start_loc)])
                        orf_lst.extend([current_c_lst])
                        current_c_lst = []
                        # print("ORF: " + str(orf_lst))
                        orf_c +=1
                        start = False
        rf += 1

    # print(orf_lst)
    print(orf_c)
    return orf_lst, orf_c


def save_file(orf_list, orf_count, k):
    print("Exporting...")
    file = open("Assignment 2 Export", "w")
    file.write(f'ORF Count (length >{k}): {orf_count} \n {orf_list}')
    print("Done!")

if __name__ == "__main__":
    k = 600  # adjust ORF length here
    sequence = import_sequence()
    orf_list, orf_count = search_orf(sequence, k)
    save_file(orf_list, orf_count, k)
