def import_genes():
    with open("Chr21.txt") as f:
        header1 = f.readline().rstrip()
        # print("Header: " + header)
        gene_data = ''
        for l, i in enumerate(f):
            gene_data += i.rstrip()
        # data_1 = {header1: data}
        # print(len(data))
        print(len(gene_data))
        print(gene_data)
    return gene_data