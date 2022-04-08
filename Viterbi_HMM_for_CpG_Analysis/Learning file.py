import numpy as np
# import io_handler


def import_sequence():
    with open("seq1.fa") as f:
        header1 = f.readline().rstrip()
        # print("Header: " + header)
        data = ''
        for l, i in enumerate(f):
            data += i.rstrip()
        # data_1 = {header1: data}
        # print(len(data))
    return data


def viterbi(emission_probs, transition_probs, initial_dist, emissions):
    probs = emission_probs[:, 0] * initial_dist
    # print("probs" + str(probs))
    stack = []
    num_states = transition_probs.shape[0]

    for emission in emissions[1:]:
        val = ['A', 'C', 'G', 'T']
        trans_probs = transition_probs * np.row_stack(probs)
        # print("trans_probs" + str(trans_probs))
        max_col_ixs = np.argmax(trans_probs, axis=0)
        # print(emission_probs[:, val.index(emission)])
        probs = emission_probs[:, val.index(emission)] * trans_probs[max_col_ixs, np.arange(num_states)]
        # print("probs"+str(probs))
        # For avoiding probability reaching 0
        # Un comment the line and apply the appropriate value
        # Refer known issues for more information
        # probs = [x*4.25 for x in probs]


        stack.append(max_col_ixs)
    print(probs)
    state_seq = [np.argmax(probs)]
    print(state_seq)
    while stack:
        max_col_ixs = stack.pop()
        state_seq.append(max_col_ixs[state_seq[-1]])

    state_seq.reverse()

    return state_seq


emission_probability = np.array([
    [0, 0, 0, 0],
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 1],
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 1]
])

transition_probability = np.array([
    [0, 0.0725193, 0.163763, 0.1788242, 0.0754545, 0.1322050, 0.1267006, 0.1226380, 0.1278950],
    [0.001, 0.1762237, 0.2682517, 0.4170629, 0.1174825, 0.0035964, 0.0054745, 0.0085104, 0.0023976],
    [0.001, 0.1672435, 0.3599201, 0.267984, 0.1838722, 0.0034131, 0.0073453, 0.005469, 0.0037524],
    [0.001, 0.1576223, 0.3318881, 0.3671328, 0.1223776, 0.0032167, 0.0067732, 0.0074915, 0.0024975],
    [0.001, 0.0773426, 0.3475514, 0.375944, 0.1781818, 0.0015784, 0.0070929, 0.0076723, 0.0036363],
    [0.001, 0.0002997, 0.0002047, 0.0002837, 0.0002097, 0.2994005, 0.2045904, 0.2844305, 0.2095804],
    [0.001, 0.0003216, 0.0002977, 0.0000769, 0.0003016, 0.3213566, 0.2974045, 0.0778441, 0.3013966],
    [0.001, 0.0001768, 0.000238, 0.0002917, 0.0002917, 0.1766463, 0.2385224, 0.2914165, 0.2914155],
    [0.001, 0.0002477, 0.0002457, 0.0002977, 0.0002077, 0.2475044, 0.2455084, 0.2974035, 0.2075844]
])

initial_dist = np.array([0, 0.0725193, 0.163763, 0.1788242, 0.0754545, 0.1322050, 0.1267006, 0.1226380, 0.1278950])

emissions = import_sequence()

cpg_islands = viterbi(emission_probability, transition_probability, initial_dist, emissions)
print(cpg_islands)
cpg_islands = ["-" if i >= 5 else "+" for i in cpg_islands]
print(cpg_islands)
# io_handler.save_cpg_islands("cpg_islands.txt", " ".join(cpg_islands))