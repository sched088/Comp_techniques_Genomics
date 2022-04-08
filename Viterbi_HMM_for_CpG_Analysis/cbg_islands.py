"""
Utilize Viterbi algorithm

Initialize start probabilities a01 + ... + a0k = 1
(equal probabilities of starting with any emission)

"""

import numpy as np
import random

e = np.array([
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 1],
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 1]])

p = 0.98

q = 0.999

# hardcoded in from CSV document
a = np.array([
    [0.180 * p, 0.274 * p, 0.426 * p, 0.120 * p, (1 - p) / 4, (1 - p) / 4, (1 - p) / 4, (1 - p) / 4],
    [0.171 * p, 0.368 * p, 0.274 * p, 0.188 * p, (1 - p) / 4, (1 - p) / 4, (1 - p) / 4, (1 - p) / 4],
    [0.161 * p, 0.339 * p, 0.375 * p, 0.125 * p, (1 - p) / 4, (1 - p) / 4, (1 - p) / 4, (1 - p) / 4],
    [0.079 * p, 0.355 * p, 0.384 * p, 0.182 * p, (1 - p) / 4, (1 - p) / 4, (1 - p) / 4, (1 - p) / 4],
    [(1 - q) / 4, (1 - q) / 4, (1 - q) / 4, (1 - q) / 4, 0.300 * q, 0.205 * q, 0.285 * q, 0.210 * q],
    [(1 - q) / 4, (1 - q) / 4, (1 - q) / 4, (1 - q) / 4, 0.322 * q, 0.298 * q, 0.078 * q, 0.302 * q],
    [(1 - q) / 4, (1 - q) / 4, (1 - q) / 4, (1 - q) / 4, 0.248 * q, 0.246 * q, 0.298 * q, 0.208 * q],
    [(1 - q) / 4, (1 - q) / 4, (1 - q) / 4, (1 - q) / 4, 0.177 * q, 0.239 * q, 0.292 * q, 0.292 * q]])

# transition_probability = np.array([
#     [0, 0.0725193, 0.163763, 0.1788242, 0.0754545, 0.1322050, 0.1267006, 0.1226380, 0.1278950],
#     [0.001, 0.1762237, 0.2682517, 0.4170629, 0.1174825, 0.0035964, 0.0054745, 0.0085104, 0.0023976],
#     [0.001, 0.1672435, 0.3599201, 0.267984, 0.1838722, 0.0034131, 0.0073453, 0.005469, 0.0037524],
#     [0.001, 0.1576223, 0.3318881, 0.3671328, 0.1223776, 0.0032167, 0.0067732, 0.0074915, 0.0024975],
#     [0.001, 0.0773426, 0.3475514, 0.375944, 0.1781818, 0.0015784, 0.0070929, 0.0076723, 0.0036363],
#     [0.001, 0.0002997, 0.0002047, 0.0002837, 0.0002097, 0.2994005, 0.2045904, 0.2844305, 0.2095804],
#     [0.001, 0.0003216, 0.0002977, 0.0000769, 0.0003016, 0.3213566, 0.2974045, 0.0778441, 0.3013966],
#     [0.001, 0.0001768, 0.000238, 0.0002917, 0.0002917, 0.1766463, 0.2385224, 0.2914165, 0.2914155],
#     [0.001, 0.0002477, 0.0002457, 0.0002977, 0.0002077, 0.2475044, 0.2455084, 0.2974035, 0.2075844]
# ])

# probs = sum(transition_probability[1:, 1])/7
# print(probs)
# print(a.shape[0])


def viterbi(a, e, seq):
    start_b = random.randint(0, 3)
    emission_prob = e[:, start_b] * a[start_b, :]
    pi = []
    # emissions = (['A+', 'C+', 'G+', 'T+', 'A-', 'C-', 'G-', 'T-'])
    # only 4 letters in alphabet, ATCG -> then possible emissions are 2 for each one of those (A+ or A-)
    # if start_k <= p / (p + q):
    #     start_state = random.choice(emissions[:4])
    # else:
    #     start_state = random.choice(emissions[4:])
    # print(emission_prob)
    # V0 = 1
    # Vi = a * np.row_stack(emission_prob)

    # print(start_state)
    # print(emissions.index(start_state))
    for x in seq[1:]:
        if x == 'N':
            continue
        val = ['A', 'C', 'G', 'T']
        # Vi = e[:, val.index(x)] * max(a)
        # Vi = np.log(a) * np.log(np.row_stack(emission_prob))
        Vi = a * np.row_stack(emission_prob)
        max_Vi = np.argmax(Vi, axis = 0)
        # print("Vi" + str(Vi))
        # print("max_Vi" + str(max_Vi))
        Vi = abs(np.ma.log(Vi))  # *-1
        # print("Vi_log" + str(Vi))
        max_Vi = np.argmax(Vi, axis = 0)
        # print("max_Vi_log" + str(max_Vi))
        # emission_prob = np.log(e[:, val.index(x)]) * a[ptri, np.linspace(0, a.shape[0]-1, num=a.shape[0]).astype(int)]
        emission_prob = e[:, val.index(x)] * Vi[max_Vi, np.linspace(0, a.shape[0]-1, num=a.shape[0]).astype(int)]
        print("e_prob" + str(emission_prob))
        # emission_prob = abs(np.ma.log(emission_prob))  # *-1
        # print("e_prob_log" + str(emission_prob))
        # print("ls" + str(np.linspace(1, a.shape[0], a.shape[0]).astype(int)))
        # print("Vi" + str(Vi))
        # print("max_vi" + str(max_Vi))
        # print("emission_prob" + str(emission_prob))

        pi.append(max_Vi)
    ptrs = [np.argmax(emission_prob)]
    # print(ptrs)
    while pi:
        max_Vi = pi.pop()
        print(max_Vi)
        ptrs.append(max_Vi[ptrs[-1]])
        print(ptrs)
    ptrs.reverse()
    print(len(ptrs))
    print(ptrs)
    return ptrs


def import_sequence():
    with open("seq1.fa") as f:
        header1 = f.readline().rstrip()
        # print("Header: " + header)
        data = ''
        for l, i in enumerate(f):
            data += i.rstrip()
        # data_1 = {header1: data}
        # print(len(data))
        print(len(data))
        print(data)
    return data


seq = import_sequence()
start = viterbi(a, e, seq)
start = ["-" if i >= 4 else "+" for i in start]
print(start)
