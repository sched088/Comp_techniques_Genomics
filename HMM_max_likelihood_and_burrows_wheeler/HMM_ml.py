"""
1) Load in Training Data
2) Calculate emission and transition probabilities
"""
import random

import numpy as np
import pandas as pd
pd.set_option('display.max_columns', None)

file_in = "protein-secondary-structure.train"
alphabet = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
states = ['_', 'e', 'h']
alpha__ = []
alpha_e = []
alpha_h = []

for l in alphabet:
    alpha__.append(l + "__")
    alpha_e.append(l + "_e")
    alpha_h.append(l + "_h")

matrix = np.zeros((20, 20))
raw_theta__ = pd.DataFrame(matrix, columns=alpha__, index=alpha__)
raw_theta_e = pd.DataFrame(matrix, columns=alpha_e, index=alpha_e)
raw_theta_h = pd.DataFrame(matrix, columns=alpha_h, index=alpha_h)


# raw_theta__ = pd.DataFrame(matrix, columns=alphabet, index=alphabet)
# raw_theta_e = pd.DataFrame(matrix, columns=alphabet, index=alphabet)
# raw_theta_h = pd.DataFrame(matrix, columns=alphabet, index=alphabet)


# Read in Sequence (x; var = seq) and Known States (pi; var = ks)
def import_training(file_in):
    start_reading = False
    seq = []
    ks = []
    seq_cur = ''
    seq_prev = ''
    ks_cur = ''
    ks_prev = ''

    # Totals
    total_len = 0
    total_trans = 0
    total__ = 0
    emission__ = []
    total_e = 0
    emission_e = []
    total_h = 0
    emission_h = []

    # Transitions
    A__e = 0
    A_e_ = 0
    A__h = 0
    A_h_ = 0
    A_eh = 0
    A_he = 0

    with open(file_in) as f:
        for l, i in enumerate(f):
            if i.rstrip() == "<>":  # look for start of relevant information
                start_reading = True
                continue
            if i.rstrip() == "<end>" or i.rstrip() == "end":  # look for end of relevant information
                start_reading = False
                continue
            if start_reading:
                total_len += 1
                seq += i[0]
                ks += i[-2]

                seq_cur = i[0]
                ks_cur = i[-2]
                # print(seq, ks)
                # Since already iterating through might as well be efficient and pull probabilities
                if ks_cur == "_":
                    total__ += 1
                    emission__ += i[0]
                if ks_cur == "e":
                    total_e += 1
                    emission_e += i[0]
                if ks_cur == "h":
                    total_h += 1
                    emission_h += i[0]

                if ks_prev == ks_cur:
                    if ks_cur == "_":
                        raw_theta__[seq_cur + "__"][seq_prev + "__"] += 1
                    if ks_cur == "e":
                        raw_theta_e[seq_cur + "_e"][seq_prev + "_e"] += 1
                    if ks_cur == "h":
                        raw_theta_h[seq_cur + "_h"][seq_prev + "_h"] += 1

                if ks_prev != ks_cur:
                    total_trans += 1
                    if ks_prev == "_":
                        if ks_cur == "e":
                            A__e += 1
                        elif ks_cur == "h":
                            A__h += 1
                    if ks_prev == "e":
                        if ks_cur == "_":
                            A_e_ += 1
                        elif ks_cur == "h":
                            A_eh += 1
                    if ks_prev == "h":
                        if ks_cur == "e":
                            A_he += 1
                        elif ks_cur == "_":
                            A_h_ += 1
            ks_prev = ks_cur
            seq_prev = seq_cur

    return seq, ks, total_len, A__e, A__h, A_e_, A_eh, A_he, A_h_, total__, total_e, total_h, \
           emission__, emission_e, emission_h, raw_theta__, raw_theta_h, raw_theta_e


def prob_table(total_len, A__e, A__h, A_e_, A_eh, A_he, A_h_, total__, total_e, total_h,
               emission__, emission_e, emission_h, raw_theta__, raw_theta_h, raw_theta_e):
    # Nomenclature: transition from state(sheet) to state(helix) : t_eh
    t__e = A__e / total_len
    t__h = A__h / total_len
    t___ = 1 - t__h - t__e

    t_e_ = A_e_ / total_len
    t_eh = A_eh / total_len
    t_ee = 1 - t_e_ - t_eh

    t_h_ = A_h_ / total_len
    t_he = A_he / total_len
    t_hh = 1 - t_h_ - t_he


    # Between State Transition matrices (row = from; column = to)
    # Improvement Opportunity: I am not tracking what they change to.
    # It is the same probability that each state change goes to a random AA.
    matrix_t__e = np.full((20, 20), t__e / 20)
    theta_t__e = pd.DataFrame(matrix_t__e, index=alpha__, columns=alpha_e)
    matrix_t__h = np.full((20, 20), t__h / 20)
    theta_t__h = pd.DataFrame(matrix_t__h, index=alpha__, columns=alpha_h)
    matrix_t___ = np.full((20, 20), t___ / 20)
    theta_t___ = pd.DataFrame(matrix_t___, index=alpha__, columns=alpha__)  # Why aren't these being used?

    matrix_t_e_ = np.full((20, 20), t_e_ / 20)
    theta_t_e_ = pd.DataFrame(matrix_t_e_, index=alpha_e, columns=alpha__)
    matrix_t_eh = np.full((20, 20), t_eh / 20)
    theta_t_eh = pd.DataFrame(matrix_t_eh, index=alpha_e, columns=alpha_h)
    matrix_t_ee = np.full((20, 20), t_ee / 20)
    theta_t_ee = pd.DataFrame(matrix_t_ee, index=alpha_e, columns=alpha_e)

    matrix_t_h_ = np.full((20, 20), t_h_ / 20)
    theta_t_h_ = pd.DataFrame(matrix_t_h_, index=alpha_h, columns=alpha__)
    matrix_t_he = np.full((20, 20), t_he / 20)
    theta_t_he = pd.DataFrame(matrix_t_he, index=alpha_h, columns=alpha_e)
    matrix_t_hh = np.full((20, 20), t_hh / 20)
    theta_t_hh = pd.DataFrame(matrix_t_hh, index=alpha_h, columns=alpha_h)

    # Normalize the Theta Matrices and account for transition to self
    theta__ = raw_theta__.div(total__)
    theta__ = theta__.multiply(theta_t___)
    theta_e = raw_theta_e.div(total_e)
    theta_e = theta_e.multiply(theta_t_ee)
    theta_h = raw_theta_h.div(total_h)
    theta_h = theta_h.multiply(theta_t_hh)


    # Build Theta Matrix
    theta1 = [theta__, theta_t__e, theta_t__h]  # I think I need to look at the theta__, theta_e, and theta_h to see if they are accurate. Why not use theta_t___? I may need to multiply them by something too?
    theta2 = [theta_t_e_, theta_e, theta_t_eh]
    theta3 = [theta_t_h_, theta_t_he, theta_h]

    theta1a = pd.concat(theta1, axis = 1)
    theta2a = pd.concat(theta2, axis = 1)
    theta3a = pd.concat(theta3, axis = 1)

    theta = pd.concat([theta1a, theta2a, theta3a])

    # Probabilities for starting
    prob__ = total__ / total_len
    prob_e = total_e / total_len
    prob_h = total_h / total_len
    starting_prob = [prob__, prob_e, prob_h]

    # print(theta)
    # print(theta["A__"]["A__"])
    return theta, starting_prob

def viterbi(theta, starting_prob, states, total_len):
    # Initialize: Multiply the initial probability of state i (1/3) [ASSUMPTION] by the emission probability from state i to observable O at time t = 1.
    start_state = ''
    # rand = random.uniform(0.0, 1.0)
    #
    # if rand <= prob__: start_state = '_'
    # if rand <= prob__ + prob_e: start_state = 'e'
    # if rand > prob__ + prob_e: start_state = 'h'
    current_prob = {}
    ptr = {state: [] for state in states}

    for state in states:
        current_prob[state] = starting_prob[states.index(state)] * theta[seq[0]+"_"+state][seq[0]+"_"+state]
    for _ in range(1, total_len):
        prev_prob = current_prob
        current_prob = {}
        for state in states:
            V, previous_state = max(((prev_prob[prev_state] * theta[seq[_]+"_"+prev_state][seq[_]+"_"+state], prev_state) for prev_state in states))
            current_prob[state] = V
            ptr[state].append(previous_state)

    V = -1
    ptr_max = None
    for state in states:
        ptr[state].append(state)
        if current_prob[state] > V:
            ptr_max = ptr[state]
            print(ptr_max)
            V = current_prob[state]

    return ptr_max



seq, ks, total_len, A__e, A__h, A_e_, A_eh, A_he, A_h_, total__, total_e, total_h, \
emission__, emission_e, emission_h, raw_theta__, raw_theta_h, raw_theta_e = import_training(file_in)

theta, starting_prob = prob_table(total_len, A__e, A__h, A_e_, A_eh, A_he, A_h_, total__, total_e, total_h,
                                           emission__, emission_e, emission_h, raw_theta__, raw_theta_h, raw_theta_e)

hidden_states = viterbi(theta, starting_prob, states, total_len)
print(hidden_states)