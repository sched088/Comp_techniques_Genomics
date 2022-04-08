"""
1) Load in Training Data
2) Calculate emission and transition probabilities

V1 = messy start
v2 = clean up emission and transition matrices. Remove large Theta table in preference of both emission and transition
"""

import random

import numpy as np
import pandas as pd
pd.set_option('display.max_columns', None)

training_file = "protein-secondary-structure.train"
query_file = "protein-secondary-structure.test"
alphabet = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
states = ['_', 'e', 'h']

""" V2 """
trans_prob = {state: {s: [] for s in states} for state in states}
emit_prob = {state: {a: [] for a in alphabet} for state in states}

# Creating Transition Matrix
trans_prob_matrix = np.zeros((3, 3))
trans_prob_matrix = pd.DataFrame(trans_prob_matrix, index=states, columns=states)

# Creating Emission Matrix
emit_prob_matrix = np.zeros((3, 20))
emit_prob_matrix = pd.DataFrame(emit_prob_matrix, index=states, columns=alphabet)

# Starting Probability Matrix
start_prob_matrix = np.zeros((3, 1))
start_prob_matrix = pd.DataFrame(start_prob_matrix, index=states, columns=['probs'])
"""    """


# Read in Sequence (x; var = seq) and Known States (pi; var = ks)
# One thing I realized that I am not doing is looking at the 'start' state of each training sequence and resetting /
     # for each restart. I am currently only reading this in as one long string which is not fully accurate.
def import_training(file_in):
    start_reading = False
    seq = []
    ks = []
    seq_cur = ''
    ks_cur = ''
    ks_prev = ''

    # Totals
    total_len = 0

    with open(file_in) as f:
        for l, i in enumerate(f):
            if i.rstrip() == "<>":  # look for start of relevant information
                start_reading = True
                continue
            if i.rstrip() == "<end>" or i.rstrip() == "end":  # look for end of relevant information
                start_reading = False
                continue
            if start_reading:
                if ks_prev == '': ks_prev = i[-2]

                total_len += 1
                seq += i[0]
                ks += i[-2]

                seq_cur = i[0]
                ks_cur = i[-2]

                # Since already iterating through might as well be efficient and pull probabilities

                """" V2 """
                # Populate Transition Matrix
                trans_prob_matrix.at[ks_prev, ks_cur] += 1

                # Populate Emission Matrix
                emit_prob_matrix.at[ks_cur, seq_cur] += 1

                # Populate Start Probability Matrix
                start_prob_matrix.at[ks_cur, 'probs'] += 1
                """    """

            ks_prev = ks_cur
            seq_prev = seq_cur

    """" V2 """
    # Normalize Transition Matrix
    trans_matrix = trans_prob_matrix.div(total_len)

    # Normalize Emission Matrix
    emit_matrix = emit_prob_matrix.div(total_len)

    # Normalize Start Probability Matrix
    start_matrix = start_prob_matrix.div(total_len)
    """    """

    return seq, ks, total_len, trans_matrix, emit_matrix, start_matrix


def import_query(file_in):
    start_reading = False
    q_seq = []
    ks_q = []
    with open(file_in) as f:
        for l, i in enumerate(f):
            if i.rstrip() == "<>":  # look for start of relevant information
                start_reading = True
                continue
            if i.rstrip() == "<end>" or i.rstrip() == "end":  # look for end of relevant information
                start_reading = False
                continue
            if start_reading:
                q_seq += i[0]
                ks_q += i[-2]

    return q_seq, ks_q

def viterbi2(seq, states, start_matrix, trans_matrix, emit_matrix):
    trellis = np.zeros((len(states), len(seq)))
    pointers = np.zeros((len(states), len(seq)))

    for s in range(len(states)):
        trellis[s, 0] = start_matrix.iloc[s].loc['probs'] * emit_matrix.iloc[s].loc[seq[0]]
    for o in range(1, len(seq)):

        for s in range(len(states)):
            k = np.argmax(trellis[k, o-1] * trans_matrix.iloc[k].iloc[s] * emit_matrix.iloc[s].loc[seq[o]] for k in trellis)
            trellis[s, o] = trellis[k, o-1] * trans_matrix.iloc[k].iloc[s] * emit_matrix.iloc[s].loc[seq[o]]
            pointers[s, o] = k
    best_path = []
    k = np.argmax(trellis[k, len(seq)-1])

    for o in range(len(seq)-1, -1, -1):

        best_path += states[int(k)]
        # best_path.insert(int(0), states[k])
        k = pointers[int(k), o]
    return best_path


def viterbi(seq, trans_matrix, emit_matrix, start_matrix):

    current_prob = {}
    ptr = {state: [] for state in states}

    for state in states:
        current_prob[state] = start_matrix.loc[state, 'probs'] * emit_matrix.loc[state, seq[0]]
    for _ in range(1, len(seq)):
        prev_prob = current_prob
        current_prob = {}
        for state in states:
            V, previous_state = max(((prev_prob[prev_state] * trans_matrix.loc[prev_state, state] * emit_matrix.at[state, seq[_]], prev_state) for prev_state in states))
            current_prob[state] = V
            ptr[state].append(previous_state)

    V = -5
    ptr_max = None
    for state in states:
        ptr[state].append(state)
        if current_prob[state] > V:
            ptr_max = ptr[state]
            V = current_prob[state]

    return ptr_max

def solution_check(seq, ks, q_seq, ks_q):
    training_match__ = 0
    training_match_e = 0
    training_match_h = 0
    total_training_match = 0
    training_match_perf = 0

    query_match__ = 0
    query_match_e = 0
    query_match_h = 0
    total_query_match = 0
    query_match_perf = 0

    for _ in range(0, len(seq)):
        if seq[_] == ks[_]: total_training_match += 1
        if seq[_] == ks[_] and ks[_] == "_":
            training_match__ += 1
        if seq[_] == ks[_] and ks[_] == "e":
            training_match_e += 1
        if seq[_] == ks[_] and ks[_] == "h":
            training_match_h += 1
    print(training_match__)
    training_match__ = (training_match__ / ks.count('_')) * 100
    training_match_e = (training_match_e / ks.count('e')) * 100
    training_match_h = (training_match_h / ks.count('h')) * 100
    training_match_perf = (total_training_match / len(seq)) * 100

    for _ in range(0, len(q_seq)):
        if q_seq[_] == ks_q[_]: total_query_match += 1
        if q_seq[_] == ks_q[_] and ks_q[_] == "_":
            query_match__ += 1
        if q_seq[_] == ks_q[_] and ks_q[_] == "e":
            query_match_e += 1
        if q_seq[_] == ks_q[_] and ks_q[_] == "h":
            query_match_h += 1

    query_match__ = (query_match__ / ks_q.count('_')) * 100
    query_match_e = (query_match_e / ks_q.count('e')) * 100
    query_match_h = (query_match_h / ks_q.count('h')) * 100
    query_match_perf = (total_query_match / len(q_seq)) * 100

    print(seq)
    print(ks)
    print("")
    print(q_seq)
    print(ks_q)

    nl = '\n'
    performance = f"HMM_ml Performance: {nl} Performance vs Training Set: {nl} %of loops(_) matched: {training_match__} " \
                  f"{nl} %of sheets(e) matched: {training_match_e} {nl} %of helixes(h) matched: {training_match_h} {nl} " \
                  f"Total % matched: {training_match_perf} {nl} {nl}Performance vs Query Set: {nl} %of loops(_) matched:" \
                  f" {query_match__} {nl} %of sheets(e) matched: {query_match_e} {nl} %of helixes(h) matched: " \
                  f"{query_match_h}{nl} Total % matched: {query_match_perf}"

    return performance

seq, ks, total_len, trans_matrix, emit_matrix, start_matrix = import_training(training_file)
q_seq, ks_q = import_query(query_file)


# test_pred_state = viterbi(seq, trans_matrix, emit_matrix, start_matrix)
# train_pred_state = viterbi(q_seq, trans_matrix, emit_matrix, start_matrix)
#
# print(solution_check(test_pred_state, ks, train_pred_state, ks_q))

test_pred_state = viterbi2(seq, states, start_matrix, trans_matrix, emit_matrix)
train_pred_state = viterbi2(q_seq, states, start_matrix, trans_matrix, emit_matrix)

print(solution_check(test_pred_state, ks, train_pred_state, ks_q))

