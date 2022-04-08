# -*- coding: utf-8 -*-
# @Author: WuLC
# @Date:   2017-04-02 08:52:24
# @Last Modified by:   WuLC
# @Last Modified time: 2017-04-02 09:50:31

###########################################################################################################
# Viterbi Algorithm for HMM
# dp, time complexity O(mn^2), m is the length of sequence of observation, n is the number of hidden states
##########################################################################################################

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

# five elements for HMM
# states = ('Healthy', 'Fever')
#
# observations = ('normal', 'cold', 'dizzy')
#
# start_probability = {'Healthy': 0.6, 'Fever': 0.4}
#
# transition_probability = {
#     'Healthy': {'Healthy': 0.7, 'Fever': 0.3},
#     'Fever': {'Healthy': 0.4, 'Fever': 0.6},
# }
#
# emission_probability = {
#     'Healthy': {'normal': 0.5, 'cold': 0.4, 'dizzy': 0.1},
#     'Fever': {'normal': 0.1, 'cold': 0.3, 'dizzy': 0.6},
# }
#

def Viterbit(obs, states, s_pro, t_pro, e_pro):
    path = {s: [] for s in states}  # init path: path[s] represents the path ends with s
    print("Path: " + str(path))
    curr_pro = {}
    for s in states:
        curr_pro[s] = sigmoid(s_pro.loc[s].loc['probs'] * e_pro.loc[s].loc[obs[0]])
        # print("e_pro[s][obs[0]] = " + str(e_pro[s][obs[0]]))
        # print("CP: " + str(curr_pro))
    for i in range(1, len(obs)):
        # print("curr_pro2: " + str(curr_pro))
        last_pro = curr_pro
        curr_pro = {}
        for curr_state in states:
            max_pro, last_sta = max(
                ((sigmoid(last_pro[last_state] * t_pro.loc[last_state].loc[curr_state] * e_pro.loc[curr_state].loc[obs[i]]), last_state)
                 for last_state in states))

            curr_pro[curr_state] = max_pro
            path[curr_state].append(last_sta)

    # find the final largest probability
    max_pro = -1
    max_path = None
    for s in states:
        path[s].append(s)
        print(path)
        if curr_pro[s] > max_pro:
            max_path = path[s]
            print("Max_path: " + str(max_path))
            max_pro = curr_pro[s]
    # print '%s: %s'%(curr_pro[s], path[s]) # different path and their probability
    return max_path


def sigmoid(x):
    if x >= 0:
        z = np.exp(-x)
        return 1 / (1 + z)
    else:

        z = np.exp(x)
        return z / (1 + z)


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
    # print(training_match__)
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
    performance = f"HMM_bw Performance: {nl} Performance vs Training Set: {nl} %of loops(_) matched: {training_match__} " \
                  f"{nl} %of sheets(e) matched: {training_match_e} {nl} %of helixes(h) matched: {training_match_h} {nl} " \
                  f"Total % matched: {training_match_perf} {nl} {nl}Performance vs Query Set: {nl} %of loops(_) matched:" \
                  f" {query_match__} {nl} %of sheets(e) matched: {query_match_e} {nl} %of helixes(h) matched: " \
                  f"{query_match_h}{nl} Total % matched: {query_match_perf}"

    return performance



def import_query(file_in):
    start_reading = False
    q_seq = []
    ks_q = []
    with open(file_in) as f:
        for l, i in enumerate(f):
            if l > 450: break

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


def import_training(file_in):
    start_reading = False
    seq = []
    ks = []
    seq_cur = ''
    ks_cur = ''
    ks_prev = ''
    seq_count = 0
    first_emit = False
    # Totals
    total_len = 0

    with open(file_in) as f:
        for l, i in enumerate(f):
            # if l > 500: break
            if i.rstrip() == "<>":  # look for start of relevant information
                start_reading = True
                first_emit = True
                seq_count += 1
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
                if first_emit == True:
                    first_emit = False
                    start_prob_matrix.at[ks_cur, 'probs'] += 1
                """    """

            ks_prev = ks_cur
            seq_prev = seq_cur
    # print(seq_count)
    # print(total_len)
    # print(trans_prob_matrix)
    """" V2 """
    # Normalize Transition Matrix
    trans_matrix = trans_prob_matrix.div(total_len)
    # Normalize Emission Matrix
    emit_matrix = emit_prob_matrix.div(total_len)

    # Normalize Start Probability Matrix
    start_matrix = start_prob_matrix.div(seq_count)
    for _ in range(len(start_prob_matrix)):
        if start_matrix.iloc[_].loc['probs'] == 0:
            start_matrix.iloc[_].loc['probs'] += .00000001

    """    """

    """ BW """
    # Note I am setting end state probabilities to start state probabilities.
    # Right now both are just the probability of a state existing.
    # Ideally would be derived differently from each training sample.
    end_matrix = start_matrix
    """    """

    return seq, ks, total_len, trans_matrix, emit_matrix, start_matrix, end_matrix

seq, ks, total_len, trans_matrix, emit_matrix, start_matrix, end_matrix = import_training("protein-secondary-structure.train")
q_seq, ks_q = import_query("protein-secondary-structure.test")
obs = seq
states = states
start_probability = start_matrix
transition_probability = trans_matrix
emission_probability = emit_matrix

bw_test_pred_state = Viterbit(seq, states, start_matrix, transition_probability, emission_probability)
bw_train_pred_state = Viterbit(q_seq, states, start_matrix, transition_probability, emission_probability)

print(solution_check(bw_test_pred_state, ks, bw_train_pred_state, ks_q))


# print(Viterbit(obs, states, start_probability, transition_probability, emission_probability))
# print(solution_check(obs,ks, q_seq, ks_q))