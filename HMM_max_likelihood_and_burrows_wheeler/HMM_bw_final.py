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
    seq_count = 0
    first_emit = False
    # Totals
    total_len = 0

    with open(file_in) as f:
        for l, i in enumerate(f):
            if l > 1000: break
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


def import_query(file_in):
    start_reading = False
    q_seq = []
    ks_q = []
    with open(file_in) as f:
        for l, i in enumerate(f):
            if l > 900: break

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


def forward_procedure(seq, q_seq, start_matrix, end_matrix, emit_matrix, trans_matrix):
    forward_prob = np.zeros((len(states), len(q_seq)))
    total_probs = []
    position = -1

    for emission in q_seq:
        position += 1
        for state_ in range(len(states)):
            if position == 0:

                forward_prob[state_, position] = start_matrix.iloc[state_].loc['probs'] * emit_matrix.iloc[state_].loc[emission]
            else:
                for s_ in range(len(states)):

                    total_probs.append(forward_prob[s_, position - 1] * emit_matrix.iloc[state_].loc[emission] * \
                                  trans_matrix.iloc[s_, state_])  # might be wrong emit_matrix.loc[state_, seq[emission]]

                forward_prob[state_, position] = sum(total_probs)

    end_prob = np.multiply(forward_prob[:, -1:], end_matrix)
    end_totals = end_prob.sum()

    return forward_prob, end_totals


def backward_procedure(seq, q_seq, start_matrix, end_matrix, emit_matrix, trans_matrix):
    backward_prob = np.zeros((len(states), len(q_seq)))
    total_probs = []
    position = 0
    for emission in range(1,len(q_seq)+1):
        position += 1
        for state_ in range(len(states)):
            if -position == -1:
                backward_prob[emission, -position] = end_matrix.iloc[state_].loc['probs']
            else:
                for s_ in range(len(states)):
                    total_probs.append(backward_prob[s_, -position + 1] * emit_matrix.iloc[s_].loc[q_seq[-position + 1]] * \
                                  trans_matrix.iloc[state_, s_])  # might be wrong emit_matrix.loc[s_, seq[-position+1]]
                backward_prob[state_, -position] = sum(total_probs)
    start_prob = []
    for _ in range(len(states)):
        start_prob.append(backward_prob[_, 0] * emit_matrix.iloc[_].loc[q_seq[0]])

    start_totals = sum(start_prob)

    # print("backward")
    # print(backward_prob)
    # print(start_totals)
    # print("----")
    # print(trans_matrix)
    # print(emit_matrix)
    # print(start_matrix)
    # print(end_matrix)

    return backward_prob, start_totals  # start_totals not needed?


def s_g_probabilities(seq, q_seq, forward_prob, backward_prob, end_totals, trans_matrix, emit_matrix):
    s_prob = np.zeros((len(states), len(q_seq) - 1, len(states)))
    g_prob = np.zeros((len(states), len(q_seq)))
    print("FP")
    print(forward_prob)
    print("BP")
    print(backward_prob)
    print("TM")
    print(trans_matrix)
    print("EM")
    print(emit_matrix)
    for _ in range(len(q_seq) - 1):
        for i in range(len(states)):
            for j in range(len(states)):
                s_prob[i, _, j] = (forward_prob[i, _] * backward_prob[j, _ + 1] * trans_matrix.iloc[i, j] * emit_matrix.iloc[
                    j].loc[q_seq[_ + 1]]) / end_totals

    for _ in range(len(q_seq)):
        for i in range(len(states)):
            g_prob[i, _] = (forward_prob[i, _] * backward_prob[i, _] / end_totals)

    print("s_probs")
    print(s_prob)
    print("g_prob")
    print(g_prob)
    # print("----")
    # print(trans_matrix)
    # print(emit_matrix)

    return s_prob, g_prob


def a_b_mxs(seq, q_seq, s_prob, g_prob):
    a_star = np.zeros((len(states), len(states)))
    b_star = np.zeros((len(states), len(alphabet)))
    denom_a = []
    # numer_b = []
    # denom_b = []
    for _ in range(len(states)):
        for i in range(len(states)):
            for j in range(len(q_seq) - 1):
                # print("HERE")
                # print(a_star)
                a_star[_, i] = a_star[_, i] + s_prob[_, j, i]
                # print("a_star[_, i]")
            # print("a_star")
            # print(a_star)
            # print("s_prob")
            # print(s_prob[_,:,:])
            # print(np.sum(s_prob[_,:,:]))
            denom_a = np.sum(s_prob[_,:,:])
            # denomenator_a = [s_prob[_, t_x, i_x] for t_x in range(len(q_seq) - 1) for i_x in range(len(states))]
            # denomenator_a = sum(denomenator_a)
            # print("denomenator_a: " + str(denomenator_a))
            # for t in range(len(q_seq) - 1):
            #     for r in range(len(states)):
            #         print(_, t, r, s_prob[_, t, r])
            #         denom_a += s_prob[_, t, r]
            #         print("denom_a")
            #         print(denom_a)
            # denom_a = denom_a.sum()
            # print("~~~~~~~~~~~")
            # print(denom_a)

            if denom_a > 0:
                a_star[_, i] = a_star[_, i] / denom_a
            else:
                denom_a = a_star[_, i] = 0

        for i in range(len(alphabet)):
            idxs = [idx for idx, val in enumerate(q_seq) if val == alphabet[i]]
            numer_b = sum(g_prob[_, idxs])
            denom_b = sum(g_prob[_, :])

            if denom_b == 0:
                b_star[_, i] = 0
            else:
                b_star[_, i] = numer_b / denom_b

    return a_star, b_star


def sigmoid(x):
    if x >= 0:
        z = np.exp(-x)
        return 1 / (1 + z)
    else:

        z = np.exp(x)
        return z / (1 + z)

def viterbi2(seq, states, start_matrix, trans_matrix, emit_matrix):
    trellis = np.zeros((len(states), len(seq)))
    pointers = np.zeros((len(states), len(seq)))

    for s in range(len(states)):
        trellis[s, 0] = sigmoid(start_matrix.iloc[s].loc['probs']) * sigmoid(emit_matrix.iloc[s].loc[seq[0]])
        print(trellis)
    for o in range(1, len(seq)):

        for s in range(len(states)):
            k = np.argmax(trellis[k, o-1] * sigmoid(trans_matrix.iloc[k].iloc[s]) * sigmoid(emit_matrix.iloc[s].loc[seq[o]]) for k in trellis)
            trellis[s, o] = trellis[k, o-1] * sigmoid(trans_matrix.iloc[k].iloc[s]) * sigmoid(emit_matrix.iloc[s].loc[seq[o]])
            # print(trellis)
            pointers[s, o] = k
    best_path = []
    print(pointers)
    k = np.argmax(trellis[k, len(seq)-1])
    print(k)
    print(trellis)
    for o in range(len(seq)-1, -1, -1):

        best_path += states[int(k)]
        # best_path.insert(int(0), states[k])
        k = pointers[int(k), o]
    return best_path


def viterbi(seq, trans_matrix, emit_matrix, start_matrix):
    current_prob = {}
    ptr = {state: [] for state in states}

    for state in states:
        # current_prob[state] = np.log(start_matrix.loc[state, 'probs'] * emit_matrix.loc[state, seq[0]])
        current_prob[state] = sigmoid(start_matrix.loc[state, 'probs'] * emit_matrix.loc[state, seq[0]])

        print(current_prob)
    for _ in range(1, len(seq)):
        prev_prob = current_prob
        current_prob = {}
        for state in states:
            # V, previous_state = max((np.log(prev_prob[prev_state] * trans_matrix.loc[prev_state, state] * emit_matrix.at[
            #     state, seq[_]]), prev_state) for prev_state in states)
            V, previous_state = max((sigmoid(prev_prob[prev_state] * trans_matrix.loc[prev_state, state] * emit_matrix.at[
                state, seq[_]]), prev_state) for prev_state in states)
            # print(V, previous_state)
            current_prob[state] = V
            ptr[state].append(previous_state)
    print(ptr)
    print(current_prob)
    V = -1
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

    print(seq[:500])
    print(ks[:500])
    print("")
    print(q_seq[:500])
    print(ks_q[:500])

    nl = '\n'
    performance = f"HMM_bw Performance: {nl} Performance vs Training Set: {nl} %of loops(_) matched: {training_match__} " \
                  f"{nl} %of sheets(e) matched: {training_match_e} {nl} %of helixes(h) matched: {training_match_h} {nl} " \
                  f"Total % matched: {training_match_perf} {nl} {nl}Performance vs Query Set: {nl} %of loops(_) matched:" \
                  f" {query_match__} {nl} %of sheets(e) matched: {query_match_e} {nl} %of helixes(h) matched: " \
                  f"{query_match_h}{nl} Total % matched: {query_match_perf}"

    return performance


def HMM_bw(iterations=5):
    print("start")
    seq, ks, total_len, trans_matrix, emit_matrix, start_matrix, end_matrix = import_training(training_file)
    q_seq, ks_q = import_query(query_file)
    print("import complete")
    # print("start")
    # print(trans_matrix)
    # print(emit_matrix)
    # print(start_matrix)
    # print(end_matrix)

    delta_prev = 0
    neg = False
    for iter in range(iterations):
        print("iteration: " + str(iter))
        print("starting forward procedure")
        forward_prob, end_totals = forward_procedure(seq, q_seq, start_matrix, end_matrix, emit_matrix, trans_matrix)
        print("starting backward procedure")
        backward_prob, start_totals = backward_procedure(seq, q_seq, start_matrix, end_matrix, emit_matrix, trans_matrix)
        print("calculating s_g probabilities")
        s_prob, g_prob = s_g_probabilities(seq, q_seq, forward_prob, backward_prob, end_totals, trans_matrix, emit_matrix)
        print("calculating a_b star")
        a_star, b_star = a_b_mxs(seq, q_seq, s_prob, g_prob)

        # with np.errstate(divide='ignore'):
        #     a_star = np.log(a_star)
        #     b_star = np.log(b_star)
        # a_star[np.isneginf(a_star)] = 0
        # b_star[np.isneginf(b_star)] = 0

        a_star = pd.DataFrame(a_star, index=states, columns=states)
        b_star = pd.DataFrame(b_star, index=states, columns=alphabet)
        trans_matrix = a_star
        emit_matrix = b_star

        print("end_totals: " + str(end_totals))
        print(trans_matrix)
        print(emit_matrix)

        if iter > 2:
            bw_delta = np.abs(prev_end_totals[:][0] - end_totals[:][0])
            print("Slope: " + str(bw_delta))
            delta_delta = bw_delta - delta_prev
            print("Slope2: " + str(delta_delta))
            delta_prev = bw_delta

            # if delta_delta < 0: neg = True
            # if neg == True and delta_delta > 0:
            #     break

        prev_end_totals = end_totals

    # bw_test_pred_state = viterbi(seq, a_star, b_star, start_matrix)
    # bw_train_pred_state = viterbi(q_seq, a_star, b_star, start_matrix)
    #
    # print(solution_check(bw_test_pred_state, ks, bw_train_pred_state, ks_q))

    test_pred_state = viterbi2(seq, states, start_matrix, trans_matrix, emit_matrix)
    train_pred_state = viterbi2(q_seq, states, start_matrix, trans_matrix, emit_matrix)

    print(solution_check(test_pred_state, ks, train_pred_state, ks_q))
    return


HMM_bw()
