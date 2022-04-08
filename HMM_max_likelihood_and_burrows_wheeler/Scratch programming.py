states = ('Healthy', 'Fever')

observations = ('normal', 'cold', 'dizzy')

start_probability = {'Healthy': 0.6, 'Fever': 0.4}

transition_probability = {
    'Healthy': {'Healthy': 0.7, 'Fever': 0.3},
    'Fever': {'Healthy': 0.4, 'Fever': 0.6},
}

emission_probability = {
    'Healthy': {'normal': 0.5, 'cold': 0.4, 'dizzy': 0.1},
    'Fever': {'normal': 0.1, 'cold': 0.3, 'dizzy': 0.6},
}


def Viterbit(obs, states, s_pro, t_pro, e_pro):
    path = {s: [] for s in states}  # init path: path[s] represents the path ends with s
    curr_pro = {}
    for s in states:
        curr_pro[s] = s_pro[s] * e_pro[s][obs[0]]
        print(curr_pro)
    for i in range(1, len(obs)):
        last_pro = curr_pro
        curr_pro = {}
        for curr_state in states:
            max_pro, last_sta = max(
                ((last_pro[last_state] * t_pro[last_state][curr_state] * e_pro[curr_state][obs[i]], last_state)
                 for last_state in states))
            print(max_pro, last_sta)
            curr_pro[curr_state] = max_pro
            path[curr_state].append(last_sta)

    # find the final largest probability
    max_pro = -1
    max_path = None
    for s in states:
        path[s].append(s)
        if curr_pro[s] > max_pro:
            max_path = path[s]
            max_pro = curr_pro[s]
    # print '%s: %s'%(curr_pro[s], path[s]) # different path and their probability
    return max_path


if __name__ == '__main__':
    obs = ['normal', 'cold', 'dizzy']
    print(Viterbit(obs, states, start_probability, transition_probability, emission_probability))

# def viterbi_algorithm(observations, states, start_p, trans_p, emit_p):
#     V = [{}]
#     for st in states:
#         V[0][st] = {"prob": start_p[st] * emit_p[st][observations[0]], "prev": None}
#
#     for t in range(1, len(observations)):
#         V.append({})
#         for st in states:
#             max_tr_prob = V[t - 1][states[0]]["prob"] * trans_p[states[0]][st]
#             prev_st_selected = states[0]
#             for prev_st in states[1:]:
#                 tr_prob = V[t - 1][prev_st]["prob"] * trans_p[prev_st][st]
#                 if tr_prob > max_tr_prob:
#                     max_tr_prob = tr_prob
#                     prev_st_selected = prev_st
#
#             max_prob = max_tr_prob * emit_p[st][observations[t]]
#             V[t][st] = {"prob": max_prob, "prev": prev_st_selected}
#
#
# for line in dptable(V):
#     print(line)
#
# opt = []
# max_prob = 0.0
# best_st = None
#
# for st, data in V[-1].items():
#     if data["prob"] > max_prob:
#         max_prob = data["prob"]
#         best_st = st
# opt.append(best_st)
# previous = best_st
#
# for t in range(len(V) - 2, -1, -1):
#     opt.insert(0, V[t + 1][previous]["prev"])
#     previous = V[t + 1][previous]["prev"]
#
# print("The steps of states are " + " ".join(opt) + " with highest probability of %s" % max_prob)
#
#
# def dptable(V):
#     yield " ".join(("%12d" % i) for i in range(len(V)))
#     for state in V[0]:
#         yield "%.7s: " % state + " ".join("%.7s" % ("%f" % v[state]["prob"]) for v in V)
