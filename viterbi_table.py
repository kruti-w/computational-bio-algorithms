#!/usr/bin/env python
# coding: utf-8

# In[12]:


import copy
from math import log
import pprint
import pandas as pd

seq = 'GGCACTGAA'

emission_dict = {'H': {'A': 0.2, 'C': 0.3, 'G': 0.3, 'T': 0.2},
                 'L': {'A': 0.3, 'C': 0.2, 'G': 0.2, 'T': 0.3}}

transition_dict = {'0': {'H': 0.5, 'L': 0.5}, 'H': {'H': 0.6, 'L': 0.4}, 'L': {'H': 0.5, 'L': 0.5}}


# part b) the complete Viterbi table
# This is currently in normal Viterbi format, not log-space format
# log format
# part c would also use viterbi to find the most probable path
def viterbi(sequence, transition_dictionary, emission_dictionary):
    # for loop for each column
    flag = False
    sequence = '0' + sequence
    p = []
    empty_col = {'0': 0, 'H': 0, 'L': 0}
    for x in sequence:
        p.append(copy.deepcopy(empty_col))
    p[0]['0'] = 1
    path_table = copy.deepcopy(p)
    # print(p)

    for i, x in enumerate(sequence):
        if i == 0:
            continue

        for k, v in p[i].items():
            if k == '0':
                continue  # skip the initial state
            temp_path = []

            v_best = 0
            v_temp = 0
            best_prev_state = ''

            for prev_k, prev_v in p[i - 1].items():
                temp_path.append(prev_v * transition_dictionary[prev_k][k] * emission_dictionary[k][x])
                p[i][k] = max(temp_path)
                v_temp = temp_path[len(temp_path)-1]
                print(v_temp, v_best)

                if v_temp > v_best:
                    v_best = v_temp
                    best_prev_state = prev_k

                for temp in temp_path:
                    if temp_path.count(temp) > 1:
                        flag = True

            path_table[i][k] = best_prev_state

    print(*p)
    pprint.pprint(p)
    a = pd.DataFrame(p).transpose()
    print(a)

    if flag == True:
        print('YES')
    else:
        print('NO')

    most_probable_path = []
    # print("The most probable final state")
    final_column = p[len(sequence) - 1]
    if final_column['H'] > final_column['L']:
        last_state = 'H'
    else:
        last_state = 'L'
    # print(last_state)

    print("The most probable path is:")
    for i in range(len(sequence)):
        most_probable_path.append('0')
    most_probable_path[len(sequence) - 1] = last_state
    i = len(sequence) - 2
    pprint.pprint(path_table)
    while i > 0:
        print(most_probable_path)
        most_probable_path[i] = path_table[i + 1][most_probable_path[i + 1]]
        i = i - 1
    print(most_probable_path)

    prob_most_probable_path = max((p[len(p) - 1]).values())
    # print(prob_most_probable_path)


# part a)  write the probability of producing x given the model
def forward(sequence, transition_dictionary, emission_dictionary):
    sequence = '0' + sequence
    p = []
    empty_col = {'0': 0, 'H': 0, 'L': 0}
    for x in sequence:
        p.append(copy.deepcopy(empty_col))
    p[0]['0'] = 1
    forward_table = []
    forward_table = copy.deepcopy(p)

    for i, x in enumerate(sequence):
        if i == 0:
            continue

        for k, v in p[i - 1].items():
            if k == '0':
                continue

            forward = 0
            for prev_k in p[i - 1].keys():
                forward = forward + (forward_table[i - 1][prev_k] * transition_dictionary[prev_k][k])
            forward_table[i][k] = forward * emission_dictionary[k][x]

    prob_of_obs = sum(forward_table[len(sequence) - 1].values())
    print(prob_of_obs)


def backward():
    pass


def baum_welch():
    pass


viterbi(seq, transition_dict, emission_dict)
# forward(seq, transition_dict, emission_dict)




