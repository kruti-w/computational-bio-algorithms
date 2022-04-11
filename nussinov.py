#!/usr/bin/env python
# coding: utf-8

# In[36]:


import numpy as np

# Nussinov 
# looking for the maximum number of basepairs that a certain sequence can
# -- theoretically -- contain according to a simple scoring model.

# Compares a sequence against itself in a dynamic programming matrix

# 4 rules for scoring the structure at a given point:
# 1. Add i,j pair onto the best structure found for subsequence i+1, j-1
# 2. Add unpaired position i onto best structure for subsequence i+1, j
# 3. Add unpaired position j onto best structure for subsequent i, j-1
# 4. Combine two optimal substructures i,k and k+1, j

# Initialization
# M[i][i] = 0 for i=1 to len(seq)
# M[i][i-1] = 0 for i=2 to len(sequence)

#M[i][j] = max of the following:
#   M[i+1][j] (i-th residue hanging off by itself)
#   M[i][j-1] (j-th residue hanging off by itself)
#   M[i+1][j-1] + S(xi, xj) (i-th and j-th residue are paired;
#if xi = complement of x, then S(xi, ji) = 1; otherwise it is 0)
#   M[i][j] = MAX (M[i][k] + M[k+1][j]) (merging two substructures)

# sequence = 'GCCACUGCC'

# initialize(len(sequence))



# get the sequence from the file
with open('5.in', 'r') as f:
    input_str = f.readlines()[1:]
    sequence = input_str[0]

# initialize the matrix
def initialize(length):
    matrix = np.full((length, length), 0)
    # print(matrix)
    return matrix
    
#s(A,U) = s(G,U) = s(U,A) = s(U,G) = -2
#s(G,C) = s(C,G) = -3
# Penalty for unpaired bases in the loop is +0.1, for unpaired bases out 
# of the loop it is 0

def scoring(seq, matrix, i, j):
    if j-i > 3:
    
        list1 = list2 = []
        if (i+1) < len(seq):
            list1.append(matrix[i+1, j]) # thing below it
        # also i+1 is the boundary condition to look out for     
            list1.append((matrix[i+1, j-1] + base_pairs(seq, i, j))) # check for base pairs

            for x in range(i+1,j):
                list1.append((matrix[i][x] + matrix[x+1][j]))

        list1.append(matrix[i, j-1]) # thing to its left
        matrix[i][j] = min(list1)
        
    else:
        return 0

# What base/boundary cases could there be? Need to look out for those
def fill_matrix(seq, matrix):
    # parses through diagonals - top-left to bottom-right
    i = j = 0
    for i in range(0, len(seq)):
        for k in range(0, len(seq)-i):
            # matrix[k, i+k] = i
            scoring(seq, matrix, k, i+k)
            
            # call scoring() for each cell of the upper triangular matrix        
            
    #print(matrix)
    return
    

def traceback(seq, matrix):
    traceback_stack = []
    base_pair_list = []
    i = 0
    j = len(seq)-1 #start in the top right corner
    seq_copy = seq
    energy = matrix[i][j]
    #print(energy)
    #checking for the smallest number
    # if there is a tie, go for the diagonal
    # diagonal = there is a pair (two bases corresponding to the square I'm coming from)
    # moving left or down - means there is an unpaired base
    # not in a loop until you've hit the diagonal - every time you move up or down, add 0.1
    #while loop - repeat until stack is empty
    traceback_stack.append([0, len(seq)-1]) #initialize
    structure = ['.' for x in seq]
    while len(traceback_stack)>0:
        h = traceback_stack.pop(len(traceback_stack)-1)
        i = h[0]
        j = h[1]
        # print('flag')
        if i >= j:
            if len(base_pair_list) > 0:
                energy += 0.1
            continue
        elif (matrix[i][j] == matrix[i+1][j]):
            traceback_stack.append([i+1, j])
            if len(base_pair_list) > 0:
                energy += 0.1
            #print(0)
        elif (matrix[i][j] == matrix[i+1][j-1] + base_pairs(seq, i, j)):
            traceback_stack.append([i+1, j-1])
            base_pair_list.append([i,j])
            structure[i] = '('
            structure[j] = ')'
            
        elif(matrix[i][j] == matrix[i][j-1]):
            traceback_stack.append([i,j-1])
            if len(base_pair_list) > 0:
                energy += 0.1
        else:
            for k in range(i+1, j-1):
                if matrix[i][j] == matrix[i][k] + matrix[k+1][j]:
                    traceback_stack.append([k+1, j])
                    traceback_stack.append([i,k])
                    break
    structure = ''.join(structure)
    energy = round(energy,2)
    print (energy)
    print(structure)
    return 0
        
        
        

def base_pairs(seq, base1, base2):
    if (seq[base1] == 'A' and seq[base2] == 'U' or seq[base1] == 'U' and seq[base2] == 'A' or
       seq[base1] == 'U' and seq[base2] == 'G' or seq[base1] == 'G' and seq[base2] == 'U'):
        return -2
    elif (seq[base1] == 'G' and seq[base2] == 'C' or seq[base1] == 'C' and seq[base2] == 'G'):
        return -3
    
    return 0
    

# sequence = 'GGCACUGAA'

mtrx = initialize(len(sequence))

print('Q1. part a)')

print('\n', sequence)
fill_matrix(sequence, mtrx)
#print('!!!!')
traceback(sequence, mtrx)
#print(0)


# In[ ]:





# In[ ]:




