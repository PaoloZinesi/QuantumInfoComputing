import os

# dimension of random matrices
size = 10000

# number of repetitions
N_rep = 10

# list of different methods
method_list = ['HS', 'DS']

for rep_ in range(N_rep):
    for method_ in method_list:
        os.system(f'./a.out {size} {method_}')
