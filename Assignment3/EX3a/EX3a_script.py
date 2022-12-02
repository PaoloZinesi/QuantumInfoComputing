import os
import math as m

# dimension grid parameters
N_min, N_max = 100, 12000
mult_factor = m.sqrt(2.0)

# number of repetitions (for statistics)
N_rep = 3

# list of dimensions to try
dim_list = [int(N_min * mult_factor**exp) for exp in range(m.floor(m.log(N_max/N_min, mult_factor))+1)]

# list of different methods
method_list = ['naive', 'opt', 'builtin']

# run for different times the desired program (a.out)
for dim_ in dim_list:
    for method_ in method_list:
        for rep_ in range(N_rep):
            os.system(f'./a.out {dim_} {dim_} {dim_} {dim_} {method_}')

