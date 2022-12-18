import os
import numpy as np

# parameters
K_max = 20
N_min, N_max = 2, 14
lambda_N = 100
lambda_min, lambda_max = 0.0, 3.0


# list of N to try
N_list = list(range(N_min, N_max+1))

# list of lambdas
lambda_list = np.linspace(lambda_min, lambda_max, num=lambda_N)

# run for different Ns and for different lambdas
for N_ in N_list:
    for lambda_ in lambda_list:

        # write lambda into the file
        os.system(f"echo {lambda_} > input/lambda.dat")

        # start computation
        os.system(f"./Ising1D.out {N_} {min(K_max,2**N_)} input/lambda.dat")

