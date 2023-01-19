import os
import numpy as np
from datetime import datetime
import pytz

os.system(f"rm -f events.log")
os.system(f"rm -f errors.log")

# parameters
N_min, N_max = 2, 7
N_iter_max = 10
lambda_N_long = 51
lambda_N_short = 11
lambda_min, lambda_max = 0.0, 3.0


# list of N to try
N_list = list(range(N_min, N_max+1))

# lists of lambdas
lambda_list_long = np.linspace(lambda_min, lambda_max, num=lambda_N_long)
lambda_list_short = np.linspace(lambda_min, lambda_max, num=lambda_N_short)

# run for different Ns and for different lambdas
for N_ in N_list:

    if N_ <= 6:
        lambda_list = lambda_list_long
    else:
        lambda_list = lambda_list_short

    for lambda_ in lambda_list:

        # write lambda into the file
        os.system(f"echo {lambda_} > input/lambda.dat")

        # start computation
        os.system(f"./RealSpaceRG.out {N_} {N_iter_max} input/lambda.dat input/eps.dat 2>> events.log")

        nowtime = datetime.now(tz=pytz.timezone("Europe/Rome")).strftime("%d/%m/%Y %H:%M:%S")
        os.system(f"echo Finished computation of RealSpaceRG.out N={N_}, maxiter={N_iter_max} at {nowtime} >> events.log")

