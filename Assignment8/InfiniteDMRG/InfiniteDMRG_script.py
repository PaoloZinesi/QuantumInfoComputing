import os
import numpy as np
from datetime import datetime
import pytz

os.system(f"rm -f events.log")
os.system(f"rm -f errors.log")

# parameters
m_min, m_max = 2, 6
N_iter_max = int(10**8)
lambda_N = 51
lambda_min, lambda_max = 0.0, 3.0


# list of m to try
m_list = list(range(m_min, m_max+1))

# lists of lambdas
lambda_list = np.linspace(lambda_min, lambda_max, num=lambda_N)

# run for different m and for different lambdas
for m_ in m_list:
    for lambda_ in lambda_list:

        # write lambda into the file
        os.system(f"echo {lambda_} > input/lambda.dat")

        # start computation
        os.system(f"./InfiniteDMRG.out {N_iter_max} {m_} input/lambda.dat input/eps.dat 2>> events.log")

        nowtime = datetime.now(tz=pytz.timezone("Europe/Rome")).strftime("%d/%m/%Y %H:%M:%S")
        os.system(f"echo Finished computation of InfiniteDRMG.out m={m_}, maxiter={N_iter_max} at {nowtime} >> events.log")

