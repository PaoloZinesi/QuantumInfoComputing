# %% [markdown]
# # Density Matrices Scalings - Many-Body Wavefunctions
# 
# Author: Paolo Zinesi

# %%
import numpy as np
from numpy import random
import pandas as pd
import os

import manybody_functions as mb

os.system("mkdir -p results")


# %% [markdown]
# ## Performance Testing

# %% [markdown]
# ### Separable and general wavefunction testing

# %%
rng = random.default_rng(123)
N_max = 8
D_max = 8
Nrep = 3
N_list = np.arange(start=2, stop=N_max+1, dtype=int)
D_list = np.arange(start=2, stop=D_max+1, step=2, dtype=int)

pd.DataFrame({"N":[], "D":[], "type":[], "dim_bytes":[], "time":[]})\
    .to_csv("results/wavefunc_perf.csv", header=True, index=False, mode='w')

# %%
for N_ in N_list:
    for D_ in D_list:
        for rep_ in range(Nrep):
            wf, dt, dim_bytes = mb.generate_separable_wf(N=N_, D=D_, rng=rng)

            pd.DataFrame({"N":[N_], "D":[D_], "type":["separable"], "dim_bytes":[dim_bytes], "time":[dt]})\
                .to_csv("results/wavefunc_perf.csv", header=False, index=False, mode='a')


            wf, dt, dim_bytes = mb.generate_general_wf(N=N_, D=D_, rng=rng)

            pd.DataFrame({"N":[N_], "D":[D_], "type":["general"], "dim_bytes":[dim_bytes], "time":[dt]})\
                .to_csv("results/wavefunc_perf.csv", header=False, index=False, mode='a')

wf = []

# %% [markdown]
# ### Density matrix testing

# %%
rng = random.default_rng(123)
N_max = 7
Nrep = 3
N_list = np.arange(start=2, stop=N_max+1, dtype=int)
D_list = np.array([2,4])

pd.DataFrame({"N":[], "D":[], "dim_bytes":[], "time":[]})\
    .to_csv("results/densitymatrix_perf.csv", header=True, index=False, mode='w')

# %%
for N_ in N_list:
    for D_ in D_list:
        for rep_ in range(Nrep):
            rho, dt, dim_bytes = mb.generate_general_densitymatrix(N=N_, D=D_, rng=rng)

            pd.DataFrame({"N":[N_], "D":[D_], "dim_bytes":[dim_bytes], "time":[dt]})\
                .to_csv("results/densitymatrix_perf.csv", header=False, index=False, mode='a')

rho = []

# %% [markdown]
# ### Partial trace testing

# %%
rng = random.default_rng(123)
N_max = 7
Nrep = 3
N_list = np.arange(start=2, stop=N_max+1, step=2, dtype=int)
D_list = np.array([2,4])

pd.DataFrame({"N":[], "D":[], "time":[]})\
    .to_csv("results/parttrace_perf.csv", header=True, index=False, mode='w')

# %%
for N_ in N_list:
    for D_ in D_list:
        for rep_ in range(Nrep):
            rho, _, _ = mb.generate_general_densitymatrix(N=N_, D=D_, rng=rng)

            rho_red, dt, _ = mb.reduced_densitymatrix(rho, N=N_, D=D_, traceout_indices=np.arange(N_//2))

            pd.DataFrame({"N":[N_], "D":[D_], "time":[dt]})\
                .to_csv("results/parttrace_perf.csv", header=False, index=False, mode='a')

rho = []
rho_red = []

# %% [markdown]
# ### Maximal entanglement test (D=2)

# %%
N_max = 14
Nrep = 3
N_list = np.arange(start=2, stop=N_max+1, dtype=int)

pd.DataFrame({"N":[], "reduced_entropy":[], "time":[]})\
    .to_csv("results/maxentang_D2_perf.csv", header=True, index=False, mode='w')

# %%
for N_ in N_list:
    for rep_ in range(Nrep):

        # create GHZ state
        wf = np.zeros(2**N_)
        wf[mb.combine_idxs([0]*N_, levels=range(N_))] = np.sqrt(0.5)
        wf[mb.combine_idxs([1]*N_, levels=range(N_))] = np.sqrt(0.5)

        # density matrix
        rho = np.outer(wf, wf)

        # partial trace over last subsystem
        rho_red, dt, _ = mb.reduced_densitymatrix(rho, N=N_, D=2, traceout_indices=np.random.choice(N_, size=1))

        # compute Von Neumann entropy of reduced matrix
        diag_rho_red = np.real(np.diag(rho_red))
        nonzero_mask = diag_rho_red > 0
        VN_entropy = -np.sum(diag_rho_red[nonzero_mask] * np.log2(diag_rho_red[nonzero_mask]))

        pd.DataFrame({"N":[N_], "reduced_entropy":[VN_entropy], "time":[dt]})\
            .to_csv("results/maxentang_D2_perf.csv", header=False, index=False, mode='a')

rho = []
rho_red = []

# %% [markdown]
# ### Scaling of partial trace on half subsystem (D=2, Nred=N/2)

# %%
rng = random.default_rng(123)
N_max = 14
Nrep = 3
N_list = np.arange(start=2, stop=N_max+1, step=2, dtype=int)

pd.DataFrame({"N":[], "time":[]})\
    .to_csv("results/parttrace_D2_perf.csv", header=True, index=False, mode='w')

# %%
for N_ in N_list:
    for rep_ in range(Nrep):
        rho, _, _ = mb.generate_general_densitymatrix(N=N_, D=2, rng=rng)

        rho_red, dt, _ = mb.reduced_densitymatrix(rho, N=N_, D=2, traceout_indices=np.arange(N_//2))

        pd.DataFrame({"N":[N_], "time":[dt]})\
            .to_csv("results/parttrace_D2_perf.csv", header=False, index=False, mode='a')

rho = []
rho_red = []
