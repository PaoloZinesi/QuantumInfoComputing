# %% [markdown]
# # Density Matrices Scalings - Many-Body Wavefunctions
# 
# Author: Paolo Zinesi

# %%
import numpy as np
from numpy import random
import time
import pandas as pd

# %% [markdown]
# ## Functions definitions

# %%
def generate_separable_wf(N=2, D=2, rng=random.default_rng()):
    """
    Function to create a separable wavefunction considering D-dimensional Hilbert spaces and N subsystems.
    The wavefunction is initialized considering a normal-distributed variable with mu=0 and sigma=1 for both the
    real and the imaginary parts of each element, EXCEPT for the first element which is real and it fixes the normalization.
    
    inputs:
    - N [integer]: number of subsystems to consider
    - D [integer]: dimension of each subsystem
    - rng [random.Generator]: random generator to produce random numbers (useful to get reproducible results)
    
    outputs:
    - WF [np.ndarray]: numpy.array with shape (N,D,) that stores the wavefunction elements
    - time [float]: time (in seconds) necessary to create wavefunction
    - dim_bytes [integer]: dimension (in bytes) necessary to store the wavefunction
    """

    if not(np.issubdtype(np.array(N).dtype, np.integer)) or N < 1:
        print("N not valid")
        return

    if not(np.issubdtype(np.array(D).dtype, np.integer)) or D < 1:
        print("D not valid")
        return
    
    ti = time.perf_counter()
    
    # the first term is real (fixes global phase)
    WF = rng.standard_normal(size=(N,D,)) + 1.0j * rng.standard_normal(size=(N,D,))
    WF[0,0] = np.real(WF[0,0])
    WF /= np.sqrt(np.sum(np.abs(WF)**2))
    
    tf = time.perf_counter()
    dim_bytes = WF.nbytes

    return WF, tf-ti, dim_bytes

# %%
def generate_general_wf(N=2, D=2, rng=random.default_rng()):
    """
    Function to create the most general wavefunction considering D-dimensional Hilbert spaces and N subsystems.
    The wavefunction is initialized considering a normal-distributed variable with mu=0 and sigma=1 for both the
    real and the imaginary parts of each element, EXCEPT for the first element which is real and it fixes the normalization.
    
    inputs:
    - N [integer]: number of subsystems to consider
    - D [integer]: dimension of each subsystem
    - rng [random.Generator]: random generator to produce random numbers (useful to get reproducible results)
    
    outputs:
    - WF [np.ndarray]: numpy.array with shape (D**N,) that stores the wavefunction elements
    - time [float]: time (in seconds) necessary to create wavefunction
    - dim_bytes [integer]: dimension (in bytes) necessary to store the wavefunction
    """


    if not(np.issubdtype(np.array(N).dtype, np.integer)) or N < 1:
        print("N not valid")
        return

    if not(np.issubdtype(np.array(D).dtype, np.integer)) or D < 1:
        print("D not valid")
        return
    
    ti = time.perf_counter()

    # the first term is real (fixes global phase)
    WF = rng.standard_normal(size=(D**N,)) + 1.0j * rng.standard_normal(size=(D**N,))
    WF[0] = np.real(WF[0])
    WF /= np.sqrt(np.sum(np.abs(WF)**2))
    
    tf = time.perf_counter()
    dim_bytes = WF.nbytes

    return WF, tf-ti, dim_bytes

# %%
def generate_general_densitymatrix(N=2, D=2, rng=random.default_rng()):
    """
    Function to create the most general density matrix considering D-dimensional Hilbert spaces and N subsystems.
    The wavefunction is initialized from the "generate_general_wf" function, which generate normal-distributed variables with 
    mu=0 and sigma=1 for both the real and the imaginary parts of each element, EXCEPT for the first element.
    The density matrix is generated by taking the outer product of the generated wf.
    
    inputs:
    - N [integer]: number of subsystems to consider
    - D [integer]: dimension of each subsystem
    - rng [random.Generator]: random generator to produce random numbers (useful to get reproducible results)
    
    outputs:
    - rho [np.ndarray]: numpy.array with shape (D**N,D**N,) that stores the density matrix elements
    - time [float]: time (in seconds) necessary to create the density matrix
    - dim_rho_bytes [integer]: dimension (in bytes) necessary to store the density matrix
    """


    if not(np.issubdtype(np.array(N).dtype, np.integer)) or N < 1:
        print("N not valid")
        return

    if not(np.issubdtype(np.array(D).dtype, np.integer)) or D < 1:
        print("D not valid")
        return
    
    
    WF, time_wf, _ = generate_general_wf(N=N, D=D, rng=rng)

    ti = time.perf_counter()

    # rho = WF @ adj(WF)
    rho = np.outer(np.conj(WF), WF)

    tf = time.perf_counter()
    dim_rho_bytes = rho.nbytes

    return rho, time_wf+(tf-ti), dim_rho_bytes
    

# %%
def combine_idxs(multi_idx, levels, multi_idx2=[], levels2=[], D=2):
    """
    Function to compact a multi-index with D possible symbols. Optionally, the indices from two different
    arrays can be combined together if the two arrays of indices are disjoint. This function can be considered
    the inverse of 'vectorize_idx'.
    
    inputs:
    - multi_idx [list of integers]: list of indices that form a multi-index
    - levels [list of integers]: levels of each index (must be of the same size of 'multi_idx')
    - multi_idx2 [list of integers]: optional list of indices to be combined with 'multi_idx'
    - levels2 [list of integers]: levels of multi_idx2 (must be of disjoint w.r.t. 'levels')
    - D [integer]: dimension of each subsystem (or equivalently number of symbols used in the multi-index)
    
    output:
    - index [int]: compacted index
    """

    if len(np.intersect1d(np.array(levels), np.array(levels2))) > 0:
        print("The levels of the two indices are not disjoint")
        return

    if np.any(np.array(multi_idx) >= D) or np.any(np.array(multi_idx2) >= D):
        print("Some multi-indices are outside the range of D")
        return

    # total length
    N = len(multi_idx) + len(multi_idx2)

    # compacted index
    index = int(np.sum(np.array(multi_idx)*(D**(N-1-np.array(levels))))) + \
            int(np.sum(np.array(multi_idx2)*(D**(N-1-np.array(levels2)))))

    return index

# %%
def vectorize_idx(idx, N=1, D=2):
    """
    Function to transform a single index into a multi-index with D possible symbols and of length N.
    This function can be considered the inverse of 'combine_idxs'.
    
    inputs:
    - idx [integer]: single index to be transformed into a multi-index
    - N [integer]: number of subsystems (or equivalently number of symbols to be used to form the multi-index)
    - D [integer]: dimension of each subsystem (or equivalently number of symbols used in the multi-index)
    
    output:
    - multi_idx [list of integers]: list of indices that form a multi-index
    """

    if idx >= D**N or idx < 0:
        print(f"The index is outside the range [0, {D**N-1}]")
        return

    digits = []
    idx_ = idx
    for i in range(N):
        digits.append(int(idx_ % D))
        idx_ //= D
    
    return np.array(digits[::-1])

# %%
def reduced_densitymatrix(rho, N=2, D=2, traceout_indices=[0]):
    """
    Function to compute the subsystem density matrix given the density matrix of the total system.
    The total system is defined on a D-dimensional Hilbert spaces and has N subsystems. The indices on which to trace out the 
    total system density matrix are given in input.
    If input dimensions do not match (for example shape of rho is not (D**N,D**N)), an error is raised.
    
    inputs:
    - rho [np.ndarray]: numpy.array with shape (D**N,D**N,) that stores the density matrix elements
    - N [integer]: number of subsystems to consider. Must be even
    - D [integer]: dimension of each subsystem
    - traceout_indices [np.ndarray]: array of indices to trace out when computing the reduced density matrix
    
    outputs:
    - rho_reduced [np.ndarray]: numpy.array with shape (D**(N/2),D**(N/2),) that stores the reduced density matrix elements
    - time [float]: time (in seconds) necessary to compute the reduced density matrix
    - dim_rhored_bytes [integer]: dimension (in bytes) necessary to store the reduced density matrix
    """

    if not(np.issubdtype(np.array(N).dtype, np.integer)) or N < 1:
        print("N not valid")
        return

    if not(np.issubdtype(np.array(D).dtype, np.integer)) or D < 1:
        print("D not valid")
        return

    if type(rho) != np.ndarray:
        print("rho is not a valid numpy.ndarray")
        return

    if rho.shape[0] != rho.shape[1] or rho.ndim != 2:
        print("rho is not a square matrix")
        return

    if rho.shape[0] != D**N:
        print("rho dimensions and input dimensions do not match")
        return
        
    
    # number of indices to trace out
    Nred = len(traceout_indices)

    # indices to keep
    keep_indices = [idx for idx in range(N) if idx not in traceout_indices]
    Nkeep = len(keep_indices)
    

    ti = time.perf_counter()
    rho_reduced = np.zeros(shape=(D**Nkeep,D**Nkeep), dtype=np.complex128)


    # loop over all the indices of the reduced density matrix
    for row in range(D**Nkeep):
        vec_row = vectorize_idx(row, N=Nkeep, D=D)

        for col in range(row, D**Nkeep):
            vec_col = vectorize_idx(col, N=Nkeep, D=D)

            # loop over all the indices to trace out, find the correct position in the
            # total density matrix rho and add element to the reduced density matrix
            for kk in range(D**Nred):
                vec_kk = vectorize_idx(kk, N=Nred, D=D)
                
                row_kk = combine_idxs(vec_kk, traceout_indices, vec_row, keep_indices, D=D)
                col_kk = combine_idxs(vec_kk, traceout_indices, vec_col, keep_indices, D=D)
                
                rho_reduced[row, col] += rho[row_kk, col_kk]

    # the lower triangular matrix is simply the adjoint of the upper triangular matrix
    rho_reduced = rho_reduced + np.conj(rho_reduced.T) - np.diag(np.diag(rho_reduced))
    
    tf = time.perf_counter()
    dim_rhored_bytes = rho_reduced.nbytes

    return rho_reduced, tf-ti, dim_rhored_bytes