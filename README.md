# EasyMPS
"A pedagogical realization of MPS method."

---------------------------------

Matrix product state (MPS) method is a series of powerful algorithms developed to solve a class of Hamiltonians with local interactions, based on the ansatz of “low entanglement”, or “area-law entanglement” precisely.



I would like to describe vMPS in a nutshell as follows:

(i) An arbitrary quantum many-body state (1-d, N sites) can be represented as a tensor on a set of certain bases. 

(ii) To solve the difficulty of the exponential growth of Hilbert space dimension with system size, we cut off the state tensor by performing Schmidt decomposition (equivalent to singular value decomposition, SVD) for any partition of the 1-d system and only keeping the data corresponding to the largest D Schmidt weights. 

(iii) In this way, we decompose the N-order state tensor as a contraction of N 3-order tensors corresponding to each site approximately, which gives us huge memory savings N*d*D^2 << d^N. Physically, we select a state most like the original one in the “low entanglement” subspace.



MPS is the 1-dimensional case of tensor network, which plays a central role in modern quantum physics and beyond.

Among a number of algorithms based on MPS, the variational MPS (vMPS) method is the most representative one, which is exactly equivalent with the famed density matrix renormalization group (DMRG) method. The vMPS method can be summarized as: optimizing the energy expectation of MPS by quadratic variation site-by-site (or 2-site) to converge to an approximate ground state.

Another typical algorithm is the infinite time-evolving block decimation (iTEBD), which evolves an initialized MPS iteratively by an Trotter-decomposed evolution operator.

---------------------------------

This repository is dedicated to provide an easy-to-understand implementation of MPS method. A number of ASCII sketches of tensors are presented in the comments to improve readability, such as

|'''|---              ---|'''|
|   |      |      |      |   |
|   |------H------H------|   |
|   |      |      |      |   |
|...|---              ---|...|
 
This repository will be updated continually in the near future (11/18/2020).
