# quantum_solver_ed

Cpp library for exact diagonalization of interacting quantum many-body systems.\
Go to "examples" to see basic routines

**What it can do**:\
*/ Constructing a Hilbert space\
*/ Using symmetries to get irreps (block-diagonalize the problem)\
*/ Creating generic operators (e.g. a Hamiltonian) acting on the Hilbert space\
*/ Getting a matrix (sparse or dense) representation of a generic operator\
*/ On-fly calculation of a few lowest-energy eigenstates of a Hermitian operator using Lanczos algorithm\
*/ Calculating desired expectation values in an arbitrary state (e.g. ground-state of the Hamiltonian)

**What it can't do for instance (work in progress)**:\
*/ Can't do on-fly "fancier" algorithms such as Jacobi–Davidson, Arnoldi, etc

**Features**:\
*/ Allows to introduce custom spatial symmetries\
*/ Hamiltonian and symmetries can be read from an external file

**Requirements**:\
*/ Uses Armadillo (http://arma.sourceforge.net/) for basic linear algebra and internal storage of vectors, matrices, etc\
*/ Armadillo can also be nicely used as a wrapper around LAPACK, ARPACK, BLAS, etc
*/ Recommended packages to install: cmake (for armadillo), libopenblas-dev, liblapack-dev, libarpack2-dev, libsuperlu-dev
