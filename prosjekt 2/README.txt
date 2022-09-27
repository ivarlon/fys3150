## problem2.cpp:
Incomplete code, couldn't get eig_sym to work. Only makes a tridiagonal matrix of size N=6


## problem3.cpp
Function that identifies largest off-diag. element in a symmetric matrix. Tests function on a matrix, giving -0.7, 1,2 as answer.


## problem4.cpp
Implements Jacobi-algorithm. Outputs estimated eigenvalues and eigenvectors.
Requires input N and maxiter. 
Tolerance set to epsilon=1e-6. 
Seems to give pretty bad results. 
Also the armadillo method sort_index gives the wrong ranking of the size of the eigenvalues, so the comparison of estimated and analytical eigenvalues is wrong.

## problem5a.py
Plots matrix size N versus no. of iterations.

## problem6.py
Plots eigenvectors corresponding to three lowest eigenvals, for nsteps=10, 100. Probably got the wrong eigenvectors, at least for nsteps=100.