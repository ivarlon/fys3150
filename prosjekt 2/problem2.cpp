// FYS4150 project 2
// problem 2
// setting up NxN tridiagonal matrix A for N = 6, solves
// eigenval. eq. and checks agreement w/ anal. result
//

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

mat create_tridiagonal(int N, double a, double d, double e){
	// creates tridiag. matrix A
	// N elements on diagonal have value d
	// N-1 elements on subdiagonal have value a
	// N-1 elements on superdiagonal have value e
	
	mat A = mat(N, N, fill::eye)*d;
	
	A(0,1) = e; // initial superdiag. element
	
	// loop for rows i=1 to penultimate row i=N-2
	for (int i = 1; i < N - 1; i++){
		A(i, i-1) = a; // subdiagonal
		A(i, i+1) = e; // superdiagonal
	}
	
	A(N-1,N-2) = a; // final subdiag. element
	
	return A;
}

int main(int argc, char **argv){
	int N = 6; // size of matrix
	int nSteps = N + 1;
	double h = 1./nSteps;
	
	double a = -1./(h*h);
	double d = 2./(h*h);
	
	mat A = create_tridiagonal(N, a, d, a);
	
	vec eigval;
	mat eigvec;
	
	eig_sym(eigval, eigvec, A);
	
	// armadillo doesn't work for me so code is incomplete.
	// compare eigval and eigvec with analytical eig.values and eig.vectors
	
	return 0;
}