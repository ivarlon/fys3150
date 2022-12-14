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
	
	//arma::mat A = arma::mat(N,N);
	//A.diag(-1).fill(a);
	//A.diag(0).fill(d);
	//A.diag(1).fill(e);
	
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
	
	// use Armadillo to find eigvalues and eigvectors
	vec eigval;
	mat eigvec;
	
	eig_sym(eigval, eigvec, A);
	
	// normalise
	eigvec = normalise(eigvec);
	
	// sort in right order
	uvec eigenvalorder = sort_index(eigval);
	eigval = eigval(eigenvalorder);
	eigvec = eigvec.cols(eigenvalorder);
	
	// calculate analytical solutions
	vec eigval_anal(N, fill::zeros);
	mat eigvec_anal(N,N, fill::zeros);
	
	for (int i = 1; i<=N; i++)
	{
		eigval_anal(i-1) = d + 2*a * cos(i*3.14159/(N+1));
		for (int jj=1; jj<=N; jj++)
		{
			eigvec_anal(jj-1, i-1) = sin(jj*i*3.14159/(N+1));
		}
		eigvec_anal.col(i-1) = normalise(eigvec_anal.col(i-1));
	}
	
	
	// print eigenvalues and eigenvectors
	
	cout << "Armadillo:" << endl;
	cout << eigval.t() << endl;
	cout << eigvec << endl;
	cout << "Analytical:" << endl;
	cout << eigval_anal.t() << endl;
	cout << eigvec_anal << endl;
	return 0;
}