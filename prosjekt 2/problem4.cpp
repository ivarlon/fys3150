// FYS4150
// project 2
// problem 4
// Jacobi's rotation algorithm

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

mat create_tridiagonal(int N, double a, double d, double e)
{
	// creates tridiag. matrix A
	// N elements on diagonal have value d
	// N-1 elements on subdiagonal have value a
	// N-1 elements on superdiagonal have value e
	
	mat A = mat(N, N, fill::eye)*d;
	
	A(0,1) = e; // initial superdiag. element
	
	// loop for rows i=1 to penultimate row i=N-2
	for (int i = 1; i < N - 1; i++)
	{
		A(i, i-1) = a; // subdiagonal
		A(i, i+1) = e; // superdiagonal
	}
	
	A(N-1,N-2) = a; // final subdiag. element
	
	return A;
}

double max_offdiag_symmetric(arma::mat& A, int&k, int& l)
{
	// returns largest off-diagonal element in symmetric matrix A
	
	int N = A.n_rows;
	
	double max_offdiag = 0.;
	
	// looping over rows and columns
	for (int i = 1; i < N; i++)
	{
		for (int j = 0; j < i; j++)
		{
			if (fabs(A(i,j)) > fabs(max_offdiag))
			{
				max_offdiag = A(i,j);
				k = i;
				l = j;
				
			}
		}
	}
	
	return max_offdiag;
}


void jacobi_rotate(mat& A, mat& R, int k, int l)
{
	// transforms symm. matrix A and rotation matrix R
	
	double tau = (A(l,l) - A(k,k)) / (2.* A(k,l));
	double t;
	if (tau <= 0.)
	{
		t = -tau - sqrt(1. + tau*tau);
	}
	else
	{
		t = -tau + sqrt(1. + tau*tau);
	}
	
	double c = 1./sqrt(1. + t*t);
	double s = t*c;
	
	// transforming A matrix
	double A_kk = A(k,k);
	double A_ll = A(l,l);
	A(k,k) = A(k,k) * c*c - 2.*A(k,l) * c*s + A(l,l) * s*s;
	A(l,l) = A(l,l) * c*c + 2.*A(k,l) * c*s + A(k,k) * s*s;
	A(k,l) = 0.;//(c*c - s*s) * A(k,l) + c*s * (A_kk - A_ll);
	A(l,k) = A(k,l);
	int N = A.n_rows;
	for (int i = 0; i < N; i++)
	{
		if (i != k and i !=l)
		{	
			double A_ik = A(i,k); // saving values from step m to be used for updating at step m+1
			double A_il = A(i,l);
			A(i,k) = A_ik * c - A_il * s;
			A(k,i) = A(i,k);
			A(i,l) = A_il * c + A_ik * s;
			A(l,i) = A(i,l);
		}
	}
	
	// updating the rotation matrix R
	for (int i = 0; i < N; i++)
	{
		double R_ik = R(i,k);
		double R_il = R(i,l);
		R(i,k) = R_ik*c - R_il*s;
		R(i,l) = R_il*c + R_ik*s;
	}
	return;
}


void jacobi_eigensolver(arma::mat& A, double eps, vec& eigenvalues, mat& eigenvectors, const int maxiter, int& iterations, bool& converged){
	// Runs Jacobi rotation algorithm until off-diag. elements of A are smaller tolerance eps
	
	int N = A.n_rows;
	mat R = mat(N, N, fill::eye); // R1 = identity
	
	
	// loop until A is sufficiently diagonal or max. no. of iterations is reached
	int k; int l;
	double max_offdiag = max_offdiag_symmetric(A, k, l);
	while (fabs(max_offdiag) > eps and iterations <= maxiter){
		
		jacobi_rotate(A, R, k, l);
		max_offdiag = max_offdiag_symmetric(A, k, l);
		iterations ++;
		// cout << "it.no. " << iterations << endl;
		
	}
	
	converged = fabs(max_offdiag) < eps;	
	if (converged){
		for (int i=0; i<N; i++){
			eigenvalues(i) = A(i,i);
		}
		eigenvectors = R;
	}
	else{
		cout << "Algorithm didn't converge after " << maxiter <<" iterations." << endl;
		cout << "Max. off-diag. element: " << max_offdiag << " at (" << k << " " << l << ")" << endl;
		//cout << A << endl;
	}
	return;
}

int main(int argc, char **argv){
	// 1. argument of cmd line: size of matrix N
	// 2. argument: max. no. of iterations maxiter
	int N = stoi(argv[1]);
	int n_steps = N+1;
	double x0 = 0.; double xn = 1.;
	double h = (xn-x0)/n_steps;
	double a = -1./(h*h);
	double d = 2./(h*h);
	mat A = create_tridiagonal(N, a, d, a);
	if (N<=10){cout << "MATRIX A \n" << A << "\n";}
	
	double eps = 1e-7;
	vec eigenvalues = zeros(N);
	mat eigenvectors = zeros(N,N);
	int maxiter = stoi(argv[2]);
	int iterations = 0;
	bool converged;
	
	jacobi_eigensolver(A, eps, eigenvalues, eigenvectors, maxiter, iterations, converged);
	if (converged)
	{
		cout << "Algorithm converged after " << iterations << " iterations\n";
		cout << "Analytical vs. estimated eigenvalues:\n";
		uvec eigenvalorder = sort_index(eigenvalues);
		eigenvalues = eigenvalues(eigenvalorder);
		eigenvectors = eigenvectors.cols(eigenvalorder);
		for (int i = 0; i <N; i++)
		{
			double lmbda = d + 2*a*cos((i+1)*3.14159/n_steps); // analytic solution
			cout << lmbda <<  " vs. " << eigenvalues(i) << " giving a rel. error of " << 100.*fabs( (lmbda - eigenvalues(i) )/eigenvalues(i) ) << "%" << endl;
		}
		
		for (int i = 0; i<=2; i++)
		{
			cout << eigenvectors.col(i) << endl;	
		}
		if (N<=10){cout << A << endl;}
		int k; int l;
		//cout << max_offdiag_symmetric(A, k, l) << endl;
		// save eigenvectors
		eigenvectors.save("eigenvecs.bin", arma::raw_binary);
	}
	return 0;
}