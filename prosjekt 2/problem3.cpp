// FYS3150
// project 2
// problem 3
// function identifying largest off-diagonal element in a matrix

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

double max_offdiag_symmetric(const mat& A, int&k, int& l)
{
	int N = A.n_rows;
	
	double max_offdiag = 0.;
	
	// looping over rows and columns
	for (int i = 1; i < N; i++){
		for (int j = 0; j < i; j++){
			cout << A(i,j) << endl;
			if (abs(A(i,j)) > abs(max_offdiag)){
				max_offdiag = A(i,j);
				k = i;
				l = j;
				
			}
		}
	}
	
	return max_offdiag;
}

int main(){
	
	// test max_offdiag_symmetric on test matrix
	mat A = mat(4, 4, fill::eye);
	A(0,3) = 0.5;
	A(3,0) = A(0,3);
	A(1,2) = -0.7;
	A(2,1) = A(1,2);
	
	int k;
	int l;
	
	cout << "Max off-diag. element in A:\n" << max_offdiag_symmetric(A, k, l) << ", (k,l) = (" << k << "," << l << ")" << endl;
	
	return 0;
}