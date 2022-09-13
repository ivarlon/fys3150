//
// FYS3150 project 1
// problem 8
// solves system of equations Av = g
//

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

vec f_func(vec x){
	// source function in Poisson problem
	return 100.*exp(-10.*x);
}

vec u_func(vec x){
	return 1 - (1 - exp(-10))*x - exp(-10*x);
}

int main(int argc, char **argv){
	// inputs: 
	// 1. int nSteps
	// optional:
	// 2. arbitrary input MaxRelErrorOnly (for problem 8c when we only want max. rel. error)
	
	int nSteps = stoi(argv[1]); // no. of steps
	int n = nSteps - 1; // length of vectors in matrix problem
	int m = nSteps + 1; // length of complete solution
	
	double x0 = 0.; // boundaries
	double xN = 1.;
	
	double h = (xN - x0)/nSteps; //steplength
	
	int v0 = 0; // boundary values
	int vN = 0;
	
	vec v = zeros(n); // inner solution v(1),...,v(N-1)
	
	vec a = -1.0*ones(n - 1); // subdiagonal
	vec b = 2.0*ones(n); // diagonal
	vec c = -1.0*ones(n - 1); //superdiagonal
	vec x = linspace(x0 + h, xN - h, n);
	vec g = f_func(x) * h*h;
	g(0) += v0;
	g(n-1) += vN;
	
	
	vec b_tilde = zeros(n);
	b_tilde(0) = b(0);
	vec g_tilde = zeros(n);
	g_tilde(0) = g(0);
	
	
	// forward substitution
	
	for (int i = 1; i < n; i++){
		double a_over_b = a(i-1)/b_tilde(i-1);
		b_tilde(i) = b(i) - a_over_b * c(i-1);
		g_tilde(i) = g(i) - a_over_b * g_tilde(i-1);
	}
	
	
	// backward substitution
	
	v(n-1) = g_tilde(n-1)/b_tilde(n-1);
	
	for (int i = n-2; i>=0; i--){
		v(i) = (g_tilde(i) - c(i)*v(i+1))/b_tilde(i);
	}
	
	
	// find error between approx. and exact solution
	
	vec u = u_func(x);
	vec delta = zeros(n); // list of absolute error, |u_i - v_i|
	vec epsilon = zeros(n); // list of relative error, |(u_i - v_i)/u_i|
	for (int i=0; i<n; i++){
		delta(i) = abs(u(i) - v(i));
		epsilon(i) = abs(1 - v(i)/u(i));
	}
	
	// print solution
	
	int nDigits = 6; // no. of significant digits in output
	cout.precision(nDigits);
	scientific; // set output to scientific notation
	
	if (argc == 3){ // if we only want to print max rel. error
		
		cout << "nSteps = " << nSteps << "  max. rel. error = " << epsilon.max() << endl;
		
	} else {
		
		cout << "x_i" << "    " << "v_i" << "    delta_i    eps_i" ;
		cout << x0 << "    " << v0 << endl;
		for (int i = 0; i < n; i++){
			cout << x(i) << "    " << v(i) << "    " << delta(i) << "    " << epsilon(i) << endl;
		}
		cout << xN << "    " << vN;
	}
	
	return 0;
}