//
// FYS3150 project 1
// problem 9
// implementing special algorithm
//

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

vec f_func(vec x){
	// source function in Poisson problem
	return 100.*exp(-10.*x);
}

vec special_algorithm(int nSteps){
	int n = nSteps - 1; // length of vectors in matrix problem
	int m = nSteps + 1; // length of complete solution
	
	double x0 = 0.; // boundaries
	double xN = 1.;
	
	double h = (xN - x0)/nSteps; //steplength
	
	int v0 = 0; // boundary values
	int vN = 0;
	
	vec v = zeros(n); // inner solution v(1),...,v(N-1)
	
	vec x = linspace(x0 + h, xN - h, n);
	vec g = f_func(x) * h*h;
	g(0) += v0;
	g(n-1) += vN;
	
	
	vec b_tilde = zeros(n);
	b_tilde(0) = 2.;
	vec g_tilde = zeros(n);
	g_tilde(0) = g(0);
	
	// forward substitution
	
	for (int i = 1; i < n; i++){
		b_tilde(i) = 2 - 1/b_tilde(i-1);
		g_tilde(i) = g(i) + g_tilde(i-1)/b_tilde(i-1);
	}
	
	// backward substitution
	
	v(n-1) = g_tilde(n-1)/b_tilde(n-1);
	
	for (int i = n-2; i>=0; i--){
		v(i) = (g_tilde(i) + v(i+1))/b_tilde(i);
	}
	
	return v;
}

int main(int argc, char **argv){
	int nSteps = stoi(argv[1]); // no. of steps
	
	vec v = special_algorithm(nSteps);
	double v0 = 0.;
	double vN = 0.;
	int n = nSteps - 1; 
	int m = nSteps + 1; 
	double x0 = 0.; 
	double xN = 1.;
	double h = (xN - x0)/nSteps;
	vec x = linspace(x0 + h, xN - h, n);
	
	bool printSol = false;
	
	if (printSol){
		
		// print solution
		
		int nDigits = 6; // no. of significant digits in output
		cout.precision(nDigits);
		scientific; // set output to scientific notation
		
		cout << "x_i" << "    " << "v_i\n";
		cout << x0 << "    " << v0 << endl;
		for (int i = 0; i < n; i++){
			cout << x(i) << "    " << v(i) << endl;
		}
		cout << xN << "    " << vN;
	}
	return 0;
}