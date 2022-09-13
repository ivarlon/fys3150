//
// FYS3150 project 1
// problem 10
// timing general vs. special algorithm
//

#include <iostream>
#include <armadillo>
#include <chrono>

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

vec general_algorithm(int nSteps){
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
	return v;	
}

int main(int argc, char **argv){
	int nSteps = stoi(argv[1]); // no. of steps
	
	int nSample = 5;
	
	// time general algorithm
	vec timesGeneral = zeros(nSample);
	for (int i=0; i<nSample; i++){
		auto t1 = std::chrono::high_resolution_clock::now();
		vec v = general_algorithm(nSteps);
		auto t2 = std::chrono::high_resolution_clock::now();
		double duration_seconds = std::chrono::duration<double>(t2-t1).count();
		timesGeneral(i) = duration_seconds;
	}
	double meanGeneral = mean(timesGeneral);
	double sdGeneral = stddev(timesGeneral);
	cout << "General algorithm: " << meanGeneral << " +/- " << sdGeneral << " s" <<endl;
	
	// time special algorithm
	vec timesSpecial = zeros(nSample);
	for (int i=0; i<nSample; i++){
		auto t1 = std::chrono::high_resolution_clock::now();
		vec v = special_algorithm(nSteps);
		auto t2 = std::chrono::high_resolution_clock::now();
		double duration_seconds = std::chrono::duration<double>(t2-t1).count();
		timesSpecial(i) = duration_seconds;
	}
	double meanSpecial = mean(timesSpecial);
	double sdSpecial = stddev(timesSpecial);
	cout << "Special algorithm: " << meanSpecial << " +/- " << sdSpecial << " s" <<endl;
	
	return 0;
}