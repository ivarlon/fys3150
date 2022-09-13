#include <iostream>
#include <armadillo>

// 
// FYS4150 Project 1
// Problem 2
// writing x and u(x) to text file
//

using namespace std;

arma::vec u_func(arma::vec x){
	return 1 - (1 - exp(-10))*x - exp(-10*x);
}

int main(){
	int nSteps = 1000;
	arma::vec x = arma::linspace(0., 1., nSteps + 1); // input values
	arma::vec u = u_func(x); // output values
	
	int nDigits = 6; // no. of significant digits in output
	cout.precision(nDigits);
	scientific; // set output to scientific notation
	
	cout << "x    u(x) \n";
	for (int i = 0; i < nSteps + 1; i++){
		cout << x(i) << "   " << u(i) << "\n";
	}
	return 0;
}