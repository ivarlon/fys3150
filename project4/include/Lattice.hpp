/*
Defines a Lattice class
parameters:
    L (int) - length of lattice
    T (double) - temperature
*/

# include <armadillo>

class Lattice
{
    public:
    
    //constructor
    Lattice(int L_in, bool ordered);

    // lattice size
    int L, N;

    // matrix of spins
    arma::mat S;

    // calculates total energy E
    int total_energy();

    // calculates magnetisation M
    int magnetisation();

    // flips a spin
    void flip_spin(int& i, int& j);

    // sums neighbour spins
    int sum_neighbours(int& i, int& j);

};