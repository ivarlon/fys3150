/*
Defines a Lattice class
parameters:
    L (int) - length of lattice
    ordered (bool) - whether to start from state of all spin-up or random orientations
*/

# include <iostream>
# include <armadillo>
# include "Lattice.hpp"
# include <random>

using namespace arma;

//constructor
Lattice::Lattice(int L_in, bool ordered)
{
    L = L_in; // lattice size
    N = L*L; // no. of spins
    if (ordered)
    {
        S = mat(L, L, fill::ones);
    }
    else
    {
        S = randi<mat>(L, L, distr_param(0,1));
        S = S*2;
        S = S - ones(L, L);

    }
}

// calculates total energy E
int Lattice::total_energy()
{
    int E = 0;
    for (int i=0; i<L; i++)
    {
        for (int j=0; j<L; j++)
        {
            E -= S(i,j)*(S(i,(j+1)%L) + S((i+1)%L,j));
        }
    }
    return E;
}

// calculates magnetisation M
int Lattice::magnetisation()
{
    int M = 0;
    for (int i=0; i<L; i++)
    {
        for (int j=0; j<L; j++)
        {
            M += S(i,j);
        }
    }
    return M;
}

// flips a spin
void Lattice::flip_spin(int& i, int& j)
{
    //int i = draw from distribution;
    //int j = draw from distribution;
    S(i,j) *= -1;
}

// sums neighbour spins
int Lattice::sum_neighbours(int& i, int& j)
{
    // armadillo is stupid and doesn't like -1 to index final element
    // so need to do some computations to get right indices
    int spin_sum = S(i,(j+1)%L) + S(i, (j-1+L)%L) + S((i+1)%L,j) + S((i-1+L)%L,j);
    return spin_sum;
}