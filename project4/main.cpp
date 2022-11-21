# include <iostream>
# include "Lattice.hpp"
# include <vector>
# include <fstream>
# include <random>

using namespace std;
using namespace arma;

void Monte_Carlo_cycle(
    Lattice& lattice, 
    vector<double>& exp_list, 
    int& E, 
    int& M, 
    bool& print_bool, 
    mt19937& generator, 
    uniform_int_distribution<int>& pdf_int,
    uniform_real_distribution<double>& pdf_uniform)
{
    int i, j;
    
    for (int n=0; n<lattice.N; n++)
    {
        i = pdf_int(generator);
        j = pdf_int(generator);
        lattice.flip_spin(i,j);
        int neighbour_sum = lattice.sum_neighbours(i,j);
        // sum = -4,-2,0,2,4 ===> (sum + 4)/2 = 0,1,2,3,4
        double boltzmann_factor = exp_list[(neighbour_sum + 4)/2];
        double r = pdf_uniform(generator);
        cout << "r " << r << " p " << boltzmann_factor << endl;
        if (r>boltzmann_factor) 
        {
            lattice.S(i,j) *= -1;
            print_bool = false;
        }
        else 
        {
            // update observables
            E -= 2*lattice.S(i,j)*neighbour_sum;
            M += 2*lattice.S(i,j);
            cout << "updated spin" << endl;
            print_bool = true;
        }


    }
}

int main(int argc, char* argv[])
{
    if (argc != 3)
    {
        cout << "You need to input lattice size L and no. of MC cycles..." << endl;
        exit(1);
    }

    int L = stoi(argv[1]);
    int n_cycles = stoi(argv[2]);

    double T = 1.; // temprecher

    // setting up list containing Boltzmann factors
    vector<double> exp_list;
    for (int delta_E = -8; delta_E < 9; delta_E += 4)
    {
        exp_list.push_back(exp(-delta_E/T));
    }
    
    bool ordered = false;
    // creating spin lattice
    Lattice lattice = Lattice(L, ordered);


    string filename = "output_L_" + to_string(L) + "_T_" + to_string(T);
    if (ordered) 
    {filename += "_ordered.csv";}
    else
    {filename += "_disordered.csv";}

    // creating RNG
    mt19937 generator;
    int seed = 202;
    generator.seed(seed);
    uniform_int_distribution<int> pdf_int(0,L-1);
    uniform_real_distribution<double> pdf_uniform(0.,1.);

    // running MCMC
    int cycle = 0;
    //outfile << "cycle" << "E" << "M";
    //outfile << cycle , lattice.total_energy, lattice.magnetisation;
    int E = lattice.total_energy();
    int M = lattice.magnetisation();
    bool print_bool = false;
    for (cycle = 1; cycle < n_cycles; cycle++)
    {
        //cout << i << " " << j << endl;
        Monte_Carlo_cycle(lattice, exp_list, E, M, print_bool, generator, pdf_int, pdf_uniform);
        if (print_bool) {cout << cycle << " " << E << " " << M << endl;}
    }
    cout << lattice.S << endl;
    return 0;
}