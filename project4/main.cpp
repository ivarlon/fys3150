# include <iostream>
# include "Lattice.hpp"
# include <vector>
# include <fstream>
# include <random>
# include "omp.h"

using namespace std;
using namespace arma;

void Monte_Carlo_cycle(
    Lattice& lattice, 
    vector<double>& exp_list, 
    int& E, 
    int& M, 
    mt19937& generator, 
    uniform_int_distribution<int>& pdf_int,
    uniform_real_distribution<double>& pdf_uniform)
{
    int i, j, delta_E_half;
    double r, boltzmann_factor;
    
    for (int n=0; n<lattice.N; n++)
    {
        // pick random spin
        i = pdf_int(generator);
        j = pdf_int(generator);
        // flip
        lattice.flip_spin(i,j);
        // delta_E = -2 * s_ij * sum(neighbours) = -8, -4, 0, 4, 8
        delta_E_half = -lattice.S(i,j)*lattice.sum_neighbours(i,j);
        // sum = -4,-2,0,2,4 ===> (sum + 4)/2 = 0,1,2,3,4
        boltzmann_factor = exp_list[delta_E_half/2 + 2];
        r = pdf_uniform(generator);
        //cout << "r " << r << " p " << boltzmann_factor << endl;
        if (r>boltzmann_factor) 
        {
            lattice.S(i,j) *= -1;
        }
        else 
        {
            // update observables
            E += 2*delta_E_half;
            M += 2*lattice.S(i,j);
        }


    }
}

int main(int argc, char* argv[])
{
    if (argc == 1)
    {
        cout << "You need to input:" << endl;
        cout << "lattice size L" << endl;
        cout << "no. of MC cycles" << endl;
        cout << "ordered lattice (bool = 0,1)" << endl;
        cout << "save_results (0 or 1)" << endl;
        cout << "temperature T" << endl;
        cout << "OR, instead," << endl;
        cout << "T_min, T_max and no. of steps between them" << endl;
        exit(1);
    }

    int L = stoi(argv[1]); // lattice size
    int n_cycles = stoi(argv[2]); // no. of MCMC iterations
    bool ordered = stoi(argv[3]);

    bool save_results = stoi(argv[4]); // sets whether to save results to file


    // from here on:
    // program will run using one T value if only one value has been input
    // OR
    // program will use parallelisation looping over desired T values

    
    if (argc == 6)
    {
        double T = stod(argv[5]);
        // setting up list containing Boltzmann factors
        vector<double> exp_list;
        for (int delta_E = -8; delta_E < 9; delta_E += 4)
        {
            exp_list.push_back(exp(-delta_E/T));
        }
        
        // creating spin lattice
        Lattice lattice = Lattice(L, ordered);

        ofstream outfile;
        if (save_results)
        {   // creating data file
            string filename = "output_L_" + to_string(L) + "_T_" + to_string(T);
            if (ordered) 
            {filename += "_ordered.csv";}
            else
            {filename += "_disordered.csv";}
            outfile.open(filename);
        }

        // creating RNG
        mt19937 generator;
        int seed = 0;
        generator.seed(seed);

        // creating relevant pdfs
        uniform_int_distribution<int> pdf_int(0,L-1);
        uniform_real_distribution<double> pdf_uniform(0.,1.);

        // running MCMC
        int cycle = 0;
        int E = lattice.total_energy();
        int M = lattice.magnetisation();
        outfile << "cycle," << "E," << "M";
        outfile << "\n" << cycle << "," << E << "," << M;
        for (cycle = 1; cycle <= n_cycles; cycle++)
        {
            //cout << i << " " << j << endl;
            Monte_Carlo_cycle(lattice, exp_list, E, M, generator, pdf_int, pdf_uniform);
            outfile << "\n" << cycle << "," << E << "," << M;
        }
    }
    else
    {
        double T_min = stod(argv[5]);
        double T_max = stod(argv[6]); // temperature
        int n_T = stoi(argv[7]);
        double delta_T;
        delta_T = (T_max - T_min)/n_T;
        // creating RNG
        mt19937 generator; // Mersienne-Twister RNG
        int seed = 0;
        generator.seed(seed);

        // creating relevant pdfs
        uniform_int_distribution<int> pdf_int(0,L-1);
        uniform_real_distribution<double> pdf_uniform(0.,1.);
        
        # pragma omp parallel // start parallel region
        {
            # pragma omp for
            for (int i_T = 0; i_T <= n_T; i_T++)
            {
                double T = T_min + i_T*delta_T;
                // setting up list containing Boltzmann factors
                vector<double> exp_list;
                for (int delta_E = -8; delta_E < 9; delta_E += 4)
                {
                    exp_list.push_back(exp(-delta_E/T));
                }
                
                // creating spin lattice
                Lattice lattice = Lattice(L, ordered);

                ofstream outfile;
                if (save_results)
                {   // creating data file
                    string filename = "output_L_" + to_string(L) + "_T_" + to_string(T);
                    if (ordered) 
                    {filename += "_ordered.csv";}
                    else
                    {filename += "_disordered.csv";}
                    outfile.open(filename);
                }

                // running MCMC
                int cycle = 0;
                int E = lattice.total_energy();
                int M = lattice.magnetisation();
                outfile << "cycle," << "E," << "M";
                outfile << "\n" << cycle << "," << E << "," << M;
                for (cycle = 1; cycle <= n_cycles; cycle++)
                {
                    //cout << i << " " << j << endl;
                    Monte_Carlo_cycle(lattice, exp_list, E, M, generator, pdf_int, pdf_uniform);
                    outfile << "\n" << cycle << "," << E << "," << M;
                }
                # pragma omp critical
                {
                    cout << "Thread no. " << omp_get_thread_num() << "/" << omp_get_num_threads() << " finished T = " << T << "." << endl;
                }
            }
        }
    }
    return 0;
}