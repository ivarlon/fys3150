#include <iostream>
#include <armadillo>
#include <vector>
#include <fstream>
#include "Particle.hpp"
#include "PenningTrap.hpp"

using namespace std;

// physical properties of Ca+
int q = 1; // charge [e]
double m = 40.08; // mass [u]

/*
void print_particle_positions(std::vector<Particle>& particles)
{
    // prints position x_i y_i z_i and velocity vx_i vy_i vz_i for each particle i
    for (Particle& p_i : particles)
    {
        cout << p_i.r(0) << " " << p_i.r(1) << " " << p_i.r(2) << " "
        << p_i.v(0) << " " << p_i.v(1) << " " << p_i.v(2) << " ";
    }
    cout << endl;
}
*/

void save_data(PenningTrap& trap, std::ofstream& r_file, std::ofstream& v_file)
{
    // saves position and velocity data to input files
    for (Particle p_i : trap.particles)
    {
        for (int i = 0; i <= 2; i++)
        {
            r_file << p_i.r(i) << ",";
            v_file << p_i.v(i) << ",";
        }
    }
    r_file << "\n";
    v_file << "\n";
}

void test_run_one_particle(PenningTrap& trap)
{
    // add one particle to trap
    arma::vec r1("20 0 20");
    arma::vec v1("0 25 0");
    Particle p1(q, m, r1, v1);
    trap.add_particle(p1);
    
    std::ofstream r_file;
    r_file.open("position.csv");
    r_file.precision(10); // print sufficiently precise output
    
    std::ofstream v_file;
    v_file.open("velocity.csv");
    v_file.precision(10);
    
    save_data(trap, r_file, v_file);
    
    // integrate eqs of motion
    double T = 50.; //us
    int nSteps = 4000;
    double dt = T/nSteps;
    
    for (int i = 0; i < nSteps; i++)
    {
        trap.evolve_RK4(dt);
        save_data(trap, r_file, v_file);
    }
    
}

void simulate_two_particles(PenningTrap& trap)
{
    for (int particleInteractions = 0; particleInteractions <= 1; particleInteractions++)
    {
        // reset trap
        trap.reset_trap();
        
        // turn on/off Coulomb interactions
        trap.particleInteractions = particleInteractions;
        
        // adding two particles
        arma::vec r1("20 0 20");
        arma::vec v1("0 25 0");
        Particle p1(q, m, r1, v1);
        trap.add_particle(p1);
        
        arma::vec r2("25 25 0");
        arma::vec v2("0 40 5");
        Particle p2(q, m, r2, v2);
        trap.add_particle(p2);
        
        // complicated way of formatting filenames in the absence of a std::format method
        std::string r_filename, v_filename;
        if (particleInteractions)
        {
            r_filename = "position_partints_on.csv";
            v_filename = "velocity_partints_on.csv";
        }
        else
        {
            r_filename = "position_partints_off.csv";
            v_filename = "velocity_partints_off.csv";
        }
        
        // create datafiles
        std::ofstream r_file;
        r_file.open(r_filename);
        r_file.precision(10);
        
        std::ofstream v_file;
        v_file.open(v_filename);
        v_file.precision(10);
        
        save_data(trap, r_file, v_file);
        
        double T = 50.; //us
        int nSteps = 4000;
        double dt = T/nSteps;
        
        for (int i = 0; i < nSteps; i++)
        {
            trap.evolve_RK4(dt);
            save_data(trap, r_file, v_file);
        }
    }
}

void performance_tests(PenningTrap& trap)
{
    // check numerical accuracy for various timesteps for RK4 and FE
    std::vector<int> n_list;
    n_list.push_back(4000);
    n_list.push_back(8000);
    n_list.push_back(16000);
    n_list.push_back(32000);
    
    std::vector<std::string> methods;
    methods.push_back("RK4");
    methods.push_back("FE");
    
    for (std::string methodname : methods)
    {
        for (int n_i : n_list)
        {
            // reset trap
            trap.reset_trap();
            
            // adding particle
            arma::vec r1("20 0 20");
            arma::vec v1("0 25 0");
            Particle p1(q, m, r1, v1);
            trap.add_particle(p1);
            
            // create datafiles
            std::string r_filename, v_filename;
            r_filename = "position_" + methodname + "_" + std::to_string(n_i) + ".csv";
            v_filename = "velocity_" + methodname + "_" + std::to_string(n_i) + ".csv";
            std::ofstream r_file;
            r_file.open(r_filename);
            r_file.precision(14); // need fairly high precision esp. for RK4
            std::ofstream v_file;
            v_file.open(v_filename);
            v_file.precision(14);
            
            save_data(trap, r_file, v_file);
            
            double T = 50.; //us
            double dt = T/n_i;
            
            for (int i = 0; i < n_i; i++)
            {
                if (methodname=="RK4")
                {trap.evolve_RK4(dt);}
                else
                {trap.evolve_forward_Euler(dt);}
                save_data(trap, r_file, v_file);
            }
            
        }
        
        
    }
}

void initialise_particles(PenningTrap& trap, int nParticles, bool notRandom)
{
    // initialises particle pos. and vel. from normal distribution
    // if notRandom, pos. and vel. are not random X D
    if (notRandom)
    {arma::arma_rng::set_seed(2022);}
    
    arma::vec r, v;
    for (int i = 0; i < nParticles; i++)
    {
        r = arma::vec(3).randn() * 0.1 * trap.d;
        v = arma::vec(3).randn() * 0.1 * trap.d;
        Particle p_i(q, m, r, v);
        trap.particles.push_back(p_i);
    }
}

void resonance_exploration(PenningTrap& trap)
{
    // runs simulation with a periodically varying el. field
    
    // turn off particle interactions for faster performance
    trap.particleInteractions = false;
    
    // create list of amplitudes (using string most convenient for later creating filename string)
    std::vector<std::string> f_list;
    f_list.push_back("0.1"); f_list.push_back("0.4"); f_list.push_back("0.7");
    
    // create array of frequencies
    int nFreq = 116;
    arma::vec omega_array = arma::linspace(0.2,2.5, nFreq);
    
    int nParticles = 100;
    
    // set integration parameters
    double T = 500.; // us
    int nSteps = 25000;
    double dt = T/nSteps;
    
    for (std::string f : f_list)
    {
        cout << "f = " << f << endl;
        
        // creating a data file
        std::string filename;
        filename = "particles_remaining_" + f + ".csv";
        std::ofstream outfile;
        outfile.open(filename);
        
        // looping over frequencies
        for (int j = 0; j < nFreq; j++)
        {
            double omega = omega_array(j);
            
            // print loop progress
            if (j%10 == 0)
            {cout << 100.0*j / nFreq << "%" << endl;}
            
            // reset and reinitialise trap
            trap.reset_trap();
            initialise_particles(trap, nParticles, true);
            trap.f = std::stod(f);
            trap.omega = omega;
            
            // integration loop
            for (int i = 0; i < nSteps; i++)
            {
                trap.evolve_RK4(dt);
            }
            outfile << trap.particles_remaining() << "\n";
            
        }
    }
}

void resonance_finestructure(PenningTrap& trap)
{
    // looking for resonances with better resolution in frequency
    // with and without Coulomb interactions
    
    // create array of frequencies
    int nFreq = 151;
    arma::vec omega_array = arma::linspace(1.3,1.45, nFreq);
    
    int nParticles = 100;
    
    double f = 0.1;
    trap.f = f;
    
    // set integration parameters
    double T = 500.; // us
    int nSteps = 20000;
    double dt = T/nSteps;
    for (int particleInteractions = 0; particleInteractions < 2; particleInteractions++)
    {
        // toggle particle interactions
        trap.particleInteractions = particleInteractions;
        
        // creating a data file
        std::string filename;
        
        if (particleInteractions==0)
        {filename = "particles_remaining_finestruct_partint_off.csv";}
        else
        {filename = "particles_remaining_finestruct_partint_on.csv";}
        
        std::ofstream outfile;
        outfile.open(filename);
        
        // looping over frequencies
        for (int j = 0; j < nFreq; j++)
        {
            double omega = omega_array(j);
            
            // print loop progress
            if (j%5 == 0)
            {cout << 100.0*j / nFreq << "%" << endl;}
            
            // reset and reinitialise trap
            trap.reset_trap();
            initialise_particles(trap, nParticles, true);
            trap.omega = omega;
            
            // integration loop
            for (int i = 0; i < nSteps; i++)
            {
                trap.evolve_RK4(dt);
            }
            outfile << trap.particles_remaining() << "\n";
        }
    }
}

int main(int argc, char **argv)
{
    // Initialise trap
    double B0 = 9.65e1; // mag. field strength [u us^(-1) e^(-1)]
    double V0 = 2.41e6; // el. potential [u (um)^2 (us)^(-2) e^(-1)]
    double d = 500.; // char. scale [um]
    bool particleInteractions = true;
    PenningTrap trap(B0, V0, d, particleInteractions);
    
    /* 
    Uncomment line below to run desired function
    */
    
    
    test_run_one_particle(trap);
    //simulate_two_particles(trap);
    //performance_tests(trap);
    //resonance_exploration(trap);
    //resonance_finestructure(trap);
    
    return 0;
}