#include <iostream>
#include <armadillo>
#include <vector>
#include "Particle.hpp"
#include "PenningTrap.hpp"

using namespace std;

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

void test_run( bool particleInteractions )
{
    return;
}

void initialise_particles(PenningTrap& trap, int nParticles)
{
    arma::arma_rng::set_seed(2022);
    int q = 1;
    double m = 20.;
    arma::vec r, v;
    for (int i = 0; i < nParticles; i++)
    {
        r = arma::vec(3).randn() * 0.1 * trap.d;
        v = arma::vec(3).randn() * 0.1 * trap.d;
        Particle p_i(q, m, r, v);
        trap.particles.push_back(p_i);
    }
}


int main()
{
    double B0 = 9.65e1; // u us^(-1) e^(-1)
    double V0 = 2.41e6; // u (um)^2 (us)^(-2) e^(-1)
    double d = 500.; // um
    bool particleInteractions = false;
    PenningTrap trap(B0, V0, d, 0., 0., particleInteractions);
    int q = 1; // e
    double m = 20.; // u
    double omega_z = sqrt( 2*q*V0/(m*d*d) ); // (us)^(-1) = MHz
    //std::cout << omega_z << std::endl;
    arma::vec r1("20 0 20");
    arma::vec v1("0 25 0");
    Particle p1(q, m, r1, v1);
    //trap.add_particle(p1);
    
    arma::vec r2("0 25 0");
    arma::vec v2("0 40 5");
    Particle p2(q, m, r2, v2);
    
    //trap.add_particle(p2);
    
    arma::vec r3("2 2 0");
    arma::vec v3("4 0 5");
    Particle p3(q, m, r3, v3);
    
    //trap.add_particle(p3);
    
    int nParticles = 100;
    initialise_particles(trap, nParticles);
    
    double T = 50.; //us
    int nSteps = 2000;
    double dt = T/nSteps;
    //std::vector<arma::vec> r_list;
    //r_list.push_back(p1.r);
    //std::vector<arma::vec> v_list;
    //v_list.push_back(p1.v);
    
    print_particle_positions(trap.particles);
    
    for (int i = 0; i < nSteps; i++)
    {
        trap.evolve_RK4(dt);
        print_particle_positions(trap.particles);
        //if (i%100 == 0){cout << i << trap.particles.size() << endl;}
    }
    
    return 0;
}