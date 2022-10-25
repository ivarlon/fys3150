// Penning trap class

#include <iostream>
#include <vector>
#include <armadillo>
#include "Particle.hpp"

class PenningTrap
{
    public:
    
    // Constructor for a constant applied potential
    PenningTrap(double B0_in, double V0_in, double d_in, bool particleInteractions_in);
    
    // Constructor for a time varying applied potential
    PenningTrap(double B0_in, double V0_in, double d_in, double f_in, double omega_in, bool particleInteractions_in);
    
    double B0, V0, d; // mag. field strength, el. potential and characteristic scale
    
    double V0_d2; // V0/d^2
    
    double t; // time
    
    double f, omega; // amplitude and ang. frequency of applied el. potential
    
    bool particleInteractions; // turn on (true) or off (false) part. interactions
    
    std::vector<Particle> particles; // list containing each Particle
    
    // Reset trap (t=0, no particles)
    void reset_trap();
    
    // Add a particle to the trap
    void add_particle(Particle p_in);
    
    // External el. field at point r=(x,y,z)
    arma::vec external_E_field(arma::vec r);
    
    // External mag. field at point r=(x,y,z)
    arma::vec external_B_field(arma::vec r);
    
    // Force on particle i from particle j
    arma::vec force_particle(int i, int j);
    
    // Force on particle i from the external fields
    arma::vec total_force_external(int i);
    
    // Total force on particle i from other particles
    arma::vec total_force_particles(int i);
    
    // Total force on particle i from both ext. fields and other particles
    arma::vec total_force(int i);
    
    // Evolve the system one time step dt using Runge-Kutta 4th order
    void evolve_RK4(double dt);
    
    // Evolve the system one time step dt using Forward Euler
    void evolve_forward_Euler(double dt);
    
    // Calculate no. of particles remaining in trap (r < d)
    int particles_remaining();

};