// Penning trap class

#include <iostream>
#include <vector>
#include <armadillo>
#include "Particle.hpp"
#include "PenningTrap.hpp"

using namespace std;

// Constructor for a constant applied potential
PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in, bool particleInteractions_in)
{
    B0 = B0_in;
    V0 = V0_in;
    d = d_in;
    V0_d2 = V0/(d*d);
    t = 0.;
    f = 0.;
    omega = 0.;
    particleInteractions = particleInteractions_in;
}

// Constructor for a time varying applied potential
PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in, 
    double f_in, double omega_in, bool particleInteractions_in)
{
    B0 = B0_in;
    V0 = V0_in;
    d = d_in;
    V0_d2 = V0/(d*d);
    t = 0.;
    f = f_in;
    omega = omega_in;
    particleInteractions = particleInteractions_in;
}

// Add a particle to the trap
void PenningTrap::add_particle(Particle p_in)
{
    particles.push_back(p_in);
}

// External el. field at point r
arma::vec PenningTrap::external_E_field(arma::vec r)
{
    // E = -V0(1 + fcos(wt))/2d^2 * (x, y, -2z)
    return V0_d2 * (1. + f*cos(omega*t) ) * arma::vec({r(0), r(1), -2.*r(2)});
}

// External mag. field at point r
arma::vec PenningTrap::external_B_field(arma::vec r)
{
    // B = (0,0,B0)
    return {0., 0., B0};
}

// Force on particle i from particle j
arma::vec PenningTrap::force_particle(int i, int j)
{
    double k_e = 1.38935333e5;
    int qi = particles[i].q;
    int qj = particles[j].q;
    arma::vec ri = particles[i].r;
    arma::vec rj = particles[j].r;
    return k_e * qi * qj * normalise(ri - rj) / dot(ri - rj, ri - rj);
}


// Force on particle i from the external fields
arma::vec PenningTrap::total_force_external(int i)
{
    int qi = particles[i].q;
    arma::vec ri = particles[i].r;
    
    // Check if particle is outside trap
    if (arma::norm(ri) > d)
    {
        return arma::vec(3, arma::fill::zeros);
    }
    else
    {
        arma::vec vi = particles[i].v;
        arma::vec E = external_E_field(ri);
        arma::vec B = external_B_field(ri);
        return qi * ( E + arma::cross(vi,B) );
    }
}

// Total force on particle i from other particles
arma::vec PenningTrap::total_force_particles(int i)
{
    arma::vec F(3, arma::fill::zeros);
    for (int j = 0; j < particles.size(); j++)
    {
        if (j != i)
        {
            F += force_particle(i, j);
        }
        else {continue;}
    }
    return F;
}

// Total force on particle i from both ext. fields and other particles
arma::vec PenningTrap::total_force(int i)
{
    if (particleInteractions)
    {
        return total_force_external(i) + total_force_particles(i);
    }
    else
    {
        return total_force_external(i);
    }
}

// Evolve the system one time step dt using Runge-Kutta 4th order
void PenningTrap::evolve_RK4(double dt)
{
    // create a copy of particles at time t_i
    vector<Particle> particles_copy = particles;
    
    // create a list to contain k1, k2, k3, k4
    vector<vector<arma::vec>> k_r;
    vector<vector<arma::vec>> k_v;
    
    // declare k_r(i), k_v(i)
    arma::vec k_r_(3);
    arma::vec k_v_(3);
    
    // calculating k_1 for each particle
    vector<arma::vec> k_r1;
    k_r.push_back(k_r1);
    vector<arma::vec> k_v1;
    k_v.push_back(k_v1);
    
    for (int i = 0; i < particles.size(); i ++)
    {
        k_v_ = total_force(i) / particles[i].m * dt; // dv = a dt = F/m dt
        k_v[0].push_back(k_v_);
        k_r_ = particles[i].v * dt;
        k_r[0].push_back(k_r_);
    }
    
    // updating particle pos. and vel. to half timestep using k_1
    for (int i = 0; i < particles.size(); i ++)
    {
        particles[i].r += k_r[0][i] * 0.5;
        particles[i].v += k_v[0][i] * 0.5;
    }
    t += 0.5*dt;
    
    // calculating k_2 (using y = y_i + 0.5 k1)
    vector<arma::vec> k_r2;
    k_r.push_back(k_r2);
    vector<arma::vec> k_v2;
    k_v.push_back(k_v2);
    
    for (int i = 0; i < particles.size(); i ++)
    {
        k_v_ = total_force(i) / particles[i].m * dt; // dv = a dt = F/m dt
        k_v[1].push_back(k_v_);
        k_r_ = particles[i].v * dt;
        k_r[1].push_back(k_r_);
    }
    
    // updating particle pos. and vel. to half timestep using k_2
    for (int i = 0; i < particles.size(); i ++)
    {
        particles[i].r = particles_copy[i].r + k_r[1][i] * 0.5;
        particles[i].v = particles_copy[i].v + k_v[1][i] * 0.5;
    }
    
    // calculating k_3 (using y = y_i + 0.5 k2)
    vector<arma::vec> k_r3;
    k_r.push_back(k_r3);
    vector<arma::vec> k_v3;
    k_v.push_back(k_v3);
    
    for (int i = 0; i < particles.size(); i ++)
    {
        k_v_ = total_force(i) / particles[i].m * dt; // dv = a dt = F/m dt
        k_v[2].push_back(k_v_);
        k_r_ = particles[i].v * dt;
        k_r[2].push_back(k_r_);
    }
    
    // updating particle pos. and vel. to half timestep using k_2
    for (int i = 0; i < particles.size(); i ++)
    {
        particles[i].r = particles_copy[i].r + k_r[2][i];
        particles[i].v = particles_copy[i].v + k_v[2][i];
    }
    t += 0.5*dt;
    
    // calculating k_4 (using y = y_i + k3)
    vector<arma::vec> k_r4;
    k_r.push_back(k_r4);
    vector<arma::vec> k_v4;
    k_v.push_back(k_v4);
    
    for (int i = 0; i < particles.size(); i ++)
    {
        k_v_ = total_force(i) / particles[i].m * dt; // dv = a dt = F/m dt
        k_v[3].push_back(k_v_);
        k_r_ = particles[i].v * dt;
        k_r[3].push_back(k_r_);
    }
    
    // final step: updating particle pos. and vel.
    for (int i = 0; i < particles.size(); i ++)
    {
        particles[i].r = particles_copy[i].r + (k_r[0][i] + 2.* k_r[1][i] + 2*k_r[2][i] + k_r[3][i]) / 6.;
        particles[i].v = particles_copy[i].v + (k_v[0][i] + 2.* k_v[1][i] + 2*k_v[2][i] + k_v[3][i]) / 6.;
    }
    
}

// Evolve the system one time step dt using Forward Euler
void PenningTrap::evolve_forward_Euler(double dt)
{
    for (int i = 0; i < particles.size(); i ++)
    {
        arma::vec a = total_force(i) / particles[i].m;
        particles[i].r += particles[i].v * dt;
        particles[i].v += a * dt;
    }
    
    t += dt;
}

int PenningTrap::particles_remaining()
{
    int nRemaining = 0;
    for (Particle& p_i : particles)
    {
        if (arma::norm(p_i.r) < d)
        {
            nRemaining += 1;
        }
    }
    
    return nRemaining;
}