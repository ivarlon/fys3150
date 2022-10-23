// Particle class

#ifndef __Particle_hpp__
#define __Particle_hpp__

#include <armadillo>

class Particle
{
public:
    // Constructor
    Particle(int charge, double mass, arma::vec position, arma::vec velocity);
    
    int q; // el. charge
    double m; // mass
    arma::vec r; // 3D position
    arma::vec v; // 3D velocity
};

#endif