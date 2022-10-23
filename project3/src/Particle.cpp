// Particle class

#include "Particle.hpp"

// Constructor
Particle::Particle(int charge, double mass, arma::vec position, arma::vec velocity)
{
    q = charge;
    m = mass;
    r = position;
    v = velocity;
}