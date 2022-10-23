#include <iostream>
#include <armadillo>
#include <vector>
#include "Particle.hpp"
#include "PenningTrap.hpp"

using namespace std;


int main()
{
    double B0 = 9.65e1; // u us^(-1) e^(-1)
    double V0 = 2.41e6; // u (um)^2 (us)^(-2) e^(-1)
    double d = 500.; // um
    bool particleInteractions = true;
    PenningTrap pt(B0, V0, d, particleInteractions);
    int q = 1; // e
    double m = 20.; // u
    double omega_z = sqrt( 2*q*V0/(m*d*d) ); // (us)^(-1) = MHz
    //std::cout << omega_z << std::endl;
    arma::vec r1("20 0 20");
    arma::vec v1("0 25 0");
    Particle p1(q, m, r1, v1);
    pt.add_particle(p1);
    
    arma::vec r2("25 25 0");
    arma::vec v2("0 40 5");
    Particle p2(q, m, r2, v2);
    
    pt.add_particle(p2);
    
    double T = 50.; //us
    int nSteps = 1000;
    double dt = T/nSteps;
    std::vector<arma::vec> r_list;
    r_list.push_back(p1.r);
    std::vector<arma::vec> v_list;
    v_list.push_back(p1.v);
    
    for (Particle& pi : pt.particles)
    {
        cout << pi.r(0) << " " << pi.r(1) << " " << pi.r(2) << " ";
    }
    cout << endl;
    
    for (int i = 0; i < nSteps; i++)
    {
        pt.evolve_RK4(dt);
        for (Particle& pi : pt.particles)
        {
            cout << pi.r(0) << " " << pi.r(1) << " " << pi.r(2) << " ";
        }
        cout << endl;
    }
    
    return 0;
}