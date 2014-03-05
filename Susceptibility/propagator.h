#ifndef PROPAGATOR_H
#define PROPAGATOR_H

#include <libconfig.h++>
#include <cmath>
#include <complex>
#include <iostream>
#include "armadillo"

using namespace std;
using namespace libconfig;
using namespace arma;

class propagator
{
protected:

double PI;
int ndim;

// initial value for up-/down- part of the staggered magnetization m_{s,\sigma}, here m,
// without distinguishing the up and down case, e.g. m_u = m_d = m = \frac12 m_{\mathrm{s}}
double m;

// here U corresponds to \frac{U}{t} with t normalized to 1
double U;
// chemical potential mu (half filling at mu=U/2?) n should be calculated from mu or the other way round.
double mu;

// temperature T and inverse temperature beta=1/T
double T;
double beta;

// sites in one dimension (2D square grid)
// number of sites must be even
int N,N_x,N_xh,N_y,N_yh;


// filling, n \in [0,2]; n_u,n_d \in [0,1]
double n_u,n_d,n;

mat epsilon; // defined over the whole Brillouine zone, (0,N_x) x (0,N_y) corresponding (-pi,pi)x(-pi,pi)
mat Ep;
mat Em;


public:
    propagator();
    propagator(Config *parameters);

    void set_U(double newU);
    void set_U(double newU,double new_mu);
    double get_U(); // updates Ep,Em and m.

    void set_T(double newT); // updates m
    double get_T();

    double get_beta();

    double get_m();
    void calc_m();
    double calc_n();


    int get_N(int i);
    int get_N();

    double get_n();

    double get_Ep(Col<int> p);
    double get_Em(Col<int> p);
    void set_E();
    double get_epsilon(Col<int>p);
    void set_epsilon();

    complex<double> G(Col<int> p, complex<double> w);
    complex<double> F(Col<int> p, int sigma, complex<double> w);

    double ResG(Col<int> p, int pm);
    double ResF(Col<int> p, int sigma, int pm);

    double dwResG(Col<int> p, int pm);
    double dwResF(Col<int> p, int sigma, int pm);



};

#endif // PROPAGATOR_H
