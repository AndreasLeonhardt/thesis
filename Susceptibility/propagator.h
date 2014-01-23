#ifndef PROPAGATOR_H
#define PROPAGATOR_H

#include libconfig.h++

using namespace std;
using namespace libconfig;

class propagator
{
protected:

double PI;

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
int N_x,N_y,N,N_xh,N_yh;

// filling, n \in [0,2]; n_u,n_d \in [0,1]
double n_u,n_d,n;

vector< vector <double> > epsilon;
vector< vector <double> > Ep;
vector< vector <double> > Em;


public:
    propagator();
    propagator(Config *parameters);

    void set_U(double newU);
    double get_U(); // updates Ep,Em and m.

    void set_T(double newT); // updates m
    double get_T();

    double get_beta();

    double get_m();
    double calc_m();

    double get_Ep(int x, int y);
    double get_Em(int x, int y);
    void set_E();
    double get_epsilon(int x, int y);
    void set_epsilon();

    double G(int x, int y, double w);
    double F(int x, int y, double w);

    double ResG(int x, int y, int pm);
    double ResF(int x, int y, int pm);



};

#endif // PROPAGATOR_H
