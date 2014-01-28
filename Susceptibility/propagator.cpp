#include "propagator.h"




propagator::propagator()
{
    PI= 3.14159265358979323846264;

    N_x=100;
    N_y=100;
    N_xh=N_x/2;
    N_yh=N_yh;
    N=N_x*N_y;

    n_u=0.5;
    n_d=0.5;
    n = n_u+n_d;

    U=2.0;
    mu=U/2;
    m=0.25;

    T=0.0;
    beta=1.0/0.0;

    Ep.resize(N_x);
    Em.resize(N_x);
    epsilon.resize(N_x);
    for (int i=0;i<N_x;i++)
    {
        Ep[i].resize(N_y);
        Em[i].resize(N_y);
        epsilon[i].resize(N_y);
    }

    set_epsilon();
    set_E();
    calc_m();
}



propagator::propagator(Config * parameters)
{

    PI = parameters->lookup("PI");
    N_x = parameters->lookup("N_x");
    N_xh = N_x/2;
    N_y = parameters->lookup("N_y");
    N_yh = N_y/2;
    N = N_x*N_y;
    U = parameters->lookup("U"); // default 2.0

    m = parameters->lookup("m_initial"); // default .45
    n_u= parameters->lookup("n_u"); // default .5
    n_d= parameters->lookup("n_d"); // default .5
    n = n_u + n_d;
    mu = parameters->lookup("mu"); // default .5*U

    T = parameters->lookup("T");
    beta = 1.0/T;


    Ep.resize(N_x);
    Em.resize(N_x);
    epsilon.resize(N_x);
    for (int i=0;i<N_x;i++)
    {
        Ep[i].resize(N_y); // extend to reduced BZ
        Em[i].resize(N_y);
        epsilon[i].resize(N_y); // extend to full BZ
    }

    set_epsilon();
    set_E();
    calc_m();
}



double propagator::G(int x, int y, double w)
{
    return (w-epsilon[(x+N_xh)%N_x][(y+N_yh)%N_y]+mu-U*.5*n) // no difference between up/down, assuming n_u=n_d=n/2
            /((w-Ep[x%N_x][y%N_y])*(w-Em[x%N_x][y%N_y]));
}

double propagator::ResG(int x, int y, int pm)
{
    double res = 0.0;
    switch(pm)
    {
    case 1:
        res = (Ep[x%N_x][y%N_y]-epsilon[(x+N_xh)%N_x][(y+N_yh)%N_y]+mu-U*.5*n)
                /(Ep[x%N_x][y%N_y]-Em[x%N_x][y%N_y]);
        break;

    case -1:
        res = (Em[x%N_x][y%N_y]-epsilon[(x+N_xh)%N_x][(y+N_yh)%N_y]+mu-U*.5*n)
                /(Em[x%N_x][y%N_y]-Ep[x%N_x][y%N_y]);
        break;

    default:
        cout<<"pm must be +1 or -1"<<endl;

    }
    return res;
}


double propagator::F(int x, int y, double w)
{
    return U*m/((w-Ep[x%N_x][y%N_y])*(w-Em[x%N_x][y%N_y]));
}


double propagator::ResF(int x, int y, int pm)
{
    return pm*U*m/(Ep[x%N_x][y%N_y]-Em[x%N_x][y%N_y]);
}



double propagator::get_Ep(int x, int y)
{
    return Ep[x%N_x][y%N_y];
}

double propagator::get_Em(int x, int y)
{
    return Em[x%N_x][y%N_y];
}

double propagator::get_epsilon(int x, int y)
{
    return epsilon[x%N_x][y%N_y];
}


void propagator::calc_m() // using Newton-Raphson
{
    double fm;
    double dfdm;

    // inverse f(m)=1; (Newton Raphson, weighted)
    int abort = 0;
    double mdiff=0.25;

    while (mdiff>1e-7)
    {
        // calculate m_u and m_d
        fm=0.0;
        dfdm=0.0;

        if (T==0.0)
        {
            for (int i=0;i<N_x;i++)
            {
                for (int j=0;j<N_y;j++)
                {
                    if (Em[i][j]<0.0 && Ep[i][j]>0.0)
                    {
                        fm+=1.0/(Ep[i][j]-Em[i][j]);
                        dfdm-=2.0/(pow(Ep[i][j]-Em[i][j],3));
                    }
                    else
                    {
                        // skip corresponding momentum
                    }

                }
            }
        }
        else // T!=0
        {
            for (int i=0;i<N_x;i++)
            {
                for (int j=0;j<N_y;j++)
                {
                    fm+= 1.0/(Ep[i][j]-Em[i][j])*
                            ( -1.0/(1+exp(beta*Ep[i][j]))
                              +1.0/(1+exp(beta*Em[i][j]))
                             );

                    dfdm+= 2.0/pow(Ep[i][j]-Em[i][j],3)
                          *(
                                1.0/( 1.0+exp(beta*Ep[i][j]) ) -1.0/( 1+exp(beta*Em[i][j]) )
                            )
                           +beta/(( Ep[i][j]-Em[i][j])*( Ep[i][j]-Em[i][j]))
                           *(
                               0.25/pow(cosh(beta*Ep[i][j]/2),2) + 0.25/pow(cosh(beta*Em[i][j]/2),2)
                             );
                }
            }
        }

        // calculate new m
        fm   *= U/(double)N;
        dfdm *= 2.0*U*U*U*m/(double)N;

        // next m with Newton-Raphson
        // check difference, with infinite recursive filter, to not only rely on a single lucky jump
        mdiff=0.1*mdiff+abs((fm-1)/dfdm);
        // set new m, weighted with A/(A+abort) to avoid jumping forth and back
        m-=(fm-1)/dfdm*1.0e4/(1.0e4+(double)abort);
        set_E();

        // check m for a start (remove later).
        //cout<<m<<endl;

        // break loop after a maximum limit of iterations
        abort++;
        if(abort>1e4)
        {
            cout<<"optimization of m reached limit"<<endl;
            break;
        }
        else
        {
            // continue
        }


}
}


void propagator::set_epsilon()
{
    // calculate energies \varepsilon_{\vec{k}} of the unpertubed system
    for (int i=0;i<N_x;i++)
    {
        for(int j=0;j<N_y;j++)
        {
            epsilon[i][j]= -2*cos(2*PI* (double(i)/N_x-.5))-2*cos(2*PI* (double(j)/N_y-.5));
        }
    }
}

void propagator::set_E()
{
    double E_Sh, E_Sd;

        for (int i=0;i<N_x;i++)
        {
            for (int j=0;j<N_y;j++)
            {
                E_Sh = .5*(epsilon[i][j]+epsilon[(i+N_xh)%N_x][(j+N_yh)%N_y]) -mu +U*n_u;
                E_Sd = sqrt( .25*pow((epsilon[i][j]-epsilon[(i+N_xh)%N_x][(j+N_yh)%N_y]),2) + U*U*m*m );

                Ep[i][j] = E_Sh+E_Sd;
                Em[i][j] = E_Sh-E_Sd;
            }
        }
}


//-------------------------------------------------------------------------------------------
// get (and sometimes set) values of non-vector class members ------------------------------------------

double propagator::get_m()
{
    return m;
}



void propagator::set_U(double newU)
{
    U=newU;
    set_E();
    calc_m();
}

double propagator::get_U()
{
    return U;
}


void propagator::set_T(double newT)
{
    T=newT;
    beta=1.0/T;
    calc_m();

}

double propagator::get_T()
{
    return T;
}


double propagator::get_beta()
{
    return beta;
}


int propagator::get_N_x()
{
    return N_x;
}


int propagator::get_N_y()
{
    return N_y;
}

int propagator::get_N()
{
    return N;
}

double propagator::get_n()
{
    return n;
}
