#include "propagator.h"




propagator::propagator()
{
    PI= 3.14159265358979323846264;

    N_y=100;
    N_x=100;
    N_xh=N_x/2;
    N_yh=N_y/2;
    N=N_x*N_y;

    n_u=0.5;
    n_d=0.5;
    n = n_u+n_d;

    U=2.0;
    mu=U/2;
    m=0.25;

    T=0.0;
    beta=1.0/0.0;

    Ep.set_size(N_x,N_y);
    Em.set_size(N_x,N_y);
    epsilon.set_size(N_x,N_y);
    t.set_size(3);
    t(0)=1.0;
    t(1)=.2326;
    t(2)=.1163;


    set_epsilon();
    set_E();
    calc_m();
}



propagator::propagator(Config * parameters)
{

    PI = parameters->lookup("PI");
    ndim = parameters->lookup("ndim");
    N_x=parameters->lookup("N_x");
    N_y = parameters->lookup("N_y");
    N_xh = N_x/2;
    N_yh = N_y/2;
    N=N_x*N_y;

    U = parameters->lookup("U"); // default 2.0

    m = parameters->lookup("m_initial"); // default .45
    n_u= parameters->lookup("n_u"); // default .5
    n_d= parameters->lookup("n_d"); // default .5
    n = n_u + n_d;
    mu = parameters->lookup("mu"); // default .5*U

    T = parameters->lookup("T");
    beta = 1.0/T;


    Ep.set_size(N_x,N_y);
    Em.set_size(N_x,N_y);
    epsilon.set_size(N_x,N_y);


    t.set_size(3);
    t(0)=parameters->lookup("t")[0];
    t(1)=parameters->lookup("t")[1];
    t(2)=parameters->lookup("t")[2];

    set_epsilon();
    set_E();
    calc_m();

    // write the energie dispersion as well as E^{\pm}
    ofstream Epm;
    Epm.open("energies.txt");
    Epm.precision(17);
    for (int i=0;i<N_x;i++)
    {
        for (int j =0;j<N_y;j++)
        {
            Epm<<i<<"\t"<<j<<"\t"<<2*PI*(double(i)/N_x-.5)<<"\t"<<2*PI*(double(j)/N_y-.5)<<"\t"<<epsilon(i,j)<<"\t"<<Ep(i,j)<<"\t"<<Em(i,j)<<endl;
        }
    }
    Epm.close();
}



complex<double> propagator::G(Col<int> p, complex<double> w)
{
    return (w-epsilon((p(0)+N_xh)%N_x,(p(1)+N_yh)%N_y)+mu-U*.5*n) // no difference between up/down, assuming n_u=n_d=n/2
            /((w- Ep(p(0)%N_x, p(1)%N_y) )*(w-Em(p(0)%N_x,p(1)%N_y)));
}

double propagator::ResG(Col<int> p, int pm)
{
    double res = 0.0;
    switch(pm)
    {
    case 1:
        res = (Ep(p(0)%N_x, p(1)%N_y)-epsilon((p(0)+N_xh)%N_x, (p(1)+N_yh)%N_y)+mu-U*.5*n)
                /(Ep(p(0)%N_x, p(1)%N_y)-Em(p(0)%N_x, p(1)%N_y));
        break;

    case -1:
        res = (Em(p(0)%N_x, p(1)%N_y)-epsilon((p(0)+N_xh)%N_x, (p(1)+N_yh)%N_y)+mu-U*.5*n)
                /(Em(p(0)%N_x, p(1)%N_y)-Ep(p(0)%N_x, p(1)%N_y));
        break;

    default:
        cout<<"pm must be +1 or -1"<<endl;

    }
    return res;
}

double propagator::dwResG(Col<int> p, int pm)
{
    double res = pm*(1-ResG(p,pm))/(Ep(p(0)%N_x,p(1)%N_y)-Em(p(0)%N_x,p(1)%N_y));
    return res;
}


complex<double> propagator::F(Col<int> p, int sigma, complex<double> w)
{
    return -sigma*U*m/((w-Ep(p(0)%N_x, p(1)%N_y))*(w-Em(p(0)%N_x, p(1)%N_y)));
}


double propagator::ResF(Col<int> p, int sigma, int pm)
{
    return -sigma*pm*U*m/(Ep(p(0)%N_x, p(1)%N_y)-Em(p(0)%N_x, p(1)%N_y));
}
double propagator::dwResF(Col<int> p, int sigma, int pm)
{
    return sigma*U*m/(pow(Ep(p(0)%N_x,p(1)%N_y)-Em(p(0)%N_x,p(1)%N_y),2));
}

double propagator::get_Ep(Col<int> p)
{
    return Ep(p(0)%N_x, p(1)%N_y);
}

double propagator::get_Em(Col<int> p)
{
    return Em(p(0)%N_x, p(1)%N_y);
}

double propagator::get_epsilon(Col<int> p)
{
    return epsilon(p(0)%N_x, p(1)%N_y);
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
                        if (Em(i,j)<0.0 && Ep(i,j)>0.0)
                        {
                            fm+=1.0/(Ep(i,j)-Em(i,j));
                            dfdm-=2.0/(pow(Ep(i,j)-Em(i,j),3));
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
                        fm+= 1.0/(Ep(i,j)-Em(i,j))*
                                ( -1.0/(1+exp(beta*Ep(i,j)))
                                  +1.0/(1+exp(beta*Em(i,j)))
                                 );

                        dfdm+= 2.0/pow(Ep(i,j)-Em(i,j),3)
                            *(
                                    1.0/( 1.0+exp(beta*Ep(i,j)) ) -1.0/( 1+exp(beta*Em(i,j)) )
                                )
                            +beta/(( Ep(i,j)-Em(i,j))*( Ep(i,j)-Em(i,j)))
                            *(
                                   0.25/pow(cosh(beta*Ep(i,j)/2),2) + 0.25/pow(cosh(beta*Em(i,j)/2),2)
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
            m-=(fm-1)/dfdm*100./(100.+(double)abort);
            if(m>0.5 || m<0.0)
            {
                m=0.01;
            }
            set_E();

            // check m for a start (remove later).
            //cout<<m<<endl;

            // break loop after a maximum limit of iterations
            abort++;
            if(abort>1e3)
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


double propagator::calc_n()
{
    // calculate n_{\sigma} (relative number of spin \sigma particles (0.5 per spin at half filling)
    // the number is assumed to be 1/2 (using n_u). E^{\pm} depend already on n_{\sigma}
    // so this only checks consistency and accuracy, no optimization.
    // Furthermore it is assumed, that E^+>0 and E^-<0
    double n=.0;

    if (T==0.0)
    {
        for (int i=0;i<N_x;i++)
        {
            for (int j=0;j<N_y;j++)
            {
                n+= -1/(Ep(i,j)-Em(i,j))*( Em(i,j) - epsilon(i,j)+mu-U*n_u );
            }
        }
    }
    else
    {
        for (int i=0;i<N_x;i++)
        {
            for (int j=0;j<N_y;j++)
            {
                n+= -1/(Ep(i,j)-Em(i,j))*( Em(i,j) - epsilon(i,j)+mu-U*n_u )/(1+exp(beta*Em(i,j)))
                        +1/(Ep(i,j)-Em(i,j))*( Ep(i,j) - epsilon(i,j)+mu-U*n_u )/(1+exp(beta*Ep(i,j)));
            }
        }
    }
    // check the difference to 1/2,
    // The result should be a last column basically of zeros (up to double precision)
    n=n/N-.5;

    return n;
}


// energy dispersion is defined here
void propagator::set_epsilon()
{
    // calculate energies \varepsilon_{\vec{k}} of the unpertubed system
    for (int i=0;i<N_x;i++)
    {
        for(int j=0;j<N_y;j++)
        {
            // t-t'-t"
            epsilon(i,j) =  -2*t(0)*( cos(2*PI*(double(i)/N_x-.5))+cos(2*PI*(double(j)/N_y-.5)) )
                            -2*t(1)*( cos(4*PI*(double(i)/N_x-.5))+cos(4*PI*(double(j)/N_y-.5)) )
                            -4*t(2)*  cos(2*PI*(double(i)/N_x-.5))*cos(2*PI*(double(j)/N_y-.5)) ;
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
                E_Sh = .5*(epsilon(i,j)+epsilon((i+N_xh)%N_x,(j+N_yh)%N_y)) -mu +U*n_u;
                E_Sd = sqrt( .25*pow((epsilon(i,j)-epsilon((i+N_xh)%N_x,(j+N_yh)%N_y)),2) + U*U*m*m );

                Ep(i,j) = E_Sh+E_Sd;
                Em(i,j) = E_Sh-E_Sd;
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

void propagator::set_U(double newU, double new_mu)
{
    U=newU;
    mu=new_mu;
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


int propagator::get_N(int i)
{
    // quick and dirty
    if (i==0)
        return N_x;

    else if (i==1)
        return N_y;
    else
    {
        cout<<"wrong dimension in calling get_N(i) (only 0,1 are valid)"<<endl;
        return 0;
    }
}

int propagator::get_N()
{
    return N;
}

double propagator::get_n()
{
    return n;
}
