// Susceptibility
// description: Calculates first the staggered magnetization according to the given parameters U,t,\beta.
//              Proceeds in calculating \Xi^{+-}(q,\omega) (later \Xi^{zz} as well)
//              and writes the result in a file.


#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <complex>
#include <stdlib.h>
#include <cstring>
#include <libconfig.h++>
#include <vector>
#include "armadillo"


#include "propagator.h"


using namespace std;
using namespace libconfig;
using namespace arma;

// compare to doubles in a meaningful way.
bool is_same(double a, double b);

int main()
{


    // create instance of config and open config file to get parameters
    Config conf_parameters;
    Config * parameters = &conf_parameters;
    conf_parameters.readFile("/home/andi/Studium/masteroppgave/parameters.cfg");

    // create instance of class propagator, setting energies, magnetization and so on
    propagator Prop(parameters);
    propagator * prop = &Prop;
    double PI = parameters->lookup("PI");
    int ndim = parameters->lookup("ndim");
    double eta = parameters->lookup("eta");

    // open result file and write overall parameters
    ofstream results;
    results.open(parameters->lookup("chi"));
    results<<"U= "<<prop->get_U()<<"\tT= "<<prop->get_T()<<"\tn= "<<prop->get_n()<<"\tm= "<<prop->get_m()
           <<"\tN= "<<prop->get_N()<<"\teta="<<eta<<endl;
    results.precision(17);

    // write column title "t q_x q_y w   Re(x)   Im(x)   Re(xb)  Im(xb)  Re(y)   Im(y)  Re(z1)  Im(z1)  Re(z2)  Im(z2)" to the file
    results<<"t\tq_x\tq_y\tw\tRe(x)\tIm(x)\tRe(xb)\tIm(xb)\tRe(y)\tIm(y)\tRe(z1)\tIm(z1)\tRe(z2)\tIm(z2)"<<endl;

    // set q ( from (0,0) to (N_x,N_y) corresponding to to be in (-\pi,\pi)\times (-\pi,\pi)
    // and w, which should be choose in a clever way to give a meaningful result
    // following a path through the Brillouin zone
    // only use existing values for (q_x,q_y) (according to discrete momenta)
    Col<int> q(ndim);
    Col<int> Q(ndim);
    Q(0)=prop->get_N(0)/2; // this is not ready for flexibel dimensions
    Q(1)=prop->get_N(1)/2;

    // range of w
    double w_min, w_max,w_step;
    w_min=parameters->lookup("w_min");
    w_max=parameters->lookup("w_max");
    w_step=parameters->lookup("w_step");

    double beta = prop->get_beta();
    // set q to start values

    // old path starting value
    q[0]=Q[0]/2;
    q[1]=3*Q[1]/2;

    // new path starting value
    q[0]=3*Q[0]/2;
    q[1]=3*Q[1]/2;

    for (int t=0;t<=7*Q[0]/2;t++)
    {
        for(double w=w_min;w<=w_max+w_step;w+=w_step) // in units of t
        {
            complex<double> cw(w,eta); // complex w, that is cw = w + i\eta
            complex<double> x, xb, y, z1, z2; // complex initializes to zero.
            Col<int> p(ndim);

            // calculate x(q,\omega), y, z_1,z_2 (\bar{x}(q,\omega) = x(q+Q,\omega)
            if (prop->get_T()!=0.0)
            {
                for (int px=0;px<2*Q[0];px++)
                {
                    p[0]=px;
                    for (int py=0;py<2*Q[1];py++)
                    {
                        p[1]=py;

                            // variables of often used values to speed up the following calculations
                            double Epp = prop->get_Ep(p);
                            double Emp = prop->get_Em(p);
                            double Eppq = prop->get_Ep(p+q);
                            double Empq = prop->get_Em(p+q);

                         /*   // check order of poles (single poles, double poles)
                            // ASSUMING THAT 0 <= w < U*m, therefore we have ALWAYS Epp > Empq-w; Eppq-w > Emp and Emp > Epp.
                            // leaving us to check that Epp != Eppq-w and Emp != Empq-w.

                            if (is_same(Epp,Eppq-w))
                            {
                                cout <<" double pole (+) at"<< q<<endl;
                                // double pole at Epp
                                x[v] += (prop->dwResG(p,+1)*prop->ResG(p+q+v*Q,+1)
                                      +prop->ResG(p,+1)*prop->dwResG(p+q+v*Q,+1))/(1+exp(beta*Epp))
                                        + prop->ResG(p,+1)*prop->ResG(p+q+v*Q,+1)*(-beta)/(2*cosh(beta*Epp));

                                y[v] += (prop->dwResF(p,-1,+1)*prop->ResF(p+q+v*Q,+1,+1)
                                      +prop->ResF(p,-1,+1)*prop->dwResF(p+q+v*Q,+1,+1))/(1+exp(beta*Epp))
                                        + prop->ResF(p,-1,+1)*prop->ResF(p+q+v*Q,+1,+1)*(-beta)/(2*cosh(beta*Epp));

                                z1[v] += (prop->dwResG(p,+1)*prop->ResF(p+q+v*Q,+1,+1)
                                       +prop->ResG(p,+1)*prop->dwResF(p+q+v*Q,+1,+1))/(1+exp(beta*Epp))
                                        + prop->ResG(p,+1)*prop->ResF(p+q+v*Q,+1,+1)*(-beta)/(2*cosh(beta*Epp));

                                z2[v] += (prop->dwResF(p,-1,+1)*prop->ResG(p+q+v*Q,+1)
                                       +prop->ResF(p,-1,+1)*prop->dwResG(p+q+v*Q,+1))/(1+exp(beta*Epp))
                                        + prop->ResF(p,-1,+1)*prop->ResG(p+q+v*Q,+1)*(-beta)/(2*cosh(beta*Epp));
                            }
                            else
                            {*/
                                // single poles at Epp, Eppq-w
                                x += prop->ResG(p  ,+1)*prop->G(p+q,Epp +cw)/(1+exp(beta*Epp ));
                                x += prop->ResG(p+q,+1)*prop->G(p  ,Eppq-cw)/(1+exp(beta*(Eppq)));

                                xb+= prop->ResG(p  ,+1)*prop->G(p+q+Q,Epp +cw)/(1+exp(beta*Epp ));
                                xb+= prop->ResG(p+q+Q,+1)*prop->G(p  ,Eppq-cw)/(1+exp(beta*(Eppq)));

                                y += prop->ResF(p  ,-1,+1)*prop->F(p+q,+1,Epp +cw)/(1+exp(beta*Epp ));
                                y += prop->ResF(p+q,+1,+1)*prop->F(p  ,-1,Eppq-cw)/(1+exp(beta*(Eppq)));

                                z1 += prop->ResG(p  ,+1)*prop->F(p+q,+1,Epp +cw)/(1+exp(beta*Epp ));
                                z1 += prop->ResF(p+q,+1,+1)*prop->G(p   ,Eppq-cw)/(1+exp(beta*(Eppq)));

                                z2 += prop->ResF(p  ,-1,+1)*prop->G(p+q,Epp +cw)/(1+exp(beta*Epp ));
                                z2 += prop->ResG(p+q,+1)*prop->F(p   ,-1,Eppq-cw)/(1+exp(beta*(Eppq)));
                           /* }

                            if (is_same(Emp,Empq-w))
                            {
                                cout <<" double pole (-) at"<< q<<endl;
                                // pole of second order at w=Emp
                                x[v] += (prop->dwResG(p,-1)*prop->ResG(p+q+v*Q,-1)
                                      +prop->ResG(p,-1)*prop->dwResG(p+q+v*Q,-1))/(1+exp(beta*Emp))
                                        + prop->ResG(p,-1)*prop->ResG(p+q+v*Q,-1)*(-beta)/(2*cosh(beta*Emp));

                                y[v] += (prop->dwResF(p,-1,-1)*prop->ResF(p+q+v*Q,+1,-1)
                                      +prop->ResF(p,-1,-1)*prop->dwResF(p+q+v*Q,+1,-1))/(1+exp(beta*Emp))
                                        + prop->ResF(p,-1,-1)*prop->ResF(p+q+v*Q,+1,-1)*(-beta)/(2*cosh(beta*Emp));

                                z1[v] += (prop->dwResG(p,-1)*prop->ResF(p+q+v*Q,+1,-1)
                                       +prop->ResG(p,-1)*prop->dwResF(p+q+v*Q,+1,-1))/(1+exp(beta*Emp))
                                        + prop->ResG(p,-1)*prop->ResF(p+q+v*Q,+1,-1)*(-beta)/(2*cosh(beta*Emp));

                                z2[v] += (prop->dwResF(p,-1,-1)*prop->ResG(p+q+v*Q,-1)
                                       +prop->ResF(p,-1,-1)*prop->dwResG(p+q+v*Q,-1))/(1+exp(beta*Emp))
                                        + prop->ResF(p,-1,-1)*prop->ResG(p+q+v*Q,-1)*(-beta)/(2*cosh(beta*Emp));
                            }
                            else
                            {*/
                                // single poles at Emp, Empq-w
                                x += prop->ResG(p  ,-1)*prop->G(p+q,Emp +cw)/(1+exp(beta*Emp ));
                                x += prop->ResG(p+q,-1)*prop->G(p  ,Empq-cw)/(1+exp(beta*(Empq)));

                                xb+= prop->ResG(p  ,-1)*prop->G(p+q+Q,Emp +cw)/(1+exp(beta*Emp ));
                                xb+= prop->ResG(p+q+Q,-1)*prop->G(p  ,Empq-cw)/(1+exp(beta*(Empq)));

                                y += prop->ResF(p  ,-1,-1)*prop->F(p+q,+1,Emp +cw)/(1+exp(beta*Emp ));
                                y += prop->ResF(p+q,+1,-1)*prop->F(p  ,-1,Empq-cw)/(1+exp(beta*(Empq)));

                                z1+= prop->ResG(p  ,-1)*prop->F(p+q,+1,Emp +cw)/(1+exp(beta*Emp ));
                                z1+= prop->ResF(p+q,+1,-1)*prop->G(p   ,Empq-cw)/(1+exp(beta*(Empq)));

                                z2+= prop->ResF(p  ,-1,-1)*prop->G(p+q,Emp +cw)/(1+exp(beta*Emp ));
                                z2+= prop->ResG(p+q,-1)*prop->F(p   ,-1,Empq-cw)/(1+exp(beta*(Empq)));
                            //}

                    }
                }
            }
            else // T=0 IS NOT UPDATED TO THE PROPER TREATMENT OF SECOND ORDER POLES. SEE ABOVE FOR PROPER TREATMENT.
                {
                    cout << "T=0.0 not fixed for proper pole treatment"<<endl;

                    // same but for T=0, e.g. sharp fermionic distribution (wohoo, almost only half the terms)
                    for (int px=0;px<2*Q[0];px++)
                    {
                        p[0]=px;
                        for (int py=0;py<2*Q[1];py++)
                        {
                            p[1]=py;

                            for (int v=0;v<2;v++)
                            {
                                double Epp = prop->get_Ep(p);
                                double Emp = prop->get_Em(p);
                                double Eppq = prop->get_Ep(p+q);
                                double Empq = prop->get_Em(p+q);

                                if (Epp<0.0)
                                {
                                    x +=  prop->ResG(p,+1)*prop->G(p+q,Epp +cw);
                                    y +=  prop->ResF(p,-1,+1)*prop->F(p+q,+1,Epp +cw);
                                    z1+=  prop->ResG(p,+1)*prop->F(p+q,+1,Epp +cw);
                                    z2+=  prop->ResF(p,-1,+1)*prop->G(p+q,Epp +cw);
                                }

                                if(Emp<0.0)
                                {
                                    x += prop->ResG(p,-1)*prop->G(p+q,Emp +cw);
                                    y += prop->ResF(p,-1,-1)*prop->F(p+q,+1,Emp +cw);
                                    z1+= prop->ResG(p,-1)*prop->F(p+q,+1,Emp +cw);
                                    z2+= prop->ResF(p,-1,-1)*prop->G(p+q,Emp +cw);
                                }

                                if(Eppq<0.0)
                                {
                                    x +=prop->ResG(p+q,+1)*prop->G(p   ,Eppq-cw);
                                    y +=prop->ResF(p+q,+1,+1)*prop->F(p   ,-1,Eppq-cw);
                                    z1+=prop->ResF(p+q,+1,+1)*prop->G(p   ,Eppq-cw);
                                    z2+=prop->ResG(p+q,+1)*prop->F(p   ,-1,Eppq-cw);
                                }

                                if(Empq<0.0)
                                {
                                    x +=prop->ResG(p+q,-1)*prop->G(p   ,Empq-cw);
                                    y +=prop->ResF(p+q,+1,-1)*prop->F(p   ,-1,Empq-cw);
                                    z1+=prop->ResF(p+q,+1,-1)*prop->G(p   ,Empq-cw);
                                    z2+=prop->ResG(p+q,-1)*prop->F(p   ,-1,Empq-cw);
                                }
                            }
                        }
                    }
                }



            results<<t<<"\t"<<q[0]*2*PI/prop->get_N(0)-PI<<"\t"<<q[1]*2*PI/prop->get_N(1)-PI<<"\t"<<w<<"\t"
                   <<real( x)<<"\t"<<imag( x)<<"\t"
                   <<real(xb)<<"\t"<<imag(xb)<<"\t"
                   <<real( y)<<"\t"<<imag( y)<<"\t"
                   <<real(z1)<<"\t"<<imag(z1)<<"\t"
                   <<real(z2)<<"\t"<<imag(z2)<<"\t"<<endl;
        }



        // move to next position in q, depending on the path.

//        old path (as in La_2CuO4 paper)
//        if (t<Q[0]/2)
//        {
//            q[0] -=1;
//            q[1] +=1;
//        }
//        else if (t<3*Q[0]/2)
//        {
//            q[1]-=1;
//        }
//        else if (t<4*Q[0]/2)
//        {
//            q[0]+=1;
//            q[1]+=1;
//        }
//        else if (t<5*Q[0]/2)
//        {
//            q[0]+=1;
//            q[1]-=1;
//        }
//        else if(t<=7*Q[0]/2)
//        {
//            q[0]-=1;
//        }

                if (t<Q[0]/2)
                {
                    q[0] +=1;
                    q[1] -=1;
                }
                else if (t<3*Q[0]/2)
                {
                    q[1]+=1;
                }
                else if (t<5*Q[0]/2)
                {
                    q[0]-=1;
                    q[1]-=1;
                }
                else if(t<=7*Q[0]/2)
                {
                    q[0]+=1;
                }

        // print progress (only dependent on t, e.g. the q in the walk through the BZ
        cout<<"\r"<<100*t/(7*Q[0]/2)<<"%"<<flush;


    }

    cout<<endl;
    results.close();
    return 0;
}



// This function compares two doubles in a meaningful way by
// comparing the differnce of a and b to the smallest notable difference
// between max(a,b) and the next bigger number.
// note: could be enhanced using the integer representation.
//          or simplified by fixing epsilon, since we have a quite limited data range.
bool is_same(double a, double b)
{
    // get difference between 1.0 and next bigger value
    double eps = numeric_limits<double>::epsilon()*2.0;
    // absolut value of difference between input numbers
    double diff = abs(a-b);
    // choose number with bigger absolut value
    a=abs(a);
    b=abs(b);
    double max = (a>b) ? a : b;

    // compare difference to scaled epsilon
    if( diff<max*eps)
        return true;

    // if it fails, they are unequal.
    return false;
}
