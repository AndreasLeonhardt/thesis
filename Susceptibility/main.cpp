// Susceptibility
// description: Calculates first the staggered magnetization according to the given parameters U,t,\beta.
//              Proceeds in calculating \Xi^{+-}(q,\omega) (later \Xi^{zz} as well)
//              and writes the result in a file.


#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
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

    // open result file and write overall parameters
    ofstream results;
    results.open(parameters->lookup("chi"));
    results<<"U= "<<prop->get_U()<<"\tT= "<<prop->get_T()<<"\tn= "<<prop->get_n()<<"\tm= "<<prop->get_m()<<"\tN= "<<prop->get_N()<<endl;
    results.precision(17);

    // write "t   q_x   q_y   w   x   xb   y   z1   z2" to the file
    results<<"t\tq_x\tq_y\tw\tx\txb\ty\tz1\tz2"<<endl;

    // set q ( from (0,0) to (N_x,N_y) corresponding to to be in (-\pi,\pi)\times (-\pi,\pi)
    // and w, which should be choose in a clever way to give a meaningful result
    // following a path through the Brillouin zone
    // only use existing values for (q_x,q_y) (according to discrete momenta)
    Col<int> q(ndim);
    Col<int> Q(ndim);
    Q(0)=prop->get_N(0)/2; // this is not ready for flexibel dimensions
    Q(1)=prop->get_N(1)/2;

    double beta = prop->get_beta();

    // set q to start values
    q[0]=Q[0]/2;
    q[1]=3*Q[1]/2;

    for (int t =0;t<7*Q[0]/2;t++)
    {
        for(double w=0.0;w<2.0;w+=0.05) // in units of t
        {
            double x=0.0,xb=0.0, y=0.0, z1=0.0, z2=0.0;
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

//                        double ResGpp = prop->ResG(p,+1);
//                        double ResGppq= prop->ResG(p+q,+1);
//                        double ResGmp = prop->ResG(p,-1);
//                        double ResGmpq= prop->ResG(p+q,-1);

//                        double ResFpp = prop->ResF(p,+1,+1);
//                        double ResFppq= prop->ResF(p+q,+1,+1);
//                        double ResFmp = prop->ResF(p,+1,-1);
//                        double ResFmpq= prop->ResF(p+q,+1,-1);

//                        double BoltzP = 1/(1+exp(beta*Epp));
//                        double BoltzM = 1/(1+exp(beta*Epp));

                        // check order of poles (single poles, double poles)
                        // ASSUMING THAT 0 <= w < U*m, therefore we have ALWAYS Epp > Empq-w; Eppq-w > Emp and Emp > Epp.
                        // leaving us to check that Epp != Eppq-w and Emp != Empq-w.

                        if (is_same(Epp,Eppq-w))
                        {
                            // double pole at Epp
                            x += (prop->dwResG(p,+1)*prop->ResG(p+q,+1)
                                  +prop->ResG(p,+1)*prop->dwResG(p+q,+1))/(1+exp(beta*Epp))
                                    + prop->ResG(p,+1)*prop->ResG(p+q,+1)*(-beta)/(2*cosh(beta*Epp));

                            xb+= (prop->dwResG(p,+1)*prop->ResG(p+q+Q,+1)
                                  +prop->ResG(p,+1)*prop->dwResG(p+q+Q,+1))/(1+exp(beta*Epp))
                                    + prop->ResG(p,+1)*prop->ResG(p+q+Q,+1)*(-beta)/(2*cosh(beta*Epp));

                            y += (prop->dwResF(p,-1,+1)*prop->ResF(p+q,+1,+1)
                                  +prop->ResF(p,-1,+1)*prop->dwResF(p+q,+1,+1))/(1+exp(beta*Epp))
                                    + prop->ResF(p,-1,+1)*prop->ResF(p+q,+1,+1)*(-beta)/(2*cosh(beta*Epp));

                            z1 += (prop->dwResG(p,+1)*prop->ResF(p+q,+1,+1)
                                   +prop->ResG(p,+1)*prop->dwResF(p+q,+1,+1))/(1+exp(beta*Epp))
                                    + prop->ResG(p,+1)*prop->ResF(p+q,+1,+1)*(-beta)/(2*cosh(beta*Epp));

                            z2 += (prop->dwResF(p,-1,+1)*prop->ResG(p+q,+1)
                                   +prop->ResF(p,-1,+1)*prop->dwResG(p+q,+1))/(1+exp(beta*Epp))
                                    + prop->ResF(p,-1,+1)*prop->ResG(p+q,+1)*(-beta)/(2*cosh(beta*Epp));
                        }
                        else
                        {
                            // single poles at Epp, Eppq-w
                            x += prop->ResG(p  ,+1)*prop->G(p+q,Epp +w)/(1+exp(beta*Epp ));
                            x += prop->ResG(p+q,+1)*prop->G(p  ,Eppq-w)/(1+exp(beta*(Eppq-w)));

                            xb+= prop->ResG(p    ,+1)*prop->G(p+q+Q,Epp +w)/(1+exp(beta*Epp ));
                            xb+= prop->ResG(p+q+Q,+1)*prop->G(p    ,Eppq-w)/(1+exp(beta*(Eppq-w)));

                            y += prop->ResF(p  ,-1,+1)*prop->F(p+q,+1,Epp +w)/(1+exp(beta*Epp ));
                            y += prop->ResF(p+q,+1,+1)*prop->F(p  ,-1,Eppq-w)/(1+exp(beta*(Eppq-w)));

                            z1+= prop->ResG(p  ,+1)*prop->F(p+q,+1,Epp +w)/(1+exp(beta*Epp ));
                            z1+= prop->ResF(p+q,+1,+1)*prop->G(p   ,Eppq-w)/(1+exp(beta*(Eppq-w)));

                            z2+= prop->ResF(p  ,-1,+1)*prop->G(p+q,Epp +w)/(1+exp(beta*Epp ));
                            z2+= prop->ResG(p+q,+1)*prop->F(p   ,-1,Eppq-w)/(1+exp(beta*(Eppq-w)));
                        }
                        if (is_same(Emp,Empq-w))
                        {
                            // pole of second order at w=Emp
                            x += (prop->dwResG(p,-1)*prop->ResG(p+q,-1)
                                  +prop->ResG(p,-1)*prop->dwResG(p+q,-1))/(1+exp(beta*Emp))
                                    + prop->ResG(p,-1)*prop->ResG(p+q,-1)*(-beta)/(2*cosh(beta*Emp));

                            xb+= (prop->dwResG(p,-1)*prop->ResG(p+q+Q,-1)
                                  +prop->ResG(p,-1)*prop->dwResG(p+q+Q,-1))/(1+exp(beta*Emp))
                                    + prop->ResG(p,-1)*prop->ResG(p+q+Q,-1)*(-beta)/(2*cosh(beta*Emp));

                            y += (prop->dwResF(p,-1,-1)*prop->ResF(p+q,+1,-1)
                                  +prop->ResF(p,-1,-1)*prop->dwResF(p+q,+1,-1))/(1+exp(beta*Emp))
                                    + prop->ResF(p,-1,-1)*prop->ResF(p+q,+1,-1)*(-beta)/(2*cosh(beta*Emp));

                            z1 += (prop->dwResG(p,-1)*prop->ResF(p+q,+1,-1)
                                   +prop->ResG(p,-1)*prop->dwResF(p+q,+1,-1))/(1+exp(beta*Emp))
                                    + prop->ResG(p,-1)*prop->ResF(p+q,+1,-1)*(-beta)/(2*cosh(beta*Emp));

                            z2 += (prop->dwResF(p,-1,-1)*prop->ResG(p+q,-1)
                                   +prop->ResF(p,-1,-1)*prop->dwResG(p+q,-1))/(1+exp(beta*Emp))
                                    + prop->ResF(p,-1,-1)*prop->ResG(p+q,-1)*(-beta)/(2*cosh(beta*Emp));
                        }
                        else
                        {
                            // single poles at Emp, Empq-w
                            x += prop->ResG(p  ,-1)*prop->G(p+q,Emp +w)/(1+exp(beta*Emp ));
                            x += prop->ResG(p+q,-1)*prop->G(p  ,Empq-w)/(1+exp(beta*(Empq-w)));

                            xb+= prop->ResG(p    ,-1)*prop->G(p+q+Q,Emp +w)/(1+exp(beta*Emp ));
                            xb+= prop->ResG(p+q+Q,-1)*prop->G(p    ,Empq-w)/(1+exp(beta*(Empq-w)));

                            y += prop->ResF(p  ,-1,-1)*prop->F(p+q,+1,Emp +w)/(1+exp(beta*Emp ));
                            y += prop->ResF(p+q,+1,-1)*prop->F(p  ,-1,Empq-w)/(1+exp(beta*(Empq-w)));

                            z1+= prop->ResG(p  ,-1)*prop->F(p+q,+1,Emp +w)/(1+exp(beta*Emp ));
                            z1+= prop->ResF(p+q,+1,-1)*prop->G(p   ,Empq-w)/(1+exp(beta*(Empq-w)));

                            z2+= prop->ResF(p  ,-1,-1)*prop->G(p+q,Emp +w)/(1+exp(beta*Emp ));
                            z2+= prop->ResG(p+q,-1)*prop->F(p   ,-1,Empq-w)/(1+exp(beta*(Empq-w)));
                        }
                    }
                }
            }
            else // T=0 IS NOT UPDATED TO THE PROPER TREATMENT OF SECOND ORDER POLES. SEE ABOVE FOR PROPER TREATMENT.
            {
                // same but for T=0, e.g. sharp fermionic distribution (wohoo, almost only half the terms)
                for (int px=0;px<Q[0];px++)
                {
                    p[0]=px;
                    for (int py=0;py<Q[1];py++)
                    {
                        p[1]=py;

                        double Epp = prop->get_Ep(p);
                        double Emp = prop->get_Em(p);
                        double Eppq = prop->get_Ep(p+q);
                        double Empq = prop->get_Em(p+q);

                        if (Epp<0.0)
                        {
                            x +=  prop->ResG(p,+1)*prop->G(p+q,Epp +w);
                            xb+=  prop->ResG(p,+1)*prop->G(p+q+Q,Epp +w);
                            y +=  prop->ResF(p,-1,+1)*prop->F(p+q,+1,Epp +w);
                            z1+=  prop->ResG(p,+1)*prop->F(p+q,+1,Epp +w);
                            z2+=  prop->ResF(p,-1,+1)*prop->G(p+q,Epp +w);
                        }

                        if(Emp<0.0)
                        {
                            x += prop->ResG(p,-1)*prop->G(p+q,Emp +w);
                            xb+= prop->ResG(p,-1)*prop->G(p+q+Q,Emp +w);
                            y += prop->ResF(p,-1,-1)*prop->F(p+q,+1,Emp +w);
                            z1+= prop->ResG(p,-1)*prop->F(p+q,+1,Emp +w);
                            z2+= prop->ResF(p,-1,-1)*prop->G(p+q,Emp +w);
                        }

                        if(Eppq-w<0.0)
                        {
                            x +=prop->ResG(p+q,+1)*prop->G(p   ,Eppq-w);
                            xb+=prop->ResG(p+q+Q,+1)*prop->G(p   ,Eppq-w);
                            y +=prop->ResF(p+q,+1,+1)*prop->F(p   ,-1,Eppq-w);
                            z1+=prop->ResF(p+q,+1,+1)*prop->G(p   ,Eppq-w);
                            z2+=prop->ResG(p+q,+1)*prop->F(p   ,-1,Eppq-w);
                        }

                        if(Empq-w<0.0)
                        {
                            x +=prop->ResG(p+q,-1)*prop->G(p   ,Empq-w);
                            xb+=prop->ResG(p+q+Q,-1)*prop->G(p   ,Empq-w);
                            y +=prop->ResF(p+q,+1,-1)*prop->F(p   ,-1,Empq-w);
                            z1+=prop->ResF(p+q,+1,-1)*prop->G(p   ,Empq-w);
                            z2+=prop->ResG(p+q,-1)*prop->F(p   ,-1,Empq-w);
                        }
                    }
                }
            }

            results<<t<<"\t"<<q[0]*2*PI/prop->get_N(0)-PI<<"\t"<<q[1]*2*PI/prop->get_N(1)-PI<<"\t"<<w
                  <<"\t"<<x<<"\t"<<xb<<"\t"<<y<<"\t"<<z1<<"\t"<<z2<<endl;

        }

        // move to next position in q, depending on the path.
        if (t<Q[0]/2)
        {
            q[0] -=1;
            q[1] +=1;
        }
        else if (t<3*Q[0]/2)
        {
            q[1]-=1;
        }
        else if (t<2*Q[0])
        {
            q[0]+=1;
            q[1]+=1;
        }
        else if (t<5*Q[0]/2)
        {
            q[0]+=1;
            q[1]-=1;
        }
        else if(t<7*Q[0]/2)
        {
            q[0]-=1;
        }
    }



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
