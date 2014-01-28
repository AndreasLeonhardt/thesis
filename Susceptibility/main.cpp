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


#include "propagator.h"


using namespace std;
using namespace libconfig;

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

    // open result file and write overall parameters
    ofstream results;
    results.open(parameters->lookup("chi"));
    results<<"U= "<<prop->get_U()<<"\tT= "<<prop->get_T()<<"\tn= "<<prop->get_n()<<"\tm= "<<prop->get_m()<<"\tN= "<<prop->get_N()<<endl;
    results.precision(17);

    // write "t   q_x   q_y   w   x   xb   y   z1   z2" to the file
    results<<"t\tq_x\tq_y\tw\tx\txb\ty\tz1\tz2"<<endl;

    // set q ( from (0,0) to (N_x,N_y) corresponding to to be in (-\pi,\pi)\times (-\pi,\pi)
    // and w, which should be choose in a clever way to give a meaningful result
    // so far it is only a diagonal straight line through the Brillouin zone.
    // only use existing values for (q_x,q_y) (according to discrete momenta)
    int qx,qy;
    int Qx=prop->get_N_x();
    int Qy=prop->get_N_y();
    double Epp,Emp,Eppq,Empq;
    double beta = prop->get_beta();

    // set qx and qy to start values
    qx=Qx/4;
    qy=3*Qy/4;

    for (int t =0;t<7*Qx/4;t++)
    {


        for(double w=0.0;w<5.0;w+=0.05) // what is a reasonable range in natural units?
        {


            double x=0.0,xb=0.0, y=0.0, z1=0.0, z2=0.0;

            // calculate x(q,\omega), y, z_1,z_2 (\bar{x}(q,\omega) = x(q+Q,\omega)
            if (prop->get_T()!=0.0)
            {
                for (int px=0;px<Qx;px++)
                {
                    for (int py=0;py<Qy;py++)
                    {   Epp = prop->get_Ep(px,py);
                        Emp = prop->get_Em(px,py);
                        Eppq = prop->get_Ep(px+qx,py+qy);
                        Empq = prop->get_Em(px+qx,py+qy);


                        x +=     prop->ResG(px   ,py   ,+1)*prop->G(px+qx,py+qy,Epp +w)/(1+exp(beta*Epp ))
                                +prop->ResG(px   ,py   ,-1)*prop->G(px+qx,py+qy,Emp +w)/(1+exp(beta*Emp ))
                                +prop->ResG(px+qx,py+qy,+1)*prop->G(px   ,py   ,Eppq-w)/(1+exp(beta*(Eppq-w)))
                                +prop->ResG(px+qx,py+qy,-1)*prop->G(px   ,py   ,Empq-w)/(1+exp(beta*(Empq-w)));

                        xb+=     prop->ResG(px   ,py   ,+1)*prop->G(px+qx+Qx,py+qy+Qy,Epp +w)/(1+exp(beta*Epp ))
                                +prop->ResG(px   ,py   ,-1)*prop->G(px+qx+Qx,py+qy+Qy,Emp +w)/(1+exp(beta*Emp ))
                                +prop->ResG(px+qx+Qx,py+qy+Qy,+1)*prop->G(px   ,py   ,Eppq-w)/(1+exp(beta*(Eppq-w)))
                                +prop->ResG(px+qx+Qx,py+qy+Qy,-1)*prop->G(px   ,py   ,Empq-w)/(1+exp(beta*(Empq-w)));


                        // …and y
                        y +=     prop->ResF(px   ,py   ,+1)*prop->F(px+qx,py+qy,Epp +w)/(1+exp(beta*Epp ))
                                +prop->ResF(px   ,py   ,-1)*prop->F(px+qx,py+qy,Emp +w)/(1+exp(beta*Emp ))
                                +prop->ResF(px+qx,py+qy,+1)*prop->F(px   ,py   ,Eppq-w)/(1+exp(beta*(Eppq-w)))
                                +prop->ResF(px+qx,py+qy,-1)*prop->F(px   ,py   ,Empq-w)/(1+exp(beta*(Empq-w)));

                        // and z_1 and z_2
                        z1+=     prop->ResG(px   ,py   ,+1)*prop->F(px+qx,py+qy,Epp +w)/(1+exp(beta*Epp ))
                                +prop->ResG(px   ,py   ,-1)*prop->F(px+qx,py+qy,Emp +w)/(1+exp(beta*Emp ))
                                +prop->ResF(px+qx,py+qy,+1)*prop->G(px   ,py   ,Eppq-w)/(1+exp(beta*(Eppq-w)))
                                +prop->ResF(px+qx,py+qy,-1)*prop->G(px   ,py   ,Empq-w)/(1+exp(beta*(Empq-w)));

                        z2+=     prop->ResF(px   ,py   ,+1)*prop->G(px+qx,py+qy,Epp +w)/(1+exp(beta*Epp ))
                                +prop->ResF(px   ,py   ,-1)*prop->G(px+qx,py+qy,Emp +w)/(1+exp(beta*Emp ))
                                +prop->ResG(px+qx,py+qy,+1)*prop->F(px   ,py   ,Eppq-w)/(1+exp(beta*(Eppq-w)))
                                +prop->ResG(px+qx,py+qy,-1)*prop->F(px   ,py   ,Empq-w)/(1+exp(beta*(Empq-w)));
                    }
                }
            }
            else
            {
                // same but for T=0, e.g. sharp fermionic distribution (wohoo, almost only half the terms)
                for (int px=0;px<Qx;px++)
                {
                    for (int py=0;py<Qy;py++)
                    {
                        Epp = prop->get_Ep(px,py);
                        Emp = prop->get_Em(px,py);
                        Eppq = prop->get_Ep(px+qx,py+qy);
                        Empq = prop->get_Em(px+qx,py+qy);

                        if (Epp<0.0)
                        {
                            x +=  prop->ResG(px   ,py   ,+1)*prop->G(px+qx,py+qy,Epp +w);
                            xb+=  prop->ResG(px   ,py   ,+1)*prop->G(px+qx+Qx,py+qy+Qy,Epp +w);
                            y +=  prop->ResF(px   ,py   ,+1)*prop->F(px+qx,py+qy,Epp +w);
                            z1+=  prop->ResG(px   ,py   ,+1)*prop->F(px+qx,py+qy,Epp +w);
                            z2+=  prop->ResF(px   ,py   ,+1)*prop->G(px+qx,py+qy,Epp +w);
                        }

                        if(Emp<0.0)
                        {
                            x += prop->ResG(px,py,-1)*prop->G(px+qx,py+qy,Emp +w);
                            xb+= prop->ResG(px   ,py   ,-1)*prop->G(px+qx+Qx,py+qy+Qy,Emp +w);
                            y += prop->ResF(px   ,py   ,-1)*prop->F(px+qx,py+qy,Emp +w);
                            z1+= prop->ResG(px   ,py   ,-1)*prop->F(px+qx,py+qy,Emp +w);
                            z2+= prop->ResF(px   ,py   ,-1)*prop->G(px+qx,py+qy,Emp +w);
                        }

                        if(Eppq-w<0.0)
                        {
                            x +=prop->ResG(px+qx,py+qy,+1)*prop->G(px   ,py   ,Eppq-w);
                            xb+=prop->ResG(px+qx+Qx,py+qy+Qy,+1)*prop->G(px   ,py   ,Eppq-w);
                            y +=prop->ResF(px+qx,py+qy,+1)*prop->F(px   ,py   ,Eppq-w);
                            z1+=prop->ResF(px+qx,py+qy,+1)*prop->G(px   ,py   ,Eppq-w);
                            z2+=prop->ResG(px+qx,py+qy,+1)*prop->F(px   ,py   ,Eppq-w);
                        }

                        if(Empq-w<0.0)
                        {
                            x +=prop->ResG(px+qx,py+qy,-1)*prop->G(px   ,py   ,Empq-w);
                            xb+=prop->ResG(px+qx+Qx,py+qy+Qy,-1)*prop->G(px   ,py   ,Empq-w);
                            y +=prop->ResF(px+qx,py+qy,-1)*prop->F(px   ,py   ,Empq-w);
                            z1+=prop->ResF(px+qx,py+qy,-1)*prop->G(px   ,py   ,Empq-w);
                            z2+=prop->ResG(px+qx,py+qy,-1)*prop->F(px   ,py   ,Empq-w);
                        }
                    }
                }
            }

            results<<t<<"\t"<<qx*2*PI/prop->get_N_x()-PI<<"\t"<<qy*2*PI/prop->get_N_y()-PI<<"\t"<<w
                  <<"\t"<<x<<"\t"<<xb<<"\t"<<y<<"\t"<<z1<<"\t"<<z2<<endl;

        }

        // move to next position in q, depending on the path.
        if (t<Qx/4)
        {
            qx -=1; // diagonal path through the Brillouin zone
            qy +=1;
        }
        else if (t<3*Qx/4)
        {
            qy-=1;
        }
        else if (t<Qy)
        {
            qx+=1;
            qy+=1;
        }
        else if (t<5*Qy/4)
        {
            qx+=1;
            qy-=1;
        }
        else if(t<7*Qy/4)
        {
            qx-=1;
        }
    }



    return 0;
}

