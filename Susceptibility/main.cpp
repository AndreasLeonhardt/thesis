// Susceptibility
// description: Calculates first the staggered magnetization according to the given parameters U,t,\beta.
//              Proceeds in calculating \Xi^{+-}(q,\omega) (later \Xi^{zz} as well)
//              and writes the result in a file.

// test comment
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <stdlib.h>
#include <cstring>
#include <libconfig.h++>



using namespace std;
using namespace libconfig;

int main()
{


    // config file to get parameters
    Config conf_parameters;
    Config * parameters = &conf_parameters;
    conf_parameters.readFile("/home/andi/Studium/masteroppgave/parameters.cfg");

        double PI = parameters->lookup("PI");
        // sites in one dimension (2D square grid)
        int N_x = parameters->lookup("N_x"); // default 100
        int N_xh = N_x/2;
        int N_y = parameters->lookup("N_y"); //default 100
        int N_yh = N_y/2;
        // number of sites, must be even
        int N = N_x*N_y;
        // here U corresponds to \frac{U}{t} with t normalized to 1
        double U = parameters->lookup("U"); // default 2.0
        // initial value for up-/down- part of the staggered magnetization m_{s,\sigma}, here m,
        // without distinguishing the up and down case, e.g. m_u = m_d = m = \frac12 m_{\mathrm{s}}
        double m = parameters->lookup("m_initial"); // default .45
        // filling, n \in [0,2]; n_u,n_d \in [0,1]
        double n_u= parameters->lookup("n_u"); // default .5
        double n_d= parameters->lookup("n_d"); // default .5
        // chemical potential mu (half filling at mu=U/2?) n should be calculated from mu or the other way round.
        double mu = parameters->lookup("mu"); // default .5*U

        double T = parameters->lookup("T");
        double beta = -1.0; // to be set below
        if (T==0.0)
        {
            // proceed differently, e.g. don't use beta
        }
        else
        {
            beta = 1.0/T;
        }


        // calculate energies \varepsilon_{\vec{k}} of the unpertubed system
        double epsilon [N_x][N_y];
        for (int i=0;i<N_x;i++)
        {
            for(int j=0;j<N_y;j++)
            {
                epsilon[i][j]= -2*cos(2*PI* (double(i)/N_x-.5))-2*cos(2*PI* (double(j)/N_y-.5));
            }
        }




//-------------------------------- optimization of m --------------------------------------------------------
        double E_D [N_x][N_y];
        double E_S [N_x][N_y];

        double fm;
        double dfdm;

        // inverse f(m)=1; (Newton Raphson, weighted)
        int abort = 0;
        double mdiff=0.25;
        while (mdiff>1e-7)
        {

            // calculate E_D and E_S ( E_S=E^{+} + E^{-} ; E_D = E^{+} - E^{-} )

            for (int i=0;i<N_x;i++)
            {
                for (int j=0;j<N_y;j++)
                {
                    E_S[i][j]=(epsilon[i][j]+epsilon[(i+N_x/2)%N_x][(j+N_y/2)%N_y]) -2*mu +2*U*n_u;
                    E_D[i][j]=2.0*sqrt( .25*pow((epsilon[i][j]-epsilon[(i+N_x/2)%N_x][(j+N_y/2)%N_y]),2) + U*U*m*m );
                }
            }


            // calculate m_u and m_d
            fm=0.0;
            dfdm=0.0;
            for (int i=0;i<N_x;i++)
            {
                for (int j=0;j<N_y;j++)
                {
                    if (abs(E_D[i][j])>abs(E_S[i][j]))
                    {
                        fm+=1.0/E_D[i][j];
                        dfdm+=1.0/pow(E_D[i][j],3);
                    }
                    else
                    {
                        // skip
                    }
                }
            }


            // calculate new m
            fm*=U/N;
            dfdm*=-4*U*U*U*m/N;


            // next m with Newton-Raphson
            // check difference, with infinite recursive filter, to not only rely on a single lucky jump
            mdiff=0.1*mdiff+abs((fm-1)/dfdm);
            // set new m, weighted with A/(A+abort) to avoid jumping forth and back
            m-=(fm-1)/dfdm*1.0e4/(1.0e4+abort);
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


        cout<<"m="<<m<<endl;

//------------------------- end of optimization of m --------------------------------------------------------





        // calculate mean-field band energies
        double E_p [N_xh][N_yh];
        double E_m [N_xh][N_yh];
        double E_Sh, E_Sd;
            for (int i=0;i<N_xh;i++)
            {
                for (int j=0;j<N_yh;j++)
                {
                    E_Sh = .5*(epsilon[i][j]+epsilon[(i+N_xh)][(j+N_yh)]) -mu +U*n_u;
                    E_Sd = sqrt( .25*pow((epsilon[i][j]-epsilon[i+N_xh][j+N_yh]),2) + U*U*m*m );

                    E_p[i][j] = E_Sh+E_Sd;
                    E_m[i][j] = E_Sh-E_Sd;
                    // cout<<"E_p="<<E_p[i][j]<<"\tE_m="<<E_m[i][j]<<endl;
                }
            }



        // set q ( from (0,0) to (N_x,N_y) corresponding to to be in (-\pi,\pi)\times (-\pi,\pi)
        // and w, which should be choose in a clever way to give a meaningful result
            int q_x,q_y;
            int t_steps=100;
            for (int t =0;t<t_steps;t++)
            {
                for(double w=0.0;w<50.0;w+=1.0)
                {
                    q_x = N_x*t/t_steps;
                    q_y = N_y*t/t_steps;

                    double x=0.0,xb=0.0, y=0.0, z1=0.0, z2=0.0;

                    // calculate x(q,\omega), y, z_1,z_2 (\bar{x}(q,\omega) = x(q,\omega)
                    if (T!=0.0)
                    {
                        for (int i=0;i<N_xh;i++)
                        {
                            for (int j=0;j<N_yh;j++)
                            {
                                x +=     (E_p[i][j]-epsilon[i+N_xh][j+N_yh])
                                        *(E_p[i][j]+w-epsilon[(i+q_x+N_xh)%N_x][(j+q_y+N_yh)%N_y])
                                        /(E_p[i][j]-E_m[i][j])
                                        /(E_p[i][j]+w-E_p[(i+q_x)%N_xh][(j+q_y)%N_yh])
                                        /(E_p[i][j]+w-E_m[(i+q_x)%N_xh][(j+q_y)%N_yh])
                                        /(1+exp(beta*E_p[i][j]));

                                x+=      (E_m[i][j]-epsilon[i+N_xh][j+N_yh])
                                        *(E_m[i][j]+w-epsilon[(i+q_x+N_xh)%N_x][(j+q_y+N_yh)%N_y])
                                        /(E_m[i][j]-E_p[i][j])
                                        /(E_m[i][j]+w-E_p[(i+q_x)%N_xh][(j+q_y)%N_yh])
                                        /(E_m[i][j]+w-E_m[(i+q_x)%N_xh][(j+q_y)%N_yh])
                                        /(1+exp(beta*E_m[i][j]));

                                x+=      ( E_p[(i+q_x)%N_xh][(j+q_y)%N_yh] -w -epsilon[i+N_xh][j+N_yh] )
                                        *( E_p[(i+q_x)%N_xh][(j+q_y)%N_yh] -epsilon[(i+q_x+N_xh)%N_x][(j+q_y+N_yh)%N_y] )
                                        /( E_p[(i+q_x)%N_xh][(j+q_y)%N_yh] -w -E_p[i][j] )
                                        /( E_p[(i+q_x)%N_xh][(j+q_y)%N_yh] -w -E_m[i][j] )
                                        /( E_p[(i+q_x)%N_xh][(j+q_y)%N_yh] -E_m[(i+q_x)%N_xh][(j+q_y)%N_yh] )
                                        /( 1 + exp(beta* E_p[(i+q_x)%N_xh][(j+q_y)%N_yh] ));

                                x+=      ( E_m[(i+q_x)%N_xh][(j+q_y)%N_yh] -w -epsilon[i+N_xh][j+N_yh] )
                                        *( E_m[(i+q_x)%N_xh][(j+q_y)%N_yh] -epsilon[(i+q_x+N_xh)%N_x][(j+q_y+N_yh)%N_y] )
                                        /( E_m[(i+q_x)%N_xh][(j+q_y)%N_yh] -w -E_p[i][j] )
                                        /( E_m[(i+q_x)%N_xh][(j+q_y)%N_yh] -w -E_m[i][j] )
                                        /( E_m[(i+q_x)%N_xh][(j+q_y)%N_yh] -E_p[(i+q_x)%N_xh][(j+q_y)%N_yh] )
                                        /( 1 + exp(beta* E_m[(i+q_x)%N_xh][(j+q_y)%N_yh] ));

                            }
                        }
                    }
                    else
                    {
                        // same but for T=0, e.g. sharp fermionic distribution (yaigh, only half the terms)
                        for (int i=0;i<N_xh;i++)
                        {
                            for (int j=0;j<N_yh;j++)
                            {
                                x+=      (E_m[i][j]-epsilon[i+N_xh][j+N_yh])
                                        *(E_m[i][j]+w-epsilon[(i+q_x+N_xh)%N_x][(j+q_y+N_yh)%N_y])
                                        /(E_m[i][j]-E_p[i][j])
                                        /(E_m[i][j]+w-E_p[(i+q_x)%N_xh][(j+q_y)%N_yh])
                                        /(E_m[i][j]+w-E_m[(i+q_x)%N_xh][(j+q_y)%N_yh]);

                                x+=      ( E_m[(i+q_x)%N_xh][(j+q_y)%N_yh] -w -epsilon[i+N_xh][j+N_yh] )
                                        *( E_m[(i+q_x)%N_xh][(j+q_y)%N_yh] -epsilon[(i+q_x+N_xh)%N_x][(j+q_y+N_yh)%N_y] )
                                        /( E_m[(i+q_x)%N_xh][(j+q_y)%N_yh] -w -E_p[i][j] )
                                        /( E_m[(i+q_x)%N_xh][(j+q_y)%N_yh] -w -E_m[i][j] )
                                        /( E_m[(i+q_x)%N_xh][(j+q_y)%N_yh] -E_p[(i+q_x)%N_xh][(j+q_y)%N_yh] );
                            }
                        }
                    }
                    cout <<"q=("<<q_x*2*PI/N_x-PI<<","<<q_y*2*PI/N_y-PI<<")\tw="<<w<<"\tx="<<x<<endl;

                }
            }



    return 0;
}

