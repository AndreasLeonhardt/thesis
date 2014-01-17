// Staggered Magnetization
// description: calculates the staggered magentization for different U/t
//              and the filling (more like a check). Results are stored in a result file.


#include <iostream>
#include <fstream>
#include <libconfig.h++>
#include <string>
#include <cmath>
#include <stdlib.h>

using namespace std;
using namespace libconfig;


int main()
{
        // config file to get parameters
        Config conf_parameters;
        Config * parameters = &conf_parameters;
        parameters->readFile("/home/andi/Studium/masteroppgave/parameters.cfg");


            double PI = parameters->lookup("PI");
            // sites in one dimension (2D square grid)
            int N_x = parameters->lookup("N_x"); // default 100
            int N_y = parameters->lookup("N_y"); //default 100
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


        // calculate energies \varepsilon_{\vec{k}}
        double epsilon [N_x][N_y];
        for (int i=0;i<N_x;i++)
        {
            for(int j=0;j<N_y;j++)
            {
                epsilon[i][j]= -2*cos(2*PI* (double(i)/N_x-.5))-2*cos(2*PI* (double(j)/N_y-.5));
            }
        }

        // open result file
        ofstream stagmag;
        ofstream edispersion;
        stagmag.open("../stagmagn.txt");
        edispersion.open("../edispersion.txt");


        // loop over different U
        for(U=10.0;U>0.1;U-=.1)
        {
            mu=U/2;
            double E_D [N_x][N_y];
            double E_S [N_x][N_y];

            double fm=0.0;
            double dfdm=0.0;

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


            // calculate n_{\sigma} (relative number of spin \sigma particles (0.5 per spin at half filling)
            // the number is assumed to be 1/2 (using n_u). E^{\pm} depend already on n_{\sigma}
            // so this only checks consistency and accuracy, no optimization.
            // Furthermore it is assumed, that E^+>0 and E^-<0
            double n=.0;
            for (int i=0;i<N_x;i++)
            {
                for (int j=0;j<N_y;j++)
                {
                    n+= -1/E_D[i][j]*( .5*(E_S[i][j]-E_D[i][j]) - epsilon[i][j]+mu-U*n_u );
                }
            }
            // check the difference to 1/2,
            // The result should be a last column basically of zeros (up to double precision)
            n=n/N-.5;

            // write results
            stagmag<<U<<"\t"<<m<<"\t"<<fm<<"\t"<<n<<endl;

        }



        stagmag.close();
        edispersion.close();

    return 0;
}

