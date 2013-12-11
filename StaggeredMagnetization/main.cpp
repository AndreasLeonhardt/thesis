#include <iostream>
#include <fstream>
#include <libconfig.h++>
#include <string>
#include <cmath>
#include <stdlib.h>

using namespace std;

int main()
{
    // settings (might be stored in parameter file later)
        double PI = 3.1415926535897932384626433832795028841971693993751058;
        // sites in one dimension (2D square grid)
        int N_x=100;
        int N_y=100;
        // number of sites, must be even
        int N = N_x*N_y;
        // intial value for U, here U corresponds to \frac{U}{t} with t normalized to 1
        double U = 2.0;
        // initial value for staggered magnetization m_{s,\sigma}, here m_u and m_d
        double m_u = .25;
        double m_d = .25;
        // filling, n \in [0,2]; n_u,n_d \in [0,1]
        double n_u=.5,n_d = .5;
        double n = n_u+n_d;
        // chemical potential mu (half filling at mu=0 or mu=U/2?) n should be calculated from mu.
        double mu = .5*U;


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

        for (m_u = 0.001;m_u<0.501;m_u+=0.01)
        {

            // calculate E_Du and E_Su ( E_S=E^{+} + E^{-} ; E_D = E^{+} - E^{-} )
            //double E_Du [N_x/2][N_y/2];
            //double E_Su [N_x/2][N_y/2];
            double E_Dd [N_x][N_y];
            double E_Sd [N_x][N_y];
            for (int i=0;i<N_x;i++)
            {
                for (int j=0;j<N_y;j++)
                {
                    //E_Su[i][j]=(epsilon[i][j]+epsilon[i+N_x/2][j+N_y/2]) -2*mu +2*U*n_d;
                    //E_Du[i][j]=2*sqrt( .25*pow((epsilon[i][j]-epsilon[i+N_x/2][j+N_y/2]),2) + U*U*m_d*m_d );
                    E_Sd[i][j]=(epsilon[i][j]+epsilon[(i+N_x/2)%N_x][(j+N_y/2)%N_y]) -2*mu +2*U*n_u;
                    E_Dd[i][j]=2.0*sqrt( .25*pow((epsilon[i][j]-epsilon[(i+N_x/2)%N_x][(j+N_y/2)%N_y]),2) + U*U*m_u*m_u );
                    //edispersion << m_u<<"\t"<< 2*PI*(double(i)/N_x-.5) <<"\t"<< 2*PI*(double(j)/N_y-.5) << "\t"
                    //            << .5*(E_Sd[i][j]-E_Dd[i][j])<<"\t"<<.5*(E_Sd[i][j]+E_Dd[i][j])<<endl;

                }
            }


            // calculate m_u and m_d
            double fmu=0.0,fmd=0.0;
            for (int i=0;i<N_x;i++)
            {
                for (int j=0;j<N_y;j++)
                {
                   /*
                    if (abs(E_Du[i][j])>abs(2*E_F-E_Su[i][j]))
                    {
                        fmu+=1/E_Du[i][j];
                    }


                    else
                    {
                        // skip
                    }
                    */

                    if (abs(E_Dd[i][j])>abs(E_Sd[i][j]))
                    {
                        fmd+=1.0/E_Dd[i][j];
                    }
                    else
                    {
                        // skip
                    }
                }
            }


            // calculate new m_u and m_d (HOW TO DO THAT REASONABLY?)
            //fmu*=U/N;
            fmd*=U/N;

            // show results
            stagmag<<m_u<<"\t"<<fmd<<endl;

        }

        stagmag.close();
        edispersion.close();
        system("python ../plot.py ../stagmagn.txt");

    return 0;
}

