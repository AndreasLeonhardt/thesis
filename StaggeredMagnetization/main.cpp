// Staggered Magnetization
// description: calculates the staggered magentization for different U/t
//              and the filling (more like a check). Results are stored in a result file.


#include <iostream>
#include <fstream>
#include <libconfig.h++>
#include <string>
#include <cmath>
#include <stdlib.h>

#include <../Susceptibility/propagator.h>

using namespace std;
using namespace libconfig;


int main()
{
        // config file to get parameters
        Config conf_parameters;
        Config * parameters = &conf_parameters;
        parameters->readFile("/home/andi/Studium/masteroppgave/parameters.cfg");


        // create instance of class propagator, setting energies, magnetization and so on
        propagator Prop(parameters);
        propagator * prop = &Prop;

        // open result file
        ofstream stagmag;
        stagmag.open(parameters->lookup("stag"));


        // loop over different U
        for(double U=10.0;U>=0.0;U-=.1)
        {
            prop->set_U(U,U/2.0);

            // write results
            stagmag<<U<<"\t"<<prop->get_m()<<"\t"<<"\t"<<prop->calc_n()<<endl;

        }

        stagmag.close();

    return 0;
}

