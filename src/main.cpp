#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <math.h>
#include <string>
#include <stdlib.h>

#include <mpi.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>

#include "bifurcation_diagram.h"
#include "Lyapunov.h"
#include "recurrence_relation.h"

using namespace std;

Numerical_Integration* model = NULL;
unsigned step = 1;
unsigned transient = 1000;
unsigned iterations = 1000;
std::vector<double> parameters;
std::vector<double> init;

int main( int argc , char * argv[]) {
    for (size_t i = 1; i < argc; ++i)
    {
        if(strcmp(argv[i], "--model") == 0)
            if(strcmp(argv[i + 1] ,"rossler") == 0){
                type_container variable(3),parameter(3);
                variable[0] = 2.61622;
                variable[1] = -6.32533;
                variable[2] = 0.0335135;
                parameter[0]= 0.15;
                parameter[1]= 0.2;
                parameter[2]= 10.0;
                model = new AdamsBashforth<rossler_func>(variable,parameter,0.0001);
            } 
        if(strcmp(argv[i], "--step") == 0) {
            step = atoi(argv[i + 1]);
            std::cerr << "#>> step: " << step << std::endl; 
        }
        if(strcmp(argv[i], "--transient") == 0) {
            transient = atoi(argv[i + 1]);
            std::cerr << "#>> step: " << transient << std::endl; 
        }
        if(strcmp(argv[i], "--iterations") == 0) {
            iterations = atoi(argv[i + 1]);
            std::cerr << "#>> iterations: " << iterations << std::endl; 
        }

        if(strcmp(argv[i], "--map") == 0) {
            double a = atof(argv[i + 1]);
            logistic_map m(0.3,a);
            m.next(transient);
            for(unsigned i = 0; i < iterations; i++){
                m.next(step);
                std::cout << m.get_variable(0) << std::endl;
            }
        }


    }

    if(model)
        delete model;

    return 0;
}


