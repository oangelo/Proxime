#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <random>
#include <stdlib.h>

#include <mpi.h>

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
            std::cerr << "#>> transient: " << transient << std::endl; 
        }
        if(strcmp(argv[i], "--iterations") == 0) {
            iterations = atoi(argv[i + 1]);
            std::cerr << "#>> iterations: " << iterations << std::endl; 
        }
    }
    for (size_t i = 1; i < argc; ++i)
    {
 
        if(strcmp(argv[i], "--map") == 0) {
            double a = atof(argv[i + 1]);
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<> dis(0,1);
            logistic_map m(dis(gen), a);
            m.next(transient);
            for(unsigned i = 0; i < iterations; i++){
                m.next(step);
                std::cout << m.get_variable(0) << std::endl;
            }
        }

        if(strcmp(argv[i], "--double_pendulum") == 0 or strcmp(argv[i], "-dp") == 0) {
            type_container variable(4),parameter(5);
            variable[V_THETA1] = M_PI/2.0;
            variable[V_THETA2] = M_PI/2.0;
            variable[V_OMEGA1] = 0.0;
            variable[V_OMEGA2] = 0.0;
            parameter[P_L1]= 0.30;
            parameter[P_L2]= 0.30;
            parameter[P_M1]= 0.10;
            parameter[P_M2]= 0.10;
            parameter[P_G]= 9.8;

            for (size_t j = 1; j < argc; ++j){
                if(strcmp(argv[j], "--parameter") == 0 or strcmp(argv[j], "-p") == 0) {
                    parameter[P_L1]= atof(argv[j + 1]);
                    parameter[P_L2]= atof(argv[j + 2]);
                    parameter[P_M1]= atof(argv[j + 3]);
                    parameter[P_M2]= atof(argv[j + 4]);
                    parameter[P_G]=  atof(argv[j + 5]);
                }
                if(strcmp(argv[j], "--variable") == 0 or strcmp(argv[j], "-v") == 0) {
                    variable[V_THETA1] = atof(argv[j + 1]);
                    variable[V_THETA2] = atof(argv[j + 2]);
                    variable[V_OMEGA1] = atof(argv[j + 3]);
                    variable[V_OMEGA2] = atof(argv[j + 4]);
                }
            }

            AdamsBashforth<pendulum_func> double_pendulum(variable,parameter,0.0001);
            for(unsigned i = 0; i < iterations; i++){
                double_pendulum.next();
                if(i % step == 0)
                    std::cout << double_pendulum << std::endl;
            }
        }



    }

    if(model)
        delete model;

    return 0;
}


