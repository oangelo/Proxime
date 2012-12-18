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
#include "lyapunov.h"
#include "recurrence_relation.h"

using namespace std;

NumericalIntegration* model = NULL;
unsigned step = 1;
unsigned transient = 1000;
unsigned iterations = 1000;
std::vector<double> parameters;
std::vector<double> init;

void help(){
    std::cout << std::endl;
    std::cout << "usage: " << std::endl;
    std::cout << "  ./model [options] [model]" << std::endl;
    std::cout << "models:" << std::endl;
    std::cout << "  --logistic_map" << std::endl;
    std::cout << "  --double_pendulum [model options]" << std::endl;
    std::cout << "  --rossler [model options]" << std::endl;
    std::cout << "model options:" << std::endl;
    std::cout << " --parameters <real>[n]" << std::endl;
    std::cout << " --initial_conditions <real>[n]" << std::endl;
    std::cout << "options:" << std::endl;
    std::cout << "  --step <integer>       number of steps to iterate before print the system properties" << std::endl;
    std::cout << "  --transient <integer>  number of steps to irerate before print anything" << std::endl;
    std::cout << "  --iterations <integer> total number of iterations" << std::endl;
    std::cout << "  --lyapunov             calculates the system lyapunov" << std::endl;
    std::cout << std::endl;
    std::cout << "double pendulum:" << std::endl;
    std::cout << "parameters: l1 l2 m1 m2 g" << std::endl;
    std::cout << "variables: theta1 theta2 omega1 omega2" << std::endl;
    std::cout << std::endl;
}

int main( int argc , char * argv[]) {
    if(argc > 1){
        if(strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0){
            help();
            return 0;
        }
    }else {
        help();
        return 0;
    }

    for (int i = 1; i < argc; ++i)
    {

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
    for (int i = 1; i < argc; ++i)
    {

        if(strcmp(argv[i], "--logistic_map") == 0) {
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

        if(strcmp(argv[i], "--rossler") == 0) {
            std::vector<double> variable(3),parameter(3);
            variable[0] = 2.61622;
            variable[1] = -6.32533;
            variable[2] = 0.0335135;
            parameter[0]= 0.15;
            parameter[1]= 0.2;
            parameter[2]= 10.0;
            model = new AdamsBashforth<RosslerFunction>(variable,parameter,0.0001);
            if(strcmp(argv[i+1], "--lyapunov") == 0) {
                std::string file_name(argv[i + 2]);
                std::cerr << " #>> Calculating Lyapunov exponent" << std::endl;
                lyapunov<Jacobian_RosslerFunction> (*model, iterations, transient, step, file_name);
                return 0;
            }
            for(unsigned i = 0; i < iterations; i++){
                model->next();
                if(i % step == 0)
                    std::cout << *model << std::endl;
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

            for (int j = 1; j < argc; ++j){
                if(strcmp(argv[j], "--parameters") == 0 or strcmp(argv[j], "-p") == 0) {
                    parameter[P_L1]= atof(argv[j + 1]);
                    parameter[P_L2]= atof(argv[j + 2]);
                    parameter[P_M1]= atof(argv[j + 3]);
                    parameter[P_M2]= atof(argv[j + 4]);
                    parameter[P_G]=  atof(argv[j + 5]);
                }
                if(strcmp(argv[j], "--initial_conditions") == 0 or strcmp(argv[j], "-v") == 0) {
                    variable[V_THETA1] = (M_PI / 180) * atof(argv[j + 1]);
                    variable[V_THETA2] = (M_PI / 180) * atof(argv[j + 2]);
                    variable[V_OMEGA1] = atof(argv[j + 3]);
                    variable[V_OMEGA2] = atof(argv[j + 4]);
                }
            }

            model = new AdamsBashforth<DoublePendulumFunction>(variable,parameter,0.00001);
            if(strcmp(argv[i + 1], "--lyapunov") == 0) {
                std::cerr << "#>> Calculating Lyapunov exponent" << std::endl;
                std::string file_name(argv[i + 2]);
                lyapunov<Jacobian_DoublePendulumFunction> (*model, iterations, transient, step, file_name);
                return 0;
            }
            for(unsigned i = 0; i < iterations; i++){
                model->next();
                if(i % step == 0)
                    std::cout << *model << std::endl;
            }
            
        }
    }

    if(model)
        delete model;

    return 0;
}

