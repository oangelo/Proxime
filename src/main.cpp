#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <cmath>
#include <random>
#include <stdlib.h>

//#include <mpi.h>

#include "functions/double_pendulum.h"
//#include "functions/lorenz.h"
#include "functions/rossler.h"
#include "bifurcation_diagram.h"
//#include "lyapunov.h"
#include "recurrence_relation.h"

using namespace std;

NumericalIntegration* model(NULL);
unsigned step(1);
unsigned transient(1000);
unsigned iterations(1000);
container init;

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
    std::cout << "  --print                          print the attractor on the screen" << std::endl;
    std::cout << "  --step <integer>                 number of steps to iterate before print the system properties" << std::endl;
    std::cout << "  --transient <integer>            number of steps to irerate before print anything" << std::endl;
    std::cout << "  --iterations <integer>           total number of iterations" << std::endl;
    std::cout << "  --lyapunov <string>              calculates the system lyapunov, and print on a file" << std::endl;
    std::cout << "  --max_lyapunov                   calculates the system maximum lyapunov" << std::endl;
    std::cout << "  --bifurcation <option> <integer> X <integer> Y <integer> value <real>" << std::endl;
    std::cout << "                                   option: variable, parameter, X = variable to use as coordinate x, " << std::endl;
    std::cout << "                                   Y = variable to use as coordinate y, " << std::endl;
    std::cout << "                                   value = value were the trajectory corss Y, " << std::endl;
    std::cout << std::endl;
    std::cout << "value pendulum:" << std::endl;
    std::cout << "parameters: l1 l2 m1 m2 g" << std::endl;
    std::cout << "variables: theta1 theta2 omega1 omega2  (in degrees)" << std::endl;
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
            value a = atof(argv[i + 1]);
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
            labels_values variable,parameter;
            variable["x"] = 2.61622;
            variable["y"] = -6.32533;
            variable["z"] = 0.0335135;
            parameter["a"]= 0.15;
            parameter["b"]= 0.2;
            parameter["c"]= 10.0;
            for (int j = 1; j < argc; ++j){
                if(strcmp(argv[j], "--parameters") == 0 or strcmp(argv[j], "-p") == 0) {
                    std::cerr << "#>> setting the parameters" << std::endl;
                    parameter["a"]= atof(argv[j + 1]);
                    parameter["b"]= atof(argv[j + 2]);
                    parameter["c"]= atof(argv[j + 3]);
                }
                if(strcmp(argv[j], "--initial_conditions") == 0 or strcmp(argv[j], "-v") == 0) {
                    std::cerr << "#>> setting the initial conditions " << std::endl;
                    variable["x"] = atof(argv[j + 1]);
                    variable["y"] = atof(argv[j + 2]);
                    variable["z"] = atof(argv[j + 3]);
                }
            }
            RosslerFunction function(parameter);
            model = new AdamsBashforth4Th(function, variable, 0.001);
            for (int j = 1; j < argc; ++j){
                if(strcmp(argv[j], "--lyapunov") == 0) {
                    std::string file_name(argv[i + 1]);
                    std::cerr << " #>> Calculating Lyapunov exponent" << std::endl;
                    for(unsigned i = 0; i < transient; i++){
                        ++model;
                    }
                    //Lyapunov<Jacobian_RosslerFunction> (*model, iterations, 1000, step, file_name);
                    return 0;
                }
                if(strcmp(argv[j], "--lyapunov") == 0) {
                    std::cerr << " #>> Calculating Maximum Lyapunov exponent" << std::endl;
                    for(unsigned i = 0; i < transient; i++){
                        ++model;
                    }
                    //std::cout << MaxLyapunov<Jacobian_RosslerFunction> (*model, iterations, 1000) << std::endl;
                    return 0;
                }

            }
        }

        if(strcmp(argv[i], "--double_pendulum") == 0 or strcmp(argv[i], "-dp") == 0) {

            labels_values variable;
            variable["theta1"] = M_PI / 10;
            variable["theta2"] = M_PI / 10;
            variable["omega1"] = 0.0;
            variable["omega2"] = 0.0;


            labels_values parameter;
            parameter["l1"]= 0.30;
            parameter["l2"]= 0.30;
            parameter["m1"]= 0.1;
            parameter["m2"]= 0.1;
            parameter["g"]= 9.8;


            for (int j = 1; j < argc; ++j){
                if(strcmp(argv[j], "--parameters") == 0 or strcmp(argv[j], "-p") == 0) {
                    std::cerr << "#>> setting the initial parameters" << std::endl;
                    parameter["l1"]= atof(argv[j + 1]);
                    parameter["l2"]= atof(argv[j + 2]);
                    parameter["m1"]= atof(argv[j + 3]);
                    parameter["m2"]= atof(argv[j + 4]);
                    parameter["g"]=  atof(argv[j + 5]);
                    std::cerr << "#>>  l1: " <<  parameter["l1"] << std::endl;
                    std::cerr << "#>>  l2: " <<  parameter["l2"] << std::endl;
                    std::cerr << "#>>  m1: " <<  parameter["m1"] << std::endl;
                    std::cerr << "#>>  m2: " <<  parameter["m2"] << std::endl;
                    std::cerr << "#>>  g:  " <<  parameter["g"] << std::endl;

                }
                if(strcmp(argv[j], "--initial_conditions") == 0 or strcmp(argv[j], "-v") == 0) {
                    std::cerr << "#>> setting the initial conditions " << std::endl;
                    variable["theta1"] = (M_PI / 180.0) * atof(argv[j + 1]);
                    variable["theta2"] = (M_PI / 180.0) * atof(argv[j + 2]);
                    variable["omega1"] = atof(argv[j + 3]);
                    variable["omega2"] = atof(argv[j + 4]);
                    std::cerr << "#>>  thetha1: " <<  variable["theta1"] << " rad" << std::endl;
                    std::cerr << "#>>  thetha2: " <<  variable["theta2"] << " rad" << std::endl;
                    std::cerr << "#>>  Omega1:  " <<  variable["omega1"] << " rad/[t]" << std::endl;
                    std::cerr << "#>>  Omega2:  " <<  variable["omega2"] << " rad/[t]" << std::endl;
                }
            }
            DoublePendulumFunction function(parameter);
            model = new AdamsBashforth4Th(function, variable,  0.00001);
            for (int j = 1; j < argc; ++j){
                if(strcmp(argv[j], "--lyapunov") == 0) {
                    std::cerr << "#>> Calculating Lyapunov exponent" << std::endl;
                    std::string file_name(argv[j + 1]);
                    for(unsigned i = 0; i < transient; i++){
                        ++model;
                    }
                    //Lyapunov<Jacobian_DoublePendulumFunction> (*model, iterations, 100000, step, file_name);
                    return 0;
                }
                if(strcmp(argv[j], "--max_lyapunov") == 0) {
                    std::cerr << "#>> Calculating maximum lyapunov exponent" << std::endl;
                    for(unsigned i = 0; i < transient; i++){
                        ++model;
                    }
                   // std::cout << MaxLyapunov<Jacobian_DoublePendulumFunction> (*model, iterations, 100000) << std::endl;
                    return 0;
                }
            }
        }
    }

    for (int i = 1; i < argc; ++i) {
        
        if(strcmp(argv[i], "--bifurcations") == 0) {
            std::cerr << "#>> Calculating bifurcations" << std::endl;
            value coordinate_value = 0.0;
            int coordinate_x = 0, coordinate_y = 2;

/*            value control_parameter = std::numeric_limits<value>::signaling_NaN();
            for (int j = i; j < argc; ++j) {
                if(strcmp(argv[j], "variable") == 0) 
                    control_parameter = model->get_variable(atoi(argv[i + 1]));
                if(strcmp(argv[j], "parameter") == 0)
                    control_parameter = model->get_parameter(atoi(argv[i + 1]));
                if(strcmp(argv[j], "X") == 0) 
                    coordinate_x = atoi(argv[j + 1]);
                if(strcmp(argv[j], "Y") == 0) 
                    coordinate_y = atoi(argv[j + 1]);
                if(strcmp(argv[j], "value") == 0) 
                    coordinate_value = atoi(argv[j + 1]);
            }

            auto results = PhasePlaneSection(*model, coordinate_x, coordinate_y, coordinate_value);
            for(auto i: results)
                std::cout << control_parameter << " " << i << std::endl;

                    */
        }
        if(strcmp(argv[i], "--print") == 0) {
                for(unsigned i = 0; i < transient; i++){
                    ++(*model);
                }
                for(unsigned i = 0; i < iterations; i++){
                    ++(*model);
                    if(i % step == 0)
                        std::cout << *model << std::endl;
                }
        }

    }

    if(model)
        delete model;

    return 0;
}

