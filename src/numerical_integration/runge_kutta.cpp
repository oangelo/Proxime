#include"runge_kutta.h"

RungeKutta::RungeKutta(FunctionCapsule & function, container variable,value dt)
    :NumericalIntegration(function, variable, dt) { } 

RungeKutta::~RungeKutta(){}


NumericalIntegration& RungeKutta::operator++(){
    RungeKutta_method();
}

void RungeKutta::RungeKutta_method()
{
    //value k[variable.size()][4];
    std::vector<container> k(variable.size(), container(4, 0));
    value t = time;

    unsigned i;
    container aux_argument(variable.size(), 0);


    //Generating K1*************************************************************
    function->set(time, variable);
    for (i = 0; i < variable.size(); i++) {
        k[i][0] = function->get_result(i);
        //k[i][0] = ((*functiontion[i])(time, variable, __parameter));
    }

    //Generating K2*************************************************************
    t = time + 0.5 * dt;
    for (i = 0; i < variable.size(); i++) {
        aux_argument[i] = variable[i]+(0.5 * dt * k[i][0]);
    }
    function->set(t, aux_argument);
    for (i = 0; i < variable.size(); i++) {
        k[i][1] = function->get_result(i);
        //k[i][1] = ((*functiontion[i])(t, aux_argument, __parameter));
    }
    //Generating K3*************************************************************
    t = time + 0.5 * dt;
    for (i = 0; i < variable.size(); i++) {
        aux_argument[i] = variable[i]+(0.5 * dt * k[i][1]);
    }
    function->set(t, aux_argument);
    for (i = 0; i < variable.size(); i++) {
        k[i][2] = function->get_result(i);
        //k[i][2] = ((*functiontion[i])(t, aux_argument, __parameter));
    }
    //Generating K4*************************************************************
    t = time + dt;
    for (i = 0; i < variable.size(); i++) {
        aux_argument[i] = variable[i] + dt * k[i][2];
    }
    function->set(t, aux_argument);
    for (i = 0; i < variable.size(); i++) {
        k[i][3] = function->get_result(i);
        //k[i][3] = ((*functiontion[i])(t, aux_argument, __parameter));
    }

    //Generating the final values***********************************************
   
    time = time + dt;
    for (i = 0; i < variable.size(); i++) {
        variable[i] = variable[i] + dt * (k[i][0] + 2 * k[i][1] + 2 * k[i][2] + k[i][3]) / 6;
        if (variable[i] != variable[i]) {
            throw Value_error("Value error in the Runge Kutta method");
        }
    }
   

}


