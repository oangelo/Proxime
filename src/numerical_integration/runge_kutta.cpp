#include"runge_kutta.h"

RungeKutta::RungeKutta(FunctionCapsule & function, container variable,value dt)
    :NumericalIntegration(variable, dt) {
        /*Point to the class witch encapsulate the functions*/  
        __func = &function;
    }


RungeKutta::~RungeKutta(){}

void RungeKutta::next(){
  RungeKutta_method();
}

void RungeKutta::RungeKutta_method()
{
    //value k[__variable.size()][4];
    std::vector<container> k(__variable.size(), container(4, 0));
    value t = __t;

    unsigned i;
    container aux_argument(__variable.size(), 0);


    //Generating K1*************************************************************
    __func->set(__t, __variable);
    for (i = 0; i < __variable.size(); i++) {
        k[i][0] = __func->get_result(i);
        //k[i][0] = ((*__function[i])(__t, __variable, __parameter));
    }

    //Generating K2*************************************************************
    t = __t + 0.5 * __h;
    for (i = 0; i < __variable.size(); i++) {
        aux_argument[i] = __variable[i]+(0.5 * __h * k[i][0]);
    }
    __func->set(t, aux_argument);
    for (i = 0; i < __variable.size(); i++) {
        k[i][1] = __func->get_result(i);
        //k[i][1] = ((*__function[i])(t, aux_argument, __parameter));
    }
    //Generating K3*************************************************************
    t = __t + 0.5 * __h;
    for (i = 0; i < __variable.size(); i++) {
        aux_argument[i] = __variable[i]+(0.5 * __h * k[i][1]);
    }
    __func->set(t, aux_argument);
    for (i = 0; i < __variable.size(); i++) {
        k[i][2] = __func->get_result(i);
        //k[i][2] = ((*__function[i])(t, aux_argument, __parameter));
    }
    //Generating K4*************************************************************
    t = __t + __h;
    for (i = 0; i < __variable.size(); i++) {
        aux_argument[i] = __variable[i] + __h * k[i][2];
    }
    __func->set(t, aux_argument);
    for (i = 0; i < __variable.size(); i++) {
        k[i][3] = __func->get_result(i);
        //k[i][3] = ((*__function[i])(t, aux_argument, __parameter));
    }

    //Generating the final values***********************************************
   
    __t = __t + __h;
    for (i = 0; i < __variable.size(); i++) {
        __variable[i] = __variable[i] + __h * (k[i][0] + 2 * k[i][1] + 2 * k[i][2] + k[i][3]) / 6;
        if (__variable[i] != __variable[i]) {
            throw Value_error("Value error in the Runge Kutta method");
        }
    }
   

}


