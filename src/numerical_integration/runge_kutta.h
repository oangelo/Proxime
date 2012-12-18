#ifndef RUNGE_KUTTA_H
#define RUNGE_KUTTA_H 

#include "numerical_integration.h"

template <class function>
class RungeKutta: public NumericalIntegration{
    public:  
        RungeKutta(type_container variable,type_container parameter,type_data dt)
            :NumericalIntegration(variable,parameter,dt) {
                /*Point to the class witch encapsulate the functions*/  
                __func.reset(new function);
                (*__func).set(__t, __variable, __parameter);   
            }
        ~RungeKutta();
        virtual void next();
    protected:
        void RungeKutta_method();
};

template <class function>
RungeKutta<function>::~RungeKutta(){}

template <class function>
void RungeKutta<function>::next(){
  RungeKutta_method();
}

template <class function>
void RungeKutta<function>::RungeKutta_method()
{
    //type_data k[__variable.size()][4];
    std::vector<std::vector<double>> k(__variable.size(), std::vector<double>(4, 0));
    type_data t = __t;

    unsigned i;
    type_container aux_argument(__variable.size(), 0);


    //Generating K1*************************************************************
    __func->set(__t, __variable, __parameter);
    for (i = 0; i < __variable.size(); i++) {
        k[i][0] = __func->get_result(i);
        //k[i][0] = ((*__function[i])(__t, __variable, __parameter));
    }

    //Generating K2*************************************************************
    t = __t + 0.5 * __h;
    for (i = 0; i < __variable.size(); i++) {
        aux_argument[i] = __variable[i]+(0.5 * __h * k[i][0]);
    }
    __func->set(t, aux_argument, __parameter);
    for (i = 0; i < __variable.size(); i++) {
        k[i][1] = __func->get_result(i);
        //k[i][1] = ((*__function[i])(t, aux_argument, __parameter));
    }
    //Generating K3*************************************************************
    t = __t + 0.5 * __h;
    for (i = 0; i < __variable.size(); i++) {
        aux_argument[i] = __variable[i]+(0.5 * __h * k[i][1]);
    }
    __func->set(t, aux_argument, __parameter);
    for (i = 0; i < __variable.size(); i++) {
        k[i][2] = __func->get_result(i);
        //k[i][2] = ((*__function[i])(t, aux_argument, __parameter));
    }
    //Generating K4*************************************************************
    t = __t + __h;
    for (i = 0; i < __variable.size(); i++) {
        aux_argument[i] = __variable[i] + __h * k[i][2];
    }
    __func->set(t, aux_argument, __parameter);
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

#endif /* RUNGE_KUTTA_H */
