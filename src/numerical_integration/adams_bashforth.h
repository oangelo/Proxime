#ifndef ADAMS_BASHFORTH_H
#define ADAMS_BASHFORTH_H 

#include "numerical_integration.h"
#include "runge_kutta.h"

template <class function>
class AdamsBashforth: public NumericalIntegration{
    public:
        AdamsBashforth(type_container variable,type_container parameter,type_data dt)
            :NumericalIntegration(variable,parameter,dt),
            step1(variable.size()), step2(variable.size()), step3(variable.size()), step4(variable.size()), new_step(variable.size())

    {
        /*Ponit to the class witch encapsulate the functions*/  
        this->__func.reset(new function);
        RungeKutta<function> model(variable,parameter,dt);
        step1 = model.get_variable();
        model.next();
        step2 = model.get_variable();
        model.next();
        step3 = model.get_variable();
        model.next();
        step4 = model.get_variable();
    }

        virtual void next(){
            AdamsBashforth_method();
            step1 = step2;
            step2 = step3;
            step3 = step4;
            step4 = new_step;
            __variable = step1;
            __t = __t + __h;
        }

        virtual ~AdamsBashforth(){};
    protected:
        void AdamsBashforth_method();
        type_container step1, step2, step3, step4, new_step;
};

template <class function>
void AdamsBashforth<function>::AdamsBashforth_method() {
    type_container aux(__variable.size());

    __func->set(__t, step4, __parameter);
    std::vector<double> function_step4(__func->get_result());

    __func->set(__t, step3, __parameter);
    std::vector<double> function_step3(__func->get_result());

    __func->set(__t, step2, __parameter);
    std::vector<double> function_step2(__func->get_result());

    __func->set(__t, step1, __parameter);
    std::vector<double> function_step1(__func->get_result());

    for (unsigned i = 0; i < __variable.size(); i++) {
        aux[i]  = 55 * function_step4[i];
        aux[i] -= 59 * function_step3[i];
        aux[i] += 37 * function_step2[i];
        aux[i] -= 9  * function_step1[i];
        aux[i] *= (__h / 24);
        aux[i] += step4[i];
        if (aux[i] != aux[i]) {
            throw Value_error("Value error in the Adams Bashforth method");
        }
    }
    new_step = aux;
}

#endif /* ADAMS_BASHFORTH_H */
