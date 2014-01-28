#include"adams_bashforth.h"

AdamsBashforth::AdamsBashforth(functions_capsule & function, container variable, value dt)
:NumericalIntegration(variable, dt),
step1(variable.size()), step2(variable.size()), 
step3(variable.size()), step4(variable.size()), 
new_step(variable.size())
{
    /*Ponit to the class witch encapsulate the functions*/  
    this->__func = &function;
    RungeKutta model(function, variable, dt);
    step1 = model.get_variable();
    model.next();
    step2 = model.get_variable();
    model.next();
    step3 = model.get_variable();
    model.next();
    step4 = model.get_variable();
}

void AdamsBashforth::next(){
    AdamsBashforth_method();
    step1 = step2;
    step2 = step3;
    step3 = step4;
    step4 = new_step;
    __variable = step1;
    __t = __t + __h;
}


void AdamsBashforth::AdamsBashforth_method() {
    container aux(__variable.size());

    __func->set(__t, step4);
    container function_step4(__func->get_result());

    __func->set(__t, step3);
    container function_step3(__func->get_result());

    __func->set(__t, step2);
    container function_step2(__func->get_result());

    __func->set(__t, step1);
    container function_step1(__func->get_result());

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
