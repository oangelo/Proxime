#include"adams_bashforth.h"

AdamsBashforth::AdamsBashforth(FunctionCapsule & function, labels_values variable, value dt)
:NumericalIntegration(function, variable, dt),
step1(variable.size()), step2(variable.size()), 
step3(variable.size()), step4(variable.size()), 
new_step(variable.size())
{
    /*Ponit to the class witch encapsulate the functions*/  
    RungeKutta model(function, variable, dt);
    step1 = model.get_variable();
    ++model;
    step2 = model.get_variable();
    ++model;
    step3 = model.get_variable();
    ++model;
    step4 = model.get_variable();
}

AdamsBashforth& AdamsBashforth::operator++(){
    AdamsBashforth_method();
    step1 = step2;
    step2 = step3;
    step3 = step4;
    step4 = new_step;
    variable = step1;
    time = time + dt;
    return *this;
}


void AdamsBashforth::AdamsBashforth_method() {
    container aux(variable.size());

    function->set(time, step4);
    container function_step4(function->get_result());

    function->set(time, step3);
    container function_step3(function->get_result());

    function->set(time, step2);
    container function_step2(function->get_result());

    function->set(time, step1);
    container function_step1(function->get_result());

    for (unsigned i = 0; i < variable.size(); i++) {
        aux[i]  = 55 * function_step4[i];
        aux[i] -= 59 * function_step3[i];
        aux[i] += 37 * function_step2[i];
        aux[i] -= 9  * function_step1[i];
        aux[i] *= (dt / 24);
        aux[i] += step4[i];
        if (aux[i] != aux[i]) {
            throw Value_error("Value error in the Adams Bashforth method");
        }
    }
    new_step = aux;
}
