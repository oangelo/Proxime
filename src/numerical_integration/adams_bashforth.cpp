#include"adams_bashforth.h"

AdamsBashforth::AdamsBashforth(FunctionCapsule & function, labels_values variable, value dt)
:NumericalIntegration(function, variable, dt),
step1(variable.size()), step2(variable.size()), 
step3(variable.size()), step4(variable.size()), 
new_step(variable.size()),
function_result_1(variable.size()), function_result_2(variable.size()), 
function_result_3(variable.size()), function_result_4(variable.size())
{
    //Starting the method with Runge-Kutta
    RungeKutta model(function, variable, dt);
    step1 = model.get_variable();
    ++model;
    step2 = model.get_variable();
    ++model;
    step3 = model.get_variable();
    ++model;
    step4 = model.get_variable();

    //This will be used to calculate the next step
    this->function->set(time, step4);
    function_result_4 = this->function->get_result();

    this->function->set(time, step3);
    function_result_3 = this->function->get_result();

    this->function->set(time, step2);
    function_result_2 = this->function->get_result();

    this->function->set(time, step1);
    function_result_1 = this->function->get_result();
}

AdamsBashforth* AdamsBashforth::Clone() const{
    return new AdamsBashforth(*this);
}

AdamsBashforth* AdamsBashforth::Create(FunctionCapsule & function, labels_values variable,value dt) const{
    return new AdamsBashforth(function, variable, dt);
}

AdamsBashforth& AdamsBashforth::operator++(){
    AdamsBashforth_method();

    step1 = step2;
    step2 = step3;
    step3 = step4;
    step4 = new_step;

    function_result_1 = function_result_2;
    function_result_2 = function_result_3;
    function_result_3 = function_result_4;
    function->set(time, step4);
    function_result_4 = function->get_result();

    variable = step1;
    time = time + dt;
    return *this;
}

void AdamsBashforth::AdamsBashforth_method() {
    for (unsigned i = 0; i < variable.size(); i++) {
        new_step[i]  = 55 * function_result_4[i];
        new_step[i] -= 59 * function_result_3[i];
        new_step[i] += 37 * function_result_2[i];
        new_step[i] -= 9  * function_result_1[i];
        new_step[i] *= (dt / 24);
        new_step[i] += step4[i];
        if (new_step[i] != new_step[i]) {
            throw Value_error("Value error in the Adams Bashforth method");
        }
    }
}
