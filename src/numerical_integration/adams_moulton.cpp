#include"adams_moulton.h"
    
AdamsMoulton::AdamsMoulton(FunctionCapsule & function, labels_values variable, value dt)
:AdamsBashforth(function, variable, dt)  
{
    __N=1;
    RungeKutta model(function, variable, dt);
    step1 = model.get_variable();
    ++model;
    step2 = model.get_variable();
    ++model;
    step3 = model.get_variable();
    ++model;
    step4 = model.get_variable();
    time += 3*dt; 

};


AdamsMoulton& AdamsMoulton::operator++(){
    AdamsBashforth_method();
    step1 = step2;
    step2 = step3;
    step3 = step4;
    step4 = new_step;
    for (unsigned i = 0; i < __N; i++) {
        AdamsMoulton_method();
    };  
    step4 = new_step;
    time  += dt;
    variable = step1;
}

void AdamsMoulton::AdamsMoulton_method() {
    container aux(variable.size(), 0);
    std::vector< container > results;

    function->set(time, step4);
    results.push_back(function->get_result());
    function->set(time, step3);
    results.push_back(function->get_result());
    function->set(time, step2);
    results.push_back(function->get_result());
    function->set(time, step1);
    results.push_back(function->get_result());

    for (unsigned i = 0; i < variable.size(); i++) {

        aux[i] =9 * results[0][i];
        aux[i] += 19 * results[1][i];
        aux[i] -= 5 * results[2][i];
        aux[i] += 1 * results[3][i];
        aux[i] *=(dt / 24);
        aux[i] += step3[i];

        if (aux[i] != aux[i]){
            throw Value_error("Value error in the Adams-Moulton method");
        }
    }
    new_step=aux;
}
