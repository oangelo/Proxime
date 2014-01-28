#include"adams_moulton.h"
    
AdamsMoulton::AdamsMoulton(FunctionCapsule & function, container variable, value dt)
:AdamsBashforth(function, variable, dt)  
{
    __N=1;
    //Pont to the class witch encapsulate the functions  
    __func = &function;
    (__func)->set(__t, __variable);
    RungeKutta model(function, variable, dt);
    step1 = model.get_variable();
    model.next();
    step2 = model.get_variable();
    model.next();
    step3 = model.get_variable();
    model.next();
    step4 = model.get_variable();
    __t += 3*__h; 

};

void AdamsMoulton::next() {
    AdamsBashforth_method();
    step1 = step2;
    step2 = step3;
    step3 = step4;
    step4 = new_step;
    for (unsigned i = 0; i < __N; i++) {
        AdamsMoulton_method();
    };  
    step4 = new_step;
    __t  += __h;
    __variable = step1;
};


void AdamsMoulton::AdamsMoulton_method() {
    container aux(__variable.size(), 0);
    std::vector< container > results;

    __func->set(__t, step4);
    results.push_back(__func->get_result());
    __func->set(__t, step3);
    results.push_back(__func->get_result());
    __func->set(__t, step2);
    results.push_back(__func->get_result());
    __func->set(__t, step1);
    results.push_back(__func->get_result());

    for (unsigned i = 0; i < __variable.size(); i++) {

        aux[i] =9 * results[0][i];
        aux[i] += 19 * results[1][i];
        aux[i] -= 5 * results[2][i];
        aux[i] += 1 * results[3][i];
        aux[i] *=(__h / 24);
        aux[i] += step3[i];

        if (aux[i] != aux[i]){
            throw Value_error("Value error in the Adams-Moulton method");
        }
    }
    new_step=aux;
}
