#include"adams_moulton.h"
    
AdamsMoulton4Th::AdamsMoulton4Th(FunctionCapsule & function, labels_values variable, value dt, unsigned corrections_amount)
:AdamsBashforth4Th(function, variable, dt), corrections_amount(corrections_amount)  
{
};

AdamsMoulton4Th* AdamsMoulton4Th::Clone() const{
    return new AdamsMoulton4Th(*this);
}

AdamsMoulton4Th* AdamsMoulton4Th::Create(FunctionCapsule & function, labels_values variable,
                                   value dt, unsigned corrections_amount) const{
    return new AdamsMoulton4Th(function, variable, dt, corrections_amount);
}


AdamsMoulton4Th& AdamsMoulton4Th::operator++(){
    AdamsBashforth4Th_method();
    step1 = step2;
    step2 = step3;
    step3 = step4;
    step4 = new_step;

    function_result_1 = function_result_2;
    function_result_2 = function_result_3;
    function_result_3 = function_result_4;
    function->set(time, step4);
    function_result_4 = function->get_result();

    for (unsigned i = 0; i < corrections_amount; i++) {
        AdamsMoulton4Th_method();
    };  
    step4 = new_step;

    time  += dt;
    variable = step1;
    return *this;
}

void AdamsMoulton4Th::AdamsMoulton4Th_method() {
    for (unsigned i = 0; i < variable.size(); i++) {
        new_step[i]  =  9 * function_result_4[i];
        new_step[i] += 19 * function_result_3[i];
        new_step[i] -=  5 * function_result_2[i];
        new_step[i] +=  1 * function_result_1[i];
        new_step[i] *=(dt / 24);
        new_step[i] += step3[i];
        if (new_step[i] != new_step[i]){
            throw Value_error("Value error in the Adams-Moulton method");
        }
    }
}
