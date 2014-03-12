#ifndef ADAMS_BASHFORTH_H
#define ADAMS_BASHFORTH_H 

#include "numerical_integration.h"
#include "runge_kutta.h"


class AdamsBashforth4Th: public NumericalIntegration{
    public:
        AdamsBashforth4Th(FunctionCapsule & function, labels_values variable, value dt);
        virtual AdamsBashforth4Th& operator++(); //Covariant return type
        virtual AdamsBashforth4Th* Clone() const;
        virtual AdamsBashforth4Th* Create(FunctionCapsule & function, labels_values variable,value dt) const;
        virtual ~AdamsBashforth4Th(){};
    protected:
        void AdamsBashforth4Th_method();
        container step1, step2, step3, step4, new_step;
        container function_result_1, function_result_2, function_result_3, function_result_4;
};

#endif /* ADAMS_BASHFORTH_H */
