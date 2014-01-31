#ifndef ADAMS_BASHFORTH_H
#define ADAMS_BASHFORTH_H 

#include "numerical_integration.h"
#include "runge_kutta.h"

class AdamsBashforth: public NumericalIntegration{
    public:
        AdamsBashforth(FunctionCapsule & function, labels_values variable, value dt);
        virtual AdamsBashforth& operator++(); //Covariant return type
        virtual ~AdamsBashforth(){};
    protected:
        void AdamsBashforth_method();
        container step1, step2, step3, step4, new_step;
};

#endif /* ADAMS_BASHFORTH_H */
