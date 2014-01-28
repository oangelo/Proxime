#ifndef RUNGE_KUTTA_H
#define RUNGE_KUTTA_H 

#include "numerical_integration.h"

class RungeKutta: public NumericalIntegration{
    public:  
        RungeKutta(FunctionCapsule & function, container variable,value dt);
        ~RungeKutta();
        virtual void next();
        virtual NumericalIntegration& operator++();
    protected:
        void RungeKutta_method();
};

#endif /* RUNGE_KUTTA_H */
