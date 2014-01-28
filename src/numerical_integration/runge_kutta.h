#ifndef RUNGE_KUTTA_H
#define RUNGE_KUTTA_H 

#include "numerical_integration.h"

class RungeKutta: public NumericalIntegration{
    public:  
        RungeKutta(functions_capsule & function, container variable,value dt);
        ~RungeKutta();
        virtual void next();
    protected:
        void RungeKutta_method();
};

#endif /* RUNGE_KUTTA_H */
