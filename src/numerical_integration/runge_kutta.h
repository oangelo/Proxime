#ifndef RUNGE_KUTTA_H
#define RUNGE_KUTTA_H 

#include "numerical_integration.h"

class RungeKutta: public NumericalIntegration{
    public:  
        RungeKutta(FunctionCapsule & function, labels_values variable,value dt);
        ~RungeKutta();
        virtual RungeKutta& operator++(); //Covariant return type
    protected:
        void RungeKutta_method();
};

#endif /* RUNGE_KUTTA_H */
