#ifndef RUNGE_KUTTA_H
#define RUNGE_KUTTA_H 

#include "numerical_integration.h"

class RungeKutta: public NumericalIntegration{
    public:  
        RungeKutta(FunctionCapsule & function, labels_values variable,value dt);
        virtual RungeKutta& operator++(); //Covariant return type
        virtual RungeKutta* Clone() const;
        virtual RungeKutta* Create(FunctionCapsule & function, labels_values variable,value dt) const;
    protected:
        void RungeKutta_method();
};

#endif /* RUNGE_KUTTA_H */
