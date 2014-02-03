#ifndef RUNGE_KUTTA_H
#define RUNGE_KUTTA_H 

#include "numerical_integration.h"

class RungeKutta4Th: public NumericalIntegration{
    public:  
        RungeKutta4Th(FunctionCapsule & function, labels_values variable,value dt);
        virtual RungeKutta4Th& operator++(); //Covariant return type
        virtual RungeKutta4Th* Clone() const;
        virtual RungeKutta4Th* Create(FunctionCapsule & function, labels_values variable,value dt) const;
    protected:
        void RungeKutta4Th_method();
};

#endif /* RUNGE_KUTTA_H */
