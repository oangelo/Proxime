#ifndef ADAMS_MOULTON_H
#define ADAMS_MOULTON_H 

#include "numerical_integration.h"
#include "adams_bashforth.h"

class AdamsMoulton: public AdamsBashforth {
    public:
        AdamsMoulton(FunctionCapsule & function, container variable, value dt);
        virtual void next();

    private:
        void AdamsMoulton_method();
        unsigned __N;
};



#endif /* ADAMS_MOULTON_H */
