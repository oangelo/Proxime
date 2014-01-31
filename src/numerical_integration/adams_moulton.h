#ifndef ADAMS_MOULTON_H
#define ADAMS_MOULTON_H 

#include "numerical_integration.h"
#include "adams_bashforth.h"

class AdamsMoulton: public AdamsBashforth {
    public:
        AdamsMoulton(FunctionCapsule & function, labels_values variable, value dt);
        virtual AdamsMoulton& operator++(); //Covariant return type
    private:
        void AdamsMoulton_method();
        unsigned __N;
};



#endif /* ADAMS_MOULTON_H */
