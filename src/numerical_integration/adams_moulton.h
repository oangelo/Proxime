#ifndef ADAMS_MOULTON_H
#define ADAMS_MOULTON_H 

#include "numerical_integration.h"
#include "adams_bashforth.h"


class AdamsMoulton4Th: public AdamsBashforth4Th{
    public:
        AdamsMoulton4Th(FunctionCapsule & function, labels_values variable, value dt, unsigned corrections_amount = 1);
        virtual AdamsMoulton4Th& operator++(); //Covariant return type
        virtual AdamsMoulton4Th* Clone() const;
        virtual AdamsMoulton4Th* Create(FunctionCapsule & function, labels_values variable, value dt, unsigned corrections_amount = 1) const;
    private:
        void AdamsMoulton4Th_method();
        unsigned corrections_amount;
};

#endif /* ADAMS_MOULTON_H */
