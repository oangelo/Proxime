#ifndef ADAMS_MOULTON_H
#define ADAMS_MOULTON_H 

#include "numerical_integration.h"
#include "adams_bashforth.h"

template <class function>
class AdamsMoulton: public AdamsBashforth<function> {
    public:
        AdamsMoulton(type_container variable, type_container parameter, type_data dt)
            :AdamsBashforth<function>(variable, parameter, dt)  
        {
            __N=1;
            //Pont to the class witch encapsulate the functions  
            this->__func.reset(new function);
            (this->__func)->set(this->__t, this->__variable, this->__parameter);
            RungeKutta<function> model(variable,parameter,dt);
            this->step1 = model.get_variable();
            model.next();
            this->step2 = model.get_variable();
            model.next();
            this->step3 = model.get_variable();
            model.next();
            this->step4 = model.get_variable();
            this->__t += 3*this->__h; 

        };
        virtual void next() {
            this->AdamsBashforth_method();
            this->step1 = this->step2;
            this->step2 = this->step3;
            this->step3 = this->step4;
            this->step4 = this->new_step;
            for (unsigned i = 0; i < __N; i++) {
                AdamsMoulton_method();
            };  
            this->step4 = this->new_step;
            this->__t  += this->__h;
            this->__variable = this->step1;
        };

    private:
        void AdamsMoulton_method();
        unsigned __N;
};

template <class function>
void AdamsMoulton<function>::AdamsMoulton_method() {
    type_container aux(this->__variable.size(), 0);
    std::vector< type_container > results;

    this->__func->set(this->__t, this->step4, this->__parameter);
    results.push_back(this->__func->get_result());
    this->__func->set(this->__t, this->step3, this->__parameter);
    results.push_back(this->__func->get_result());
    this->__func->set(this->__t, this->step2, this->__parameter);
    results.push_back(this->__func->get_result());
    this->__func->set(this->__t, this->step1, this->__parameter);
    results.push_back(this->__func->get_result());

    for (unsigned i = 0; i < this->__variable.size(); i++) {
        
        aux[i] =9 * results[0][i];
        aux[i] += 19 * results[1][i];
        aux[i] -= 5 * results[2][i];
        aux[i] += 1 * results[3][i];
        aux[i] *=(this->__h / 24);
        aux[i] += this->step3[i];

        if (aux[i] != aux[i]){
            throw Value_error("Value error in the Adams-Moulton method");
        }
    }
    this->new_step=aux;
}

#endif /* ADAMS_MOULTON_H */
