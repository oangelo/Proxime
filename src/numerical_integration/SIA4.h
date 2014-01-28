#ifndef SIA4_H
#define SIA4_H 

#include "numerical_integration.h"
#include "runge_kutta.h"
#include <iostream>

template <class function>
class SIA4: public NumericalIntegration{
    public:
        SIA4(container variable, container parameter, value dt)
            :NumericalIntegration(variable, parameter, dt), q(), p(), a(), b()
    {
        //Initializing the variables
        for(size_t i = 0; i < variable.size() / 2; ++i){
            q.push_back(variable[i]);
        }
        for(size_t i = variable.size() / 2; i < variable.size(); ++i){
            p.push_back(variable[i]);
        }
        this->function.reset(new function);
        function->set(time, q, __parameter);

        value a1, a2, a3, a4;
        value b1, b2, b3, b4;
        a1 = a4 = (1.0 / 6.0) * (2.0 + pow(2.0, 1.0 / 3.0) + pow(2.0,  -1.0 / 3.0));
        a2 = a3 = (1.0 / 6.0) * (1.0 - pow(2.0, 1.0 / 3.0) - pow(2.0,  -1.0 / 3.0));
        b1 = 0;
        b2 = b4 = 1 / (2.0 - pow(2.0, 1.0 / 3.0));
        b3 = 1 / (1.0 - pow(2.0, 2.0 / 3.0));

        b.push_back(b1); b.push_back(b2); b.push_back(b3); b.push_back(b4); 
        a.push_back(a1); a.push_back(a2); a.push_back(a3); a.push_back(a4); 

    }

        virtual void next(){
            SIA4_method();
            container aux;
            for(auto i: q)
                aux.push_back(i);
            for(auto i: p)
                aux.push_back(i);
            variable = aux;
        }

        virtual ~SIA4(){};
    protected:
        void H(container aux_q, container aux_p, container& Uq, container& Tp);

        void SIA4_method();
        container q, p;
        container a, b;  
};

template <class function>
void SIA4<function>::H(container aux_q, container aux_p, container& Uq, container& Tp){
        container aux;
        for(auto i: aux_q)
            aux.push_back(i);
        for(auto i: aux_p)
            aux.push_back(i);
        function->set(time, aux, __parameter);
        container H(function->get_result());
        for(size_t i = 0; i < H.size() / 2; ++i)
            Uq.push_back(H[i]);
        for(size_t i = H.size() / 2; i < H.size(); ++i)
            Tp.push_back(H[i]);

}

template <class function>
void SIA4<function>::SIA4_method() {
    container aux_q, aux_p;
    aux_q = q;
    aux_p = p;
    time = time + dt;
    for(size_t k = 0; k < 4; ++k){
        container Uq, Tp;
        H(aux_q, aux_p, Uq, Tp);
        for(size_t i = 0; i < p.size(); ++i){
            aux_p[i] = aux_p[i] - b[k] * Uq[i] * dt; 
        }
        Uq.clear(), Tp.clear();
        H(aux_q, aux_p, Uq, Tp);
        for(size_t i = 0; i < q.size(); ++i){
            aux_q[i] = aux_q[i] + a[k] * Tp[i] * dt;
        }
    }
    q = aux_q;
    p = aux_p;
}

#endif /* SIA4_H */
