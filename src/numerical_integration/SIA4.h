#ifndef SIA4_H
#define SIA4_H 

#include "numerical_integration.h"
#include "runge_kutta.h"
#include <iostream>

template <class function>
class SIA4: public NumericalIntegration{
    public:
        SIA4(type_container variable, type_container parameter, type_data dt)
            :NumericalIntegration(variable, parameter, dt), q(), p(), a(), b()
    {
        //Initializing the variables
        for(size_t i = 0; i < variable.size() / 2; ++i){
            q.push_back(variable[i]);
        }
        for(size_t i = variable.size() / 2; i < variable.size(); ++i){
            p.push_back(variable[i]);
        }
        this->__func.reset(new function);
        __func->set(__t, q, __parameter);

        double a1, a2, a3, a4;
        double b1, b2, b3, b4;
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
            std::vector<double> aux;
            for(auto i: q)
                aux.push_back(i);
            for(auto i: p)
                aux.push_back(i);
            __variable = aux;
        }

        virtual ~SIA4(){};
    protected:
        std::vector<double> a;  
        std::vector<double> b;  
        void H(std::vector<double> aux_q, std::vector<double> aux_p, std::vector<double>& Uq, std::vector<double>& Tp);

        void SIA4_method();
        type_container q, p;
};

template <class function>
void SIA4<function>::H(std::vector<double> aux_q, std::vector<double> aux_p, std::vector<double>& Uq, std::vector<double>& Tp){
        std::vector<double> aux;
        for(auto i: aux_q)
            aux.push_back(i);
        for(auto i: aux_p)
            aux.push_back(i);
        __func->set(__t, aux, __parameter);
        std::vector<double> H(__func->get_result());
        for(size_t i = 0; i < H.size() / 2; ++i)
            Uq.push_back(H[i]);
        for(size_t i = H.size() / 2; i < H.size(); ++i)
            Tp.push_back(H[i]);

}

template <class function>
void SIA4<function>::SIA4_method() {
    type_container aux_q, aux_p;
    aux_q = q;
    aux_p = p;
    __t = __t + __h;
    for(size_t k = 0; k < 4; ++k){
        for(size_t i = 0; i < p.size(); ++i){
            std::vector<double> Uq, Tp;
            H(aux_q, aux_p, Uq, Tp);
            aux_p[i] = aux_p[i] - b[k] * Uq[i] * __h; 

        }
        for(size_t i = 0; i < q.size(); ++i){
            std::vector<double> Uq, Tp;
            H(aux_q, aux_p, Uq, Tp);
            aux_q[i] = aux_q[i] + a[k] * Tp[i] * __h;
        }
    }
    q = aux_q;
    p = aux_p;
}

#endif /* SIA4_H */
