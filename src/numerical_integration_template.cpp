

template <class function>
RungeKutta<function>::~RungeKutta(){}
template <class function>
void RungeKutta<function>::next(){
  RungeKutta_method();
}

template <class function>
void RungeKutta<function>::RungeKutta_method()
{
    //type_data k[__variable.size()][4];
    std::vector<std::vector<double>> k(__variable.size(), std::vector<double>(4, 0));
    type_data t = __t;

    unsigned i;
    type_container aux_argument(__variable.size(), 0);


    //Generating K1*************************************************************
    __func->set(__t, __variable, __parameter);
    for (i = 0; i < __variable.size(); i++) {
        k[i][0] = __func->get_result(i);
        //k[i][0] = ((*__function[i])(__t, __variable, __parameter));
    }

    //Generating K2*************************************************************
    t = __t + 0.5 * __h;
    for (i = 0; i < __variable.size(); i++) {
        aux_argument[i] = __variable[i]+(0.5 * __h * k[i][0]);
    }
    __func->set(t, aux_argument, __parameter);
    for (i = 0; i < __variable.size(); i++) {
        k[i][1] = __func->get_result(i);
        //k[i][1] = ((*__function[i])(t, aux_argument, __parameter));
    }
    //Generating K3*************************************************************
    t = __t + 0.5 * __h;
    for (i = 0; i < __variable.size(); i++) {
        aux_argument[i] = __variable[i]+(0.5 * __h * k[i][1]);
    }
    __func->set(t, aux_argument, __parameter);
    for (i = 0; i < __variable.size(); i++) {
        k[i][2] = __func->get_result(i);
        //k[i][2] = ((*__function[i])(t, aux_argument, __parameter));
    }
    //Generating K4*************************************************************
    t = __t + __h;
    for (i = 0; i < __variable.size(); i++) {
        aux_argument[i] = __variable[i] + __h * k[i][2];
    }
    __func->set(t, aux_argument, __parameter);
    for (i = 0; i < __variable.size(); i++) {
        k[i][3] = __func->get_result(i);
        //k[i][3] = ((*__function[i])(t, aux_argument, __parameter));
    }

    //Generating the final values***********************************************
   
    __t = __t + __h;
    for (i = 0; i < __variable.size(); i++) {
        __variable[i] = __variable[i] + __h * (k[i][0] + 2 * k[i][1] + 2 * k[i][2] + k[i][3]) / 6;
        if (__variable[i] != __variable[i]) {
            throw Value_error("Value error in the Runge Kutta method");
        }
    }
   

} 

template <class function>
void AdamsBashforth<function>::AdamsBashforth_method() {
    type_container aux(__variable.size(), 0);
    std::vector< type_container > results;

    __func->set(__t, step4, __parameter);
    results.push_back(__func->get_result());
    __func->set(__t, step3, __parameter);
    results.push_back(__func->get_result());
    __func->set(__t, step2, __parameter);
    results.push_back(__func->get_result());
    __func->set(__t, step1, __parameter);
    results.push_back(__func->get_result());

    for (unsigned i = 0; i < __variable.size(); i++) {
        aux[i] = 55 * results[0][i];
        aux[i] -= 59 * results[1][i];
        aux[i] += 37 * results[2][i];
        aux[i] -= 9 * results[3][i];
        aux[i] *= (__h / 24);
        aux[i] += step4[i];
        if (aux[i] != aux[i]) {
            throw Value_error("Value error in the Adams Bashforth method");
        }
        new_step=aux;
    }
}

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
