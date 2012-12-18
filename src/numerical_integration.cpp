#include "numerical_integration.h"

std::ostream & operator<<(std::ostream &out, NumericalIntegration &object) {
    unsigned i;
    for (i = 0; i < object.size_variable()-1; i++)
        out << object[i] << "\t";
    out << object[i];
    return out;
}

NumericalIntegration::NumericalIntegration(type_container variable,type_container parameter,type_data dt):
    __func(),__variable(variable), __parameter(parameter),  __h(dt), __t(0), __model_name(), __method()
{
}

type_data NumericalIntegration::get_t() const{
    return (__t);
}

type_data NumericalIntegration::get_dt() const{
    return (__h);
}

type_data NumericalIntegration::get_variable(unsigned n) const{
    if (n < __variable.size()) {
        return (__variable[n]);
    } else {
        throw Index_error("Wrong access to VARIABLES");
    }
}

type_data NumericalIntegration::get_parameter(unsigned n) const{
    if (n < __parameter.size()) {
        return (__parameter[n]);
    } else {
        throw Index_error("Wrong access to PARAMETERS");
    }
}

unsigned NumericalIntegration::size_variable() const{
    return (__variable.size());
}

unsigned NumericalIntegration::size_parameter() const{
    return (__parameter.size());
}




type_data NumericalIntegration::operator[] (const unsigned nIndex) {
    if (nIndex < __variable.size()) {
        return __variable[nIndex];
    } else {
        throw Index_error("Wrong access to VARIABLES");
    }
}

type_data NumericalIntegration::operator[] (std::string nIndex) {
   if (__func->variable_name_index.count(nIndex)) {
       return __variable[__func->variable_name_index[nIndex]];
    } else {
        throw Index_error("Wrong access to VARIABLES");
    }
}


