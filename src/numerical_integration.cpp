#include "numerical_integration.h"

std::ostream & operator<<(std::ostream &out, Numerical_Integration &object) {
    unsigned i;
    for (i = 0; i < object.size_variable()-1; i++)
        out << object[i] << "\t";
    out << object[i];
    return out;
}

Numerical_Integration::Numerical_Integration(type_container variable,type_container parameter,type_data dt):
    __variable(variable), __parameter(parameter),  __h(dt), __t(0)
{
}

const type_data Numerical_Integration::get_t() const{
    return (__t);
}

const type_data Numerical_Integration::get_dt() const{
    return (__h);
}

const type_data Numerical_Integration::get_variable(unsigned n) const{
    if (n < __variable.size()) {
        return (__variable[n]);
    } else {
        throw Index_error("Wrong access to VARIABLES");
    }
}

const type_data Numerical_Integration::get_parameter(unsigned n) const{
    if (n < __parameter.size()) {
        return (__parameter[n]);
    } else {
        throw Index_error("Wrong access to PARAMETERS");
    }
}

const unsigned Numerical_Integration::size_variable() const{
    return (__variable.size());
}

const unsigned Numerical_Integration::size_parameter() const{
    return (__parameter.size());
}




type_data Numerical_Integration::operator[] (const unsigned nIndex) {
    if (nIndex < __variable.size()) {
        return __variable[nIndex];
    } else {
        throw Index_error("Wrong access to VARIABLES");
    }
}

type_data Numerical_Integration::operator[] (std::string nIndex) {
   if (__func->variable_name_index.count(nIndex)) {
       return __variable[__func->variable_name_index[nIndex]];
    } else {
        throw Index_error("Wrong access to VARIABLES");
    }
}


