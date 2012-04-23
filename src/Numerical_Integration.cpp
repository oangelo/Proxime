/*
 * File:   Numerical_Integration.cpp
 * Author: angelo
 *
 * Created on October 27, 2009, 5:01 PM
 */

#include "Numerical_Integration.h"

std::ostream & operator<<(std::ostream &out, Numerical_Integration &object) {
    unsigned i;
    for (i = 0; i < object.size_variable()-1; i++)
        out << object[i] << "\t";
    out << object[i];
    return out;
}

Numerical_Integration::Numerical_Integration(std::vector<double> variable,std::vector<double> parameter,double dt,std::string model_name) {
    //storing data//
    __variable=variable;
    __parameter=parameter;
    __h=dt;
    //set initial conditions//
    __t=0;
    __model_name=model_name;
}

const double Numerical_Integration::get_t() const{
    return (__t);
}

const double Numerical_Integration::get_dt() const{
    return (__h);
}

const double Numerical_Integration::get_variable(unsigned n) const{
    if (n < __variable.size()) {
        return (__variable[n]);
    } else {
        throw Index_error("Wrong access to VARIABLES");
    }
}

const double Numerical_Integration::get_parameter(unsigned n) const{
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




double Numerical_Integration::operator[] (const unsigned nIndex) {
    if (nIndex < __variable.size()) {
        return __variable[nIndex];
    } else {
        throw Index_error("Wrong access to VARIABLES");
    }
}
