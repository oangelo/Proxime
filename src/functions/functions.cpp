#include "functions.h"

const value functions_capsule::get_result(unsigned i) const{
    return result[i];
}

const value functions_capsule::get_result(std::string name) {
    return result[index_var[name]];
}

const container & functions_capsule::get_result() const{
    return result;
}

unsigned functions_capsule::size() const{
    return result.size();
}

functions_capsule::functions_capsule(std::string function_name, size_t variable_amount):
function_name(function_name), result(variable_amount)
{

}

functions_capsule::functions_capsule(std::string function_name, size_t variable_amount,
                                     dictionary index_variables, dictionary index_parameters):
function_name(function_name), result(variable_amount),
index_var(index_variables), index_par(index_parameters)
{

}
