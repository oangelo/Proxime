#include "functions.h"

const value functions_capsule::get_result(unsigned i) const{
    return result[i];
}

labels_and_values functions_capsule::GetLabelsValues(){
    labels_and_values aux;
    dictionary::iterator it; 
    for(it = index_var.begin(); it != index_var.end(); ++it){
        std::string name((*it).first);
        aux[name] = result[index_var[name]];
    }
    return aux;
}

const value functions_capsule::operator[](std::string name) {
    return result[index_var[name]];
}

const container & functions_capsule::get_result() const{
    return result;
}

unsigned functions_capsule::size() const{
    return result.size();
}

functions_capsule::functions_capsule(std::string function_name, 
                                     dictionary index_variables, dictionary index_parameters,
                                     labels_and_values parameters_values):
                                     function_name(function_name), result(index_variables.size()),
                                     index_var(index_variables), index_par(index_parameters)
{
//TODO: Check if the labels on parameters are consistnt with index_parameters 
}

std::string functions_capsule::GetLabel(size_t index){
    dictionary::iterator it;
    for(it = index_var.begin(); it != index_var.end(); ++it)
        if((*it).second == index)
            return (*it).first;
}
