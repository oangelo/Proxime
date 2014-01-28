#include "functions.h"

const value FunctionCapsule::get_result(unsigned i) const{
    return result[i];
}

labels_and_values FunctionCapsule::GetLabelsValues(){
    labels_and_values aux;
    dictionary::iterator it; 
    for(it = index_var.begin(); it != index_var.end(); ++it){
        std::string name((*it).first);
        aux[name] = result[index_var[name]];
    }
    return aux;
}

const value FunctionCapsule::operator[](std::string name) {
    return result[index_var[name]];
}

const container & FunctionCapsule::get_result() const{
    return result;
}

unsigned FunctionCapsule::size() const{
    return result.size();
}

FunctionCapsule::FunctionCapsule(std::string function_name, 
                                     dictionary index_variables, dictionary index_parameters,
                                     labels_and_values parameters_values):
                                     function_name(function_name), result(index_variables.size()),
                                     index_var(index_variables), index_par(index_parameters)
{
//TODO: Check if the labels on parameters are consistnt with index_parameters 
}

std::string FunctionCapsule::GetLabel(size_t index){
    dictionary::iterator it;
    for(it = index_var.begin(); it != index_var.end(); ++it)
        if((*it).second == index)
            return (*it).first;
}
