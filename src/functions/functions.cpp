#include "functions.h"

const value FunctionCapsule::get_result(unsigned i) const{
    return result[i];
}

labels_values FunctionCapsule::get_labels_values(){
    labels_values aux;
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
                                     labels_values parameters_values):
                                     index_var(index_variables), index_par(index_parameters),
                                     function_name(function_name), result(index_variables.size())
{
//TODO: Check if the labels on parameters are consistnt with index_parameters 
}

std::string FunctionCapsule::GetLabel(size_t index){
    dictionary::iterator it;
    std::string var_label("");
    for(it = index_var.begin(); it != index_var.end(); ++it){
        if((*it).second == index)
            var_label = (*it).first;
    }
    if(var_label == "")
        std::cerr << "Wrong index access to function results" << std::endl;
    //TODO: trow exception 
    return var_label;
}
