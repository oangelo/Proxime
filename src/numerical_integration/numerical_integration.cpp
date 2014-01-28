#include "numerical_integration.h"

std::ostream & operator<<(std::ostream &out, NumericalIntegration &object) {
    unsigned i;
    out << std::setprecision(12) << object.get_t() << "\t";
    for (i = 0; i < object.size_variable()-1; i++)
        out << object[i] << "\t";
    out << object[i];
    return out;
}

NumericalIntegration::NumericalIntegration(FunctionCapsule& function, container variable, value dt):
function(&function), variable(variable),  dt(dt), time(0), model_name(), method()
{}

value NumericalIntegration::get_t() const{
    return (time);
}

value NumericalIntegration::get_dt() const{
    return (dt);
}

value NumericalIntegration::get_variable(unsigned n) const{
    if (n < variable.size()) {
        return (variable[n]);
    } else {
        throw Index_error("Wrong access to VARIABLES");
    }
}

unsigned NumericalIntegration::size_variable() const{
    return (variable.size());
}


value NumericalIntegration::operator[] (const unsigned nIndex) {
    if (nIndex < variable.size()) {
        return variable[nIndex];
    } else {
        throw Index_error("Wrong access to VARIABLES");
    }
}


labels_and_values NumericalIntegration::GetLabelsValues(){
    labels_and_values aux; 
    for(container::iterator it(variable.begin()); it != variable.end(); ++it){
            std::string label((*function).GetLabel(std::distance(variable.begin(), it)));
            aux[label] = *it;
    }
    return aux;
}

value NumericalIntegration::operator[] (std::string nIndex) {
//   if (function->variable_name_index.count(nIndex)) {
//       return variable[function->variable_name_index[nIndex]];
//    } else {
//        throw Index_error("Wrong access to VARIABLES");
//    }
}

