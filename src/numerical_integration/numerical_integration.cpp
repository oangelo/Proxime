#include "numerical_integration.h"

NumericalIntegration::NumericalIntegration(NumericalIntegration const& other)
:function(other.function->Clone()), variable(other.variable),  dt(other.dt), time(other.time), method(other.method)
{
}

std::ostream & operator<<(std::ostream &out, NumericalIntegration &object) {
    labels_values aux(object.get_labels_values());
    labels_values::iterator it(aux.begin());
    out << std::setprecision(12) << object.get_t() << ",";
    for (; it != --(aux.end()); ++it)
        out << (*it).second << ",";
    ++it;
    out << (*it).second;
    return out;
}

NumericalIntegration::NumericalIntegration(FunctionCapsule& func, labels_values initial_condition, value dt):
function(func.Clone()), variable(initial_condition.size()),  dt(dt), time(0), method()
{
    for(labels_values::iterator it(initial_condition.begin()); it != initial_condition.end(); ++it){
        size_t index = function->index_var[(*it).first];
        variable[index] = (*it).second;
    }
}

value NumericalIntegration::get_t() const{
    return (time);
}

value NumericalIntegration::get_dt() const{
    return (dt);
}

unsigned NumericalIntegration::size() const{
    return (variable.size());
}

const value NumericalIntegration::operator[] (const size_t index) const{
    if (index < variable.size()) {
        return variable[index];
    } else {
        throw Index_error("Wrong access to VARIABLES");
    }
}

const value NumericalIntegration::operator[] (std::string index) const {
    dictionary::iterator it((function->index_var).find(index));
   if (it != (function->index_var).end()) {
       return variable[function->index_var[index]];
    } else {
        throw Index_error("Wrong access to VARIABLES");
    }
}

labels_values NumericalIntegration::get_labels_values() {
    labels_values aux; 
    for(container::iterator it(variable.begin()); it != variable.end(); ++it){
            std::string label((*function).GetLabel(std::distance(variable.begin(), it)));
            aux[label] = *it;
    }
    return aux;
}


