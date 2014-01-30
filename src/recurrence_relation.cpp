#include "recurrence_relation.h"

recurrence_relation::recurrence_relation(container variable,container parameter):
__variable(variable),
__parameter(parameter),
__t(0),
__model_name()
{}

recurrence_relation::recurrence_relation():
__variable(),
__parameter(),
__t(0),
__model_name()
{}
value recurrence_relation::get_t() const{ return(__t);}
value recurrence_relation::get_variable(unsigned n) const{return(__variable[n]);}
const container & recurrence_relation::get_variable() const {return(__variable);}
value recurrence_relation::get_parameter(unsigned n) const{return(__parameter[n]);}
const std::string & recurrence_relation::get_model_name() const {return(__model_name);}
unsigned recurrence_relation::size() const{return(__variable.size());}
unsigned recurrence_relation::size_parameter() const{return(__parameter.size());}

logistic_map::logistic_map(value x0,value a){
  __model_name="Logistic Equation";
  __parameter.push_back(a);
  __variable.push_back(x0);
}
void logistic_map::next(){
  __variable[0]=__parameter[0]*__variable[0]*(1-__variable[0]);
  __t++;
}
void logistic_map::next(unsigned steps){
  for(unsigned counter=0;counter<steps;counter++)
    __variable[0]=__parameter[0]*__variable[0]*(1-__variable[0]);
  __t+=steps;
}

value logistic_map::lyapunov(unsigned n){
  value sum=0;
  for(unsigned counter=0;counter<n;counter++){
    sum+=log(fabs(__parameter[0]*(1-2*__variable[0])));
    this->next();
  }
  return(sum/n);
}
