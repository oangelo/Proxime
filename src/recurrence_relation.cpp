#include "recurrence_relation.h"

recurrence_relation::recurrence_relation(std::vector<double> variable,std::vector<double> parameter){
  __variable=variable;
  __parameter=parameter;
  __t=0;
}
recurrence_relation::recurrence_relation(){__t=0;}
const double recurrence_relation::get_t() const{return(__t);};
const double recurrence_relation::get_variable(unsigned n) const{return(__variable[n]);};
const std::vector<double> & recurrence_relation::get_variable() const {return(__variable);};
const double recurrence_relation::get_parameter(unsigned n) const{return(__parameter[n]);};
const std::string & recurrence_relation::get_model_name() const {return(__model_name);};
const unsigned recurrence_relation::size_variable() const{return(__variable.size());};
const unsigned recurrence_relation::size_parameter() const{return(__parameter.size());};

logistic_map::logistic_map(double x0,double a){
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

double logistic_map::lyapunov(unsigned n){
  double sum=0;
  for(unsigned counter=0;counter<n;counter++){
    sum+=log(fabs(__parameter[0]*(1-2*__variable[0])));
    this->next();
  }
  return(sum/n);
}