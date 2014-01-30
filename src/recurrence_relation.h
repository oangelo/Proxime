#include <vector>
#include <string>
#include <math.h>

#include "functions/functions.h" 

class recurrence_relation{
public:
  recurrence_relation(container variable,container parameter);
  recurrence_relation();
  virtual ~recurrence_relation(){};
  value get_t() const;
  value get_variable(unsigned n) const;
  const container & get_variable() const;
  value get_parameter(unsigned n) const;
  const std::string & get_model_name() const;
  
  unsigned size() const;
  unsigned size_parameter() const;
  
  virtual void next()=0;
  virtual void next(unsigned steps)=0;
  
protected:

  container __variable;
  container __parameter;
  value  __t;
  std::string __model_name; 
}; 

class logistic_map: public recurrence_relation{
public:  
  logistic_map(value x0,value a);
  virtual void next();
  virtual void next(unsigned steps);
  value lyapunov(unsigned n);
};
