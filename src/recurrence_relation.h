#include <vector>
#include <string>
#include <math.h>

class recurrence_relation{
public:
  recurrence_relation(std::vector<double> variable,std::vector<double> parameter);
  recurrence_relation();
  const double get_t() const;
  const double get_variable(unsigned n) const;
  const std::vector<double> & get_variable() const;
  const double get_parameter(unsigned n) const;
  const std::string & get_model_name() const;
  
  const unsigned size_variable() const;
  const unsigned size_parameter() const;
  
  virtual void next()=0;
  virtual void next(unsigned steps)=0;
  
protected:

  std::vector<double> __variable;
  std::vector<double> __parameter;
  double  __t;
  std::string __model_name; 
}; 

class logistic_map: public recurrence_relation{
public:  
  logistic_map(double x0,double a);
  virtual void next();
  virtual void next(unsigned steps);
  double lyapunov(unsigned n);
};