#include <vector>
#include <string>
#include "Exceptions.h"
#include "functions.h"



class Numerical_Integration {
public:
  Numerical_Integration(std::vector<double> variable,std::vector<double> parameter,double dt,std::string model_name);
  const double get_dt() const;
  const double get_t() const;
  const double get_variable(unsigned n) const;
  const std::vector<double> & get_variable() const {return(__variable);};
  const double get_parameter(unsigned n) const;
  const std::string & get_model_name() const {return(__model_name);};
  const std::string & get_method_name() const {return(__method);};
  
  const unsigned size_variable() const;
  const unsigned size_parameter() const;
  
  virtual void next()=0;
  
  //return the variables
  double operator[] (const unsigned nIndex);
  //return the variable ready to print
  friend std::ostream& operator<< (std::ostream &out, Numerical_Integration &object);
  
protected:

  functions_capsule *__func;
  std::vector<double> __variable;
  std::vector<double> __parameter;
  double __h, __t;
  std::string __model_name;
  std::string __method;
  
};

template <class function>
class RungeKutta: public Numerical_Integration{
public:  
  RungeKutta(std::vector<double> variable,std::vector<double> parameter,double dt,std::string model_name)
  :Numerical_Integration(variable,parameter,dt,model_name) {
    /*Point to the class witch encapsulate the functions*/  
    __func=new function;
    (*__func).set(__t, __variable, __parameter);   
  }
  ~RungeKutta();
  virtual void next();
protected:
  void RungeKutta_method();
};

template <class function>
class AdamsBashforth: public Numerical_Integration{
public:
   AdamsBashforth(std::vector<double> variable,std::vector<double> parameter,double dt,std::string model_name="")
  :Numerical_Integration(variable,parameter,dt,model_name) {
    /*Pont to the class witch encapsulate the functions*/  
    __func=new function;
    (*__func).set(__t, __variable, __parameter);
    RungeKutta<function> model(variable,parameter,dt,model_name);
    step1 = model.get_variable();
    model.next();
    step2 = model.get_variable();
    model.next();
    step3 = model.get_variable();
    model.next();
    step4 = model.get_variable();
    __t += 3*__h;
  }
  virtual void next(){
    AdamsBashforth_method();
    step1 = step2;
    step2 = step3;
    step3 = step4;
    step4 = new_step;
    __variable=step1;
    __t = __t + __h;
  }
  
protected:
  void AdamsBashforth_method();
  std::vector<double> step1, step2, step3, step4,new_step;
};

template <class function>
class AdamsMoulton: public AdamsBashforth<function> {
public:
  AdamsMoulton(std::vector<double> variable,std::vector<double> parameter,double dt,std::string model_name="")
  :AdamsBashforth<function>(variable,parameter,dt,model_name) 
  {
    __N=1;
    //Pont to the class witch encapsulate the functions  
    this->__func=new function;
    (this->__func)->set(this->__t, this->__variable, this->__parameter);
    RungeKutta<function> model(variable,parameter,dt,model_name);
    this->step1 = model.get_variable();
    model.next();
    this->step2 = model.get_variable();
    model.next();
    this->step3 = model.get_variable();
    model.next();
    this->step4 = model.get_variable();
    this->__t += 3*this->__h; 
    
  }
  virtual void next(){
    this->AdamsBashforth_method();
    this->step1 = this->step2;
    this->step2 = this->step3;
    this->step3 = this->step4;
    this->step4 = this->new_step;
    for (unsigned i = 0; i < __N; i++) {
        AdamsMoulton_method();
    }  
    this->step4 = this->new_step;
    this->__t  += this->__h;
    this->__variable = this->step1;
  }

  
private:
  void AdamsMoulton_method();
  unsigned __N;
};

#include "Template_Numerical_integration.cpp"