#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <vector>
#include <math.h>
#include <pthread.h>
#include <iostream>
#include <map>
#include <string>


/****************************************
 * These types must be used every where,*
 * to make the precision consistent.    *
 ****************************************/
typedef double value;
typedef std::vector<value> container ;
typedef std::map<std::string, double> labels; 
typedef std::map<std::string, size_t> dictionary; 
/****************************************/

class functions_capsule {
public:
  functions_capsule(std::string function_name, size_t variable_amount, dictionary index_variables, dictionary index_parameters);
  virtual ~functions_capsule(){};

  virtual void set(value &t, container & variables, container & parameters) = 0;


  const value get_result(unsigned i) const;
  const value get_result(std::string name);
  const container & get_result() const;

  //amount of functions and variables 
  unsigned size() const;

  //Dictionry of names and number of variables anda parameters
  dictionary index_var, index_par;
protected:
  std::string function_name;
  container result;
};

#endif /* FUNCTIONS_H*/
