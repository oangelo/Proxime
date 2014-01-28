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
typedef std::map<std::string, value> labels_values; 
typedef std::map<std::string, size_t> dictionary; 
/****************************************/

class FunctionCapsule {
public:
  FunctionCapsule(std::string function_name, 
                    dictionary index_variables, dictionary index_parameters,
                    labels_values parameters_values);
  virtual ~FunctionCapsule(){};

  virtual void set(value &t, container & variables) = 0;


  const value get_result(unsigned i) const;
  const value operator[](std::string name);
  const container & get_result() const;
  labels_values get_labels_values();
  std::string GetLabel(size_t index);

  //amount of functions and variables 
  unsigned size() const;

  //Dictionry of names and number of variables anda parameters
  dictionary index_var, index_par;
protected:
  std::string function_name;
  container result;
};

#endif /* FUNCTIONS_H*/
