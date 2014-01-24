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
/****************************************/

class functions_capsule {
public:
  functions_capsule();
  virtual ~functions_capsule(){};

  virtual void set(value &t, container & variables, container & parameters) = 0;

  value get_result(unsigned i) const;
  value get_result(std::string) const;
  const container & get_result() const;

  //amount of functions and variables 
  unsigned size() const;

  typedef std::map<std::string, size_t> items; 
  items variable_name_index,parameter_name_index;
protected:
  std::string func_name;
  container result;
};

#endif /* FUNCTIONS_H*/
