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
  const container & get_result() const;
  unsigned size() const;
 
  //map of the variables and index associated to them
  items variable_name_index,parameter_name_index;     
protected:
  std::string func_name;
  container __result;
};

#endif /* FUNCTIONS_H*/
