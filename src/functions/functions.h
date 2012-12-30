#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <vector>
#include <math.h>
#include <pthread.h>
#include <iostream>
#include <map>
#include <string>

typedef double type_data;
typedef std::vector<type_data> type_container;

class functions_capsule {
public:
  functions_capsule();
  virtual ~functions_capsule(){};
  typedef std::pair<std::string, unsigned> name_item;
  typedef std::map<std::string,unsigned> items; 

  virtual void set(type_data &t, type_container & variables, type_container & parameters)=0;
  type_data get_result(unsigned i) const;
  const type_container & get_result() const;
  unsigned size() const;
 
  //map of the variables and index associated to them
  items variable_name_index,parameter_name_index;     
protected:
  std::string func_name;
  type_container __result;
};

#endif /* FUNCTIONS_H*/
