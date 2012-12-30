#ifndef ROSSLER_H
#define ROSSLER_H 

#include "functions.h"


class RosslerFunction : public functions_capsule {
public:
  RosslerFunction();
  void set(type_data &t, type_container & variables, type_container & parameters);

 protected:
  type_data dx();
  type_data dy();
  type_data dz();
  
  type_data X,Y,Z;
  type_data a,b,c;

};

class Jacobian_RosslerFunction : public functions_capsule {
public:
  Jacobian_RosslerFunction();
  void set(type_data &t, type_container & variables, type_container & parameters);
protected:
  type_data dx();
  type_data dy();
  type_data dz();
  
  type_data X,Y,Z;
  type_data a,b,c,X_fiducial,Y_fiducial,Z_fiducial;
};
#endif /* ROSSLER_H */
