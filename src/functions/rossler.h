#ifndef ROSSLER_H
#define ROSSLER_H 

#include "functions.h"


class RosslerFunction : public functions_capsule {
public:
  RosslerFunction(labels_and_values parameters);
  void set(value &t, container & variables);

 protected:
  value dx();
  value dy();
  value dz();
  
  value X,Y,Z;
  value a,b,c;

};

class Jacobian_RosslerFunction : public functions_capsule {
public:
  Jacobian_RosslerFunction(labels_and_values parameters);
  void set(value &t, container & variables);
protected:
  value dx();
  value dy();
  value dz();
  
  value X,Y,Z;
  value a,b,c,X_fiducial,Y_fiducial,Z_fiducial;
};
#endif /* ROSSLER_H */
