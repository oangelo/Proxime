#ifndef ROSSLER_H
#define ROSSLER_H 

#include "functions.h"


class RosslerFunction : public FunctionCapsule{
public:
  RosslerFunction(labels_values parameters);
  void set(value &t, container & variables);
  virtual RosslerFunction* Clone() const;
  virtual RosslerFunction* Create(labels_values parameters) const;

 protected:
  value dx();
  value dy();
  value dz();
  
  value X,Y,Z;
  value a,b,c;
};

class Jacobian_RosslerFunction : public FunctionCapsule {
public:
  Jacobian_RosslerFunction(labels_values parameters);
  void set(value &t, container & variables);
  virtual Jacobian_RosslerFunction* Clone() const;
  virtual Jacobian_RosslerFunction* Create(labels_values parameters) const;
protected:
  value dx();
  value dy();
  value dz();
  
  value X,Y,Z;
  value a,b,c,X_fiducial,Y_fiducial,Z_fiducial;
};
#endif /* ROSSLER_H */
