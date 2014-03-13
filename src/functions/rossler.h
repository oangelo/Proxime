#ifndef ROSSLER_H
#define ROSSLER_H 

#include "functions.h"


class Rossler: public FunctionCapsule{
public:
  Rossler(labels_values parameters);
  void set(value &t, container & variables);
  virtual Rossler* Clone() const;
  virtual Rossler* Create(labels_values parameters) const;

 protected:
  value dx();
  value dy();
  value dz();
  
  value X,Y,Z;
  value a,b,c;
};

class Jacobian_Rossler: public FunctionCapsule {
public:
  Jacobian_Rossler(labels_values parameters);
  void set(value &t, container & variables);
  virtual Jacobian_Rossler* Clone() const;
  virtual Jacobian_Rossler* Create(labels_values parameters) const;
protected:
  value dx();
  value dy();
  value dz();
  
  value X,Y,Z;
  value a,b,c,X_fiducial,Y_fiducial,Z_fiducial;
};
#endif /* ROSSLER_H */
