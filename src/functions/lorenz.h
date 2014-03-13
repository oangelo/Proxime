#ifndef LORENZ_H
#define LORENZ_H 

#include "functions.h"


class Lorenz: public FunctionCapsule {
public:
  Lorenz(labels_values parameters);
  void set(value &t, container & variables);
  virtual Lorenz* Clone() const;
  virtual Lorenz* Create(labels_values parameters) const;
protected:
  value dx();
  value dy();
  value dz();
  
  value X,Y,Z;
  value sigma,gamma,beta;
};

class Jacobian_Lorenz: public FunctionCapsule {
public:
  Jacobian_Lorenz(labels_values parameters);
  void set(value &t, container & variables);
  virtual Jacobian_Lorenz* Clone() const;
  virtual Jacobian_Lorenz* Create(labels_values parameters) const;
protected:
  value dx();
  value dy();
  value dz();
  
  value X,Y,Z;
  value sigma,gamma,beta,X_fiducial,Y_fiducial,Z_fiducial;
};

#endif /* LORENZ_H */
