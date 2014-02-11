#ifndef LORENZ_H
#define LORENZ_H 

#include "functions.h"


class LorenzFunction : public FunctionCapsule {
public:
  LorenzFunction(labels_values parameters);
  void set(value &t, container & variables);
  virtual LorenzFunction* Clone() const;
  virtual LorenzFunction* Create(labels_values parameters) const;
protected:
  value dx();
  value dy();
  value dz();
  
  value X,Y,Z;
  value sigma,gamma,beta;
};

class Jacobian_LorenzFunction : public FunctionCapsule {
public:
  Jacobian_LorenzFunction(labels_values parameters);
  void set(value &t, container & variables);
  virtual Jacobian_LorenzFunction* Clone() const;
  virtual Jacobian_LorenzFunction* Create(labels_values parameters) const;
protected:
  value dx();
  value dy();
  value dz();
  
  value X,Y,Z;
  value sigma,gamma,beta,X_fiducial,Y_fiducial,Z_fiducial;
};

#endif /* LORENZ_H */
