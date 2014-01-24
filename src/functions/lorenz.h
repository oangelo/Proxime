#ifndef LORENZ_H
#define LORENZ_H 

#include "functions.h"


class LorenzFunction : public functions_capsule {
public:
  LorenzFunction();
  void set(value &t, container & variables, container & parameters);
protected:
  value dx();
  value dy();
  value dz();
  
  value X,Y,Z;
  value sigma,gamma,beta;
};

class Jacobian_LorenzFunction : public functions_capsule {
public:
  Jacobian_LorenzFunction();
  void set(value &t, container & variables, container & parameters);
protected:
  value dx();
  value dy();
  value dz();
  
  value X,Y,Z;
  value sigma,gamma,beta,X_fiducial,Y_fiducial,Z_fiducial;
};

enum lorenz_enum {
  P_SIGMA_LORENZ,P_GAMMA_LORENZ,P_BETA_LORENZ,P_X,P_Y,P_Z
};


#endif /* LORENZ_H */
