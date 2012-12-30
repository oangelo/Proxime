#ifndef LORENZ_H
#define LORENZ_H 

#include "functions.h"


class LorenzFunction : public functions_capsule {
public:
  LorenzFunction();
  void set(type_data &t, type_container & variables, type_container & parameters);
protected:
  type_data dx();
  type_data dy();
  type_data dz();
  
  type_data X,Y,Z;
  type_data sigma,gamma,beta;
};

class Jacobian_LorenzFunction : public functions_capsule {
public:
  Jacobian_LorenzFunction();
  void set(type_data &t, type_container & variables, type_container & parameters);
protected:
  type_data dx();
  type_data dy();
  type_data dz();
  
  type_data X,Y,Z;
  type_data sigma,gamma,beta,X_fiducial,Y_fiducial,Z_fiducial;
};

enum lorenz_enum {
  P_SIGMA_LORENZ,P_GAMMA_LORENZ,P_BETA_LORENZ,P_X,P_Y,P_Z
};


#endif /* LORENZ_H */
