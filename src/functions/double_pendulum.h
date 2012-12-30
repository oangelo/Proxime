#ifndef DOUBLE_PENDULUM_H
#define DOUBLE_PENDULUM_H 

#include "functions.h"

class DoublePendulumFunction : public functions_capsule {
public:
  DoublePendulumFunction();
  void set(type_data &t, type_container & variables, type_container & parameters);
protected:
  type_data dTheta1();
  type_data dTheta2();
  type_data dOmega1();
  type_data dOmega2();
  
  type_data theta1,theta2,omega1,omega2;
  type_data l1,l2,m1,m2,g;

};

class Jacobian_DoublePendulumFunction : public functions_capsule {
public:
  Jacobian_DoublePendulumFunction();
  void set(type_data &t, type_container & variables, type_container & parameters);
protected:
  void Matrix_Jacob(type_data theta1, type_data theta2,type_data omega1, type_data omega2);
  type_data JdTheta1();
  type_data JdTheta2();
  type_data JdOmega1();
  type_data JdOmega2();

  type_data theta1,theta2,omega1,omega2;
  type_data l1,l2,m1,m2,g;
  type_data Jacobian[4][4];
};

enum variables {
  V_THETA1, V_THETA2, V_OMEGA1, V_OMEGA2
};

enum parameters {
  P_L1, P_L2, P_M1, P_M2, P_G, P_OMEGA0, P_A, P_SIGMA, P_THETA1, P_THETA2, P_OMEGA1, P_OMEGA2
};

#endif /* DOUBLE_PENDULUM_H */
