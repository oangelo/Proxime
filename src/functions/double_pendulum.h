#ifndef DOUBLE_PENDULUM_H
#define DOUBLE_PENDULUM_H 

#include "functions.h"

class DoublePendulum : public FunctionCapsule {
public:
  DoublePendulum(labels_values parameters);
  virtual void set(value &t, container & variables);
  virtual DoublePendulum* Clone() const; 
  virtual DoublePendulum* Create(labels_values parameters) const; 

protected:
  value dTheta1();
  value dTheta2();
  value dOmega1();
  value dOmega2();
  
  value theta1,theta2,omega1,omega2;
  value l1,l2,m1,m2,g;

};

class Jacobian_DoublePendulum: public FunctionCapsule {
public:
  Jacobian_DoublePendulum(labels_values parameters);
  void set(value &t, container & variables);
  virtual Jacobian_DoublePendulum* Clone() const; 
  virtual Jacobian_DoublePendulum* Create(labels_values parameters) const; 
protected:
  void Matrix_Jacob(value theta1, value theta2, value omega1, value omega2);
  value JdTheta1();
  value JdTheta2();
  value JdOmega1();
  value JdOmega2();

  value theta1, theta2, omega1, omega2;
  value l1, l2, m1, m2, g;
  value Jacobian[4][4];
};

enum variables {
  V_THETA1, V_THETA2, V_OMEGA1, V_OMEGA2
};

enum parameters {
  P_L1, P_L2, P_M1, P_M2, P_G, P_THETA1, P_THETA2, P_OMEGA1, P_OMEGA2
};

value DoublePendulumHamiltonian(value q1, value q2, value p1, value p2, value l1, value l2, value m1, value m2, value g);
value DoublePendulumEnergy(value theta1, value theta2, value omega1, value omega2, value l1, value l2, value m1, value m2, value g);
#endif /* DOUBLE_PENDULUM_H */
