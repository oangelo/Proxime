#ifndef SIMPLE_PENDULUM_H
#define SIMPLE_PENDULUM_H 

#include "functions.h"

class SimplePendulumFunction : public FunctionCapsule {
public:
    SimplePendulumFunction();
    void set(value &t, container & variables, container & parameters);

    enum variables {
      V_THETA,
      V_OMEGA
    };

    enum parameters {
      P_L, P_G
    };

protected:
    value dTheta();
    value dOmega();

    value theta, omega;
    value l, g;
};

//The equation for the simpletic integration of the simple pendulum
//Derivative of H in ralation to q(theta)
class SimplePendulum_H: public FunctionCapsule {
public:
    SimplePendulum_H();
    void set(value &t, container & variables, container & parameters);

    enum variables {
      V_Q, V_P
    };

    enum parameters {
      P_L, P_G, P_M
    };

protected:
    value Vq();
    value Tp();

    value p, q;
    value l, g, m;
};

value SimplePendulumEnergy(value theta, value omega, value l, value m, value g);

#endif /* SIMPLE_PENDULUM_H */
