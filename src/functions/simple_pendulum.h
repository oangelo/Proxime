#ifndef SIMPLE_PENDULUM_H
#define SIMPLE_PENDULUM_H 

#include "functions.h"

class SimplePendulumFunction : public functions_capsule {
public:
    SimplePendulumFunction();
    void set(type_data &t, type_container & variables, type_container & parameters);

    enum variables {
      V_THETA,
      V_OMEGA
    };

    enum parameters {
      P_L, P_G
    };

protected:
    type_data dTheta();
    type_data dOmega();

    type_data theta, omega;
    type_data l, g;
};

//The equation for the simpletic integration of the simple pendulum
//Derivative of H in ralation to q(theta)
class SimplePendulum_H: public functions_capsule {
public:
    SimplePendulum_H();
    void set(type_data &t, type_container & variables, type_container & parameters);

    enum variables {
      V_Q, V_P
    };

    enum parameters {
      P_L, P_G, P_M
    };

protected:
    type_data Vq();
    type_data Tp();

    type_data p, q;
    type_data l, g, m;
};

double SimplePendulumEnergy(double theta, double omega, double l, double m, double g);

#endif /* SIMPLE_PENDULUM_H */
