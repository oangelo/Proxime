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

double SimplePendulumEnergy(double theta, double omega, double l, double m, double g);

#endif /* SIMPLE_PENDULUM_H */
