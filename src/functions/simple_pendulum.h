#ifndef SIMPLE_PENDULUM_H
#define SIMPLE_PENDULUM_H 

#include "functions.h"

class SimplePendulum : public FunctionCapsule {
public:
    SimplePendulum(labels_values parameters);
    void set(value &t, container & variables);
    SimplePendulum* Clone() const;
    SimplePendulum* Create(labels_values parameters) const;

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



value SimplePendulumEnergy(value theta, value omega, value l, value m, value g);

#endif /* SIMPLE_PENDULUM_H */
