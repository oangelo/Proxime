#ifndef SIMPLE_PENDULUM_T_H
#define SIMPLE_PENDULUM_T_H 

#include "../src/functions/simple_pendulum.h"

TEST(Simple_Pendulum, Static) {
    std::vector<double> variable(2),parameter(2);
    variable[SimplePendulumFunction::V_THETA] = 0.0;
    variable[SimplePendulumFunction::V_OMEGA] = 0.0;
    parameter[SimplePendulumFunction::P_L]= 0.30;
    parameter[SimplePendulumFunction::P_G]= 9.8;

    RungeKutta<SimplePendulumFunction> model(variable, parameter, 0.0001); 
    for(size_t i=0:10){
        model.next();
        EXPECT_NEAR(model[0], 0.0, 0.0001);
    }
}

#endif /* SIMPLE_PENDULUM_T_H */
