#ifndef SIMPLE_PENDULUM_T_H
#define SIMPLE_PENDULUM_T_H 
#include <fstream>
#include "../src/functions/simple_pendulum.h"

TEST(Simple_Pendulum, Static) {
    std::vector<double> variable(2),parameter(2);
    variable[SimplePendulumFunction::V_THETA] = 0.0;
    variable[SimplePendulumFunction::V_OMEGA] = 0.0;
    parameter[SimplePendulumFunction::P_L]= 0.30;
    parameter[SimplePendulumFunction::P_G]= 9.8;
    

    RungeKutta<SimplePendulumFunction> model(variable, parameter, 0.0001); 
    for(size_t i = 0; i < 10000; ++i){
        model.next();
        EXPECT_NEAR(model[0], 0.0, 0.0001);
    }
}

TEST(Simple_Pendulum, sin) {
    double theta_0 = M_PI / 40;
    double l=0.30; double g = 9.8;
    double delta_t = 0.0001;

    std::vector<double> variable(2),parameter(2);
    variable[SimplePendulumFunction::V_THETA] = theta_0;
    variable[SimplePendulumFunction::V_OMEGA] = 0.0;
    parameter[SimplePendulumFunction::P_L]= l;
    parameter[SimplePendulumFunction::P_G]= g;

    AdamsBashforth<SimplePendulumFunction> model(variable, parameter, delta_t); 
    //RungeKutta<SimplePendulumFunction> model(variable, parameter, delta_t); 
//    std::ofstream myfile;
//    myfile.open("sp.out");
    for(size_t i = 0; i < 100000; ++i){
        double t = delta_t*i; 
        model.next();
        double test = theta_0*cos(sqrt(g/l)*(t));
        EXPECT_NEAR(model[0], test, 0.0025);
//        myfile << t << " " << model[0] << " " << model[1] << " " << theta_0*cos(sqrt(g/l)*(t)) << std::endl; 
    }
}

TEST(Simple_Pendulum, CrossModels) {
    double theta_0 = M_PI / 40;
    double l=0.30; double g = 9.8;
    double delta_t = 0.0001;

    std::vector<double> variable(2),parameter(2);
    variable[SimplePendulumFunction::V_THETA] = theta_0;
    variable[SimplePendulumFunction::V_OMEGA] = 0.0;
    parameter[SimplePendulumFunction::P_L]= l;
    parameter[SimplePendulumFunction::P_G]= g;

    AdamsBashforth<SimplePendulumFunction> model1(variable, parameter, delta_t); 
    RungeKutta<SimplePendulumFunction> model2(variable, parameter, delta_t); 
    AdamsMoulton<SimplePendulumFunction> model3(variable, parameter, delta_t); 

    for(size_t i = 0; i < 100000; ++i){
        model1.next();
        model2.next();
        model3.next();
        EXPECT_NEAR(model1[0], model2[0], 0.001);
        EXPECT_NEAR(model1[1], model2[1], 0.001);
        EXPECT_NEAR(model1[0], model3[0], 0.001);
        EXPECT_NEAR(model1[1], model3[1], 0.001);

    }
}

TEST(Simple_Pendulum, Phase_Space) {
    double theta_0 = M_PI / 40;
    double l=0.30; double g = 9.8;
    double delta_t = 0.001;

    std::vector<double> variable(2),parameter(2);
    variable[SimplePendulumFunction::V_THETA] = theta_0;
    variable[SimplePendulumFunction::V_OMEGA] = 0.0;
    parameter[SimplePendulumFunction::P_L]= l;
    parameter[SimplePendulumFunction::P_G]= g;

    std::ofstream myfile;
    myfile.open("Simple_pendulum_phase_space.out");


    RungeKutta<SimplePendulumFunction> model(variable, parameter, delta_t); 
    for(double angle = 0.001; angle < 10 ; angle += 0.5){
        variable[SimplePendulumFunction::V_OMEGA] = angle;
        RungeKutta<SimplePendulumFunction> model(variable, parameter, delta_t); 
        for(size_t i = 0; i < 2000; ++i){
            model.next();
            myfile << model << std::endl; 
        }
    }
}

#endif /* SIMPLE_PENDULUM_T_H */
