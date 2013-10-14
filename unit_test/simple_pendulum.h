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
    double l=1; double g = 1;
    double delta_t = 0.0001;

    std::vector<double> variable(2),parameter(2);
    variable[SimplePendulumFunction::V_THETA] = theta_0;
    variable[SimplePendulumFunction::V_OMEGA] = 0.0;
    parameter[SimplePendulumFunction::P_L]= l;
    parameter[SimplePendulumFunction::P_G]= g;

    AdamsBashforth<SimplePendulumFunction> model(variable, parameter, delta_t); 
    for(size_t i = 0; i < 100000; ++i){
        double t = delta_t*i; 
        model.next();
        double test = theta_0*cos(sqrt(g/l)*(t));
        EXPECT_NEAR(model[0], test, 0.0025);
    }
}

TEST(Simple_Pendulum, CrossModels) {
    double theta_0 = M_PI / 2;
    double l=0.30, g = 9.8, m = 0.1;
    double delta_t = 0.0001;

    std::vector<double> variable(2),parameter(3);
    variable[SimplePendulumFunction::V_THETA] = theta_0;
    variable[SimplePendulumFunction::V_OMEGA] = 0.0;
    parameter[SimplePendulumFunction::P_L]= l;
    parameter[SimplePendulumFunction::P_G]= g;
    parameter[SimplePendulum_H::P_M]= m;

    AdamsBashforth<SimplePendulumFunction> model1(variable, parameter, delta_t); 
    RungeKutta<SimplePendulumFunction> model2(variable, parameter, delta_t); 
    AdamsMoulton<SimplePendulumFunction> model3(variable, parameter, delta_t); 
    SIA4<SimplePendulum_H> model4(variable, parameter, delta_t); 

    for(size_t i = 0; i < 100000; ++i){
        model1.next();
        model2.next();
        model3.next();
        model4.next();
        EXPECT_NEAR(model1[0], model2[0], 0.00001);
        EXPECT_NEAR(model1[1], model2[1], 0.00001);
        EXPECT_NEAR(model1[0], model3[0], 0.00001);
        EXPECT_NEAR(model1[1], model3[1], 0.00001);
        EXPECT_NEAR(model1[0], model4[0], 0.00001);
        EXPECT_NEAR(model1[1], model4[1] / (m * pow(l, 2)), 0.00001);

    }
}


TEST(Simple_Pendulum, Energy) {
    double theta_0 = M_PI / 2;
    double l=0.30; double g = 9.8;
    double delta_t = 0.0001;

    std::vector<double> variable(2),parameter(2);
    variable[SimplePendulumFunction::V_THETA] = theta_0;
    variable[SimplePendulumFunction::V_OMEGA] = 0.0;
    parameter[SimplePendulumFunction::P_L]= l;
    parameter[SimplePendulumFunction::P_G]= g;

    RungeKutta<SimplePendulumFunction> model(variable, parameter, delta_t); 
    double energy_0 = SimplePendulumEnergy(theta_0, 0, l, 1, g);
    for(size_t i = 0; i < 100000; ++i){
        model.next();
        EXPECT_NEAR(SimplePendulumEnergy(model[0], model[1], l, 1, g), energy_0, pow(10, -10));
    }
}

TEST(Simple_Pendulum, Energy_simpletic) {
    double theta_0 = M_PI / 2;
    double l=0.30, g = 9.8, m = 0.1;
    double delta_t = 0.0001;

    std::vector<double> variable(2),parameter(3);
    variable[SimplePendulum_H::V_Q] = theta_0;
    variable[SimplePendulum_H::V_P] = 0.0;
    parameter[SimplePendulum_H::P_L]= l;
    parameter[SimplePendulum_H::P_G]= g;
    parameter[SimplePendulum_H::P_M]= m;

    SIA4<SimplePendulum_H> model(variable, parameter, delta_t); 
    double energy_0 = SimplePendulumEnergy(theta_0, 0, l, m, g);
    for(size_t i = 0; i < 100000; ++i){
        model.next();
        EXPECT_NEAR(SimplePendulumEnergy(model[0], model[1] / (m * pow(l, 2)), l, m, g), energy_0, pow(10,-10));
    }
}

TEST(Simple_Pendulum, Phase_Space) {
    double l=1; double g = 1; double m = 1;
    double delta_t = 0.001;

    std::vector<double> variable(2),parameter(3);
    variable[SimplePendulumFunction::V_THETA] = 0.0;
    variable[SimplePendulumFunction::V_OMEGA] = 0.0;
    parameter[SimplePendulumFunction::P_L]= l;
    parameter[SimplePendulumFunction::P_G]= g;
    parameter[SimplePendulum_H::P_M]= m;

    std::ofstream myfile1, myfile2;
    myfile1.open("data/Simple_pendulum_phase_space_SIA4.out");
    myfile2.open("data/Simple_pendulum_phase_space_RK.out");

    for(double angle = 0.001; angle < 1.5 ; angle += 0.1){
        variable[SimplePendulumFunction::V_OMEGA] = angle;
        SIA4<SimplePendulum_H> model1(variable, parameter, delta_t); 
        AdamsBashforth<SimplePendulumFunction> model2(variable, parameter, delta_t); 
        for(size_t i = 0; i < 7500; ++i){
            model1.next();
            model2.next();
            myfile1 << model1 << std::endl; 
            myfile2 << model2 << std::endl; 
        }
    }
}



#endif /* SIMPLE_PENDULUM_T_H */
