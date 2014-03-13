#ifndef SIMPLE_PENDULUM_T_H
#define SIMPLE_PENDULUM_T_H 
#include <fstream>
#include "../src/functions/simple_pendulum.h"

TEST(Simple_Pendulum, Static) {
    container variable(2),parameter(2);
    variable[SimplePendulum::V_THETA] = 0.0;
    variable[SimplePendulum::V_OMEGA] = 0.0;
    parameter[SimplePendulum::P_L]= 0.30;
    parameter[SimplePendulum::P_G]= 9.8;
    

    RungeKutta4Th<SimplePendulum> model(variable, parameter, 0.0001); 
    for(size_t i = 0; i < 10000; ++i){
        model.next();
        EXPECT_NEAR(model[0], 0.0, 0.0001);
    }
}

TEST(Simple_Pendulum, sin) {
    value theta_0 = M_PI / 40;
    value l=1; value g = 1;
    value delta_t = 0.0001;

    container variable(2),parameter(2);
    variable[SimplePendulum::V_THETA] = theta_0;
    variable[SimplePendulum::V_OMEGA] = 0.0;
    parameter[SimplePendulum::P_L]= l;
    parameter[SimplePendulum::P_G]= g;

    AdamsBashforth4Th<SimplePendulum> model(variable, parameter, delta_t); 
    for(size_t i = 0; i < 100000; ++i){
        value t = delta_t*i; 
        model.next();
        value test = theta_0*cos(sqrt(g/l)*(t));
        EXPECT_NEAR(model[0], test, 0.0025);
    }
}

TEST(Simple_Pendulum, CrossModels) {
    value theta_0 = M_PI / 2;
    value l=0.30, g = 9.8, m = 0.1;
    value delta_t = 0.0001;

    container variable(2),parameter(3);
    variable[SimplePendulum::V_THETA] = theta_0;
    variable[SimplePendulum::V_OMEGA] = 0.0;
    parameter[SimplePendulum::P_L]= l;
    parameter[SimplePendulum::P_G]= g;
    parameter[SimplePendulum_H::P_M]= m;

    AdamsBashforth4Th<SimplePendulum> model1(variable, parameter, delta_t); 
    RungeKutta4Th<SimplePendulum> model2(variable, parameter, delta_t); 
    AdamsMoulton4Th<SimplePendulum> model3(variable, parameter, delta_t); 
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
    value theta_0 = M_PI / 2;
    value l=0.30; value g = 9.8;
    value delta_t = 0.0001;

    container variable(2),parameter(2);
    variable[SimplePendulum::V_THETA] = theta_0;
    variable[SimplePendulum::V_OMEGA] = 0.0;
    parameter[SimplePendulum::P_L]= l;
    parameter[SimplePendulum::P_G]= g;

    RungeKutta4Th<SimplePendulum> model(variable, parameter, delta_t); 
    value energy_0 = SimplePendulumEnergy(theta_0, 0, l, 1, g);
    for(size_t i = 0; i < 100000; ++i){
        model.next();
        EXPECT_NEAR(SimplePendulumEnergy(model[0], model[1], l, 1, g), energy_0, pow(10, -10));
    }
}

TEST(Simple_Pendulum, EnergySimpletic) {
    value theta_0 = M_PI / 2;
    value l=0.30, g = 9.8, m = 0.1;
    value delta_t = 0.000001;

    container variable(2),parameter(3);
    variable[SimplePendulum_H::V_Q] = theta_0;
    variable[SimplePendulum_H::V_P] = 0.0;
    parameter[SimplePendulum_H::P_L]= l;
    parameter[SimplePendulum_H::P_G]= g;
    parameter[SimplePendulum_H::P_M]= m;

    std::ofstream myfile;
    myfile.open("data/simple_pendulum_energy_SIA4.out");
    SIA4<SimplePendulum_H> model(variable, parameter, delta_t); 
    value energy_0 = SimplePendulumEnergy(theta_0, 0, l, m, g);
    for(size_t i = 0; i < 10000000; ++i){
        model.next();
        value energy = SimplePendulumEnergy(model[0], model[1] / (m * pow(l, 2)), l, m, g);
        myfile << delta_t * i << " " << energy << std::endl;
        EXPECT_NEAR(energy, energy_0, pow(10,-10));
    }
}

TEST(Simple_Pendulum, Phase_Space) {
    value l=1; value g = 1; value m = 1;
    value delta_t = 0.001;

    container variable(2),parameter(3);
    variable[SimplePendulum::V_THETA] = 0.0;
    variable[SimplePendulum::V_OMEGA] = 0.0;
    parameter[SimplePendulum::P_L]= l;
    parameter[SimplePendulum::P_G]= g;
    parameter[SimplePendulum_H::P_M]= m;

    std::ofstream myfile1, myfile2;
    myfile1.open("data/Simple_pendulum_phase_space_SIA4.out");
    myfile2.open("data/Simple_pendulum_phase_space_RK.out");

    for(value angle = 0.001; angle < 1.5 ; angle += 0.1){
        variable[SimplePendulum::V_OMEGA] = angle;
        SIA4<SimplePendulum_H> model1(variable, parameter, delta_t); 
        AdamsBashforth4Th<SimplePendulum> model2(variable, parameter, delta_t); 
        for(size_t i = 0; i < 7500; ++i){
            model1.next();
            model2.next();
            myfile1 << model1 << std::endl; 
            myfile2 << model2 << std::endl; 
        }
    }
}



#endif /* SIMPLE_PENDULUM_T_H */
