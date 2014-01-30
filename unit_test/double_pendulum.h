#ifndef value_PENDULUM_H
#define value_PENDULUM_H 

#include "../src/functions/double_pendulum.h"
#include <fstream>
#include <iostream>

TEST(DoublePendulum, ZeroLyapunov) {
    container variable(4);
    variable[V_THETA1] = M_PI / 10;
    variable[V_THETA2] = M_PI / 10;
    variable[V_OMEGA1] = 0.0;
    variable[V_OMEGA2] = 0.0;

    labels_values parameters;
    parameters["l1"]= 0.30;
    parameters["l2"]= 0.30;
    parameters["m1"]= 0.1;
    parameters["m2"]= 0.1;
    parameters["g"]= 9.8;

    DoublePendulumFunction f(parameters);


    RungeKutta model(f, variable, 0.0001); 
    value max_lyapunov = MaxLyapunov<Jacobian_DoublePendulumFunction>(model, parameters, 2000000, 20000);
    EXPECT_NEAR(max_lyapunov, 0, 0.01);
}

TEST(DoublePendulum, PositiveLyapunov) {
    container variable(4),parameter(5);
    variable[V_THETA1] = M_PI / 1.1;
    variable[V_THETA2] = M_PI / 1.1;
    variable[V_OMEGA1] = 0.0;
    variable[V_OMEGA2] = 0.0;

    labels_values parameters;
    parameters["l1"]= 0.30;
    parameters["l2"]= 0.30;
    parameters["m1"]= 0.1;
    parameters["m2"]= 0.1;
    parameters["g"]= 9.8;

    DoublePendulumFunction f(parameters);

    RungeKutta model(f, variable, 0.0001); 
    value max_lyapunov = MaxLyapunov<Jacobian_DoublePendulumFunction>(model, parameters, 2000000, 20000);
    EXPECT_TRUE(max_lyapunov > 1);
}


TEST(DoublePendulum, Energy_CrossMethods) {
    container variable(4), parameter(5);
    variable[V_THETA1] = M_PI / 20;
    variable[V_THETA2] = M_PI / 20;
    variable[V_OMEGA1] = 0.1;
    variable[V_OMEGA2] = 0.1;
     parameter[P_L1]= 0.3;
     parameter[P_L2]= 0.3;
     parameter[P_M1]= 0.1;
     parameter[P_M2]= 0.1;
     parameter[P_G]= 9.8;


    labels_values parameters;
    parameters["l1"]= 0.30;
    parameters["l2"]= 0.30;
    parameters["m1"]= 0.1;
    parameters["m2"]= 0.1;
    parameters["g"]= 9.8;

    DoublePendulumFunction f(parameters);

    value delta_t = pow(10, -5);

    value energy0 = DoublePendulumEnergy(variable[V_THETA1],variable[V_THETA2], variable[V_OMEGA1], variable[V_OMEGA2],
            parameter[P_L1], parameter[P_L2], parameter[P_M1], parameter[P_M2], parameter[P_G]);

    RungeKutta model1(f, variable, delta_t); 
    AdamsBashforth model2(f, variable, delta_t); 
    AdamsMoulton model3(f, variable, delta_t); 

    for(size_t i = 0; i < 100000; ++i){
        ++model1;
        ++model2;
        ++model3;
    }

    value energy1 = DoublePendulumEnergy(model1[V_THETA1], model1[V_THETA2], model1[V_OMEGA1], model1[V_OMEGA2],
            parameter[P_L1], parameter[P_L2], parameter[P_M1], parameter[P_M2], parameter[P_G]);
    value energy2 = DoublePendulumEnergy(model2[V_THETA1], model2[V_THETA2], model2[V_OMEGA1], model2[V_OMEGA2],
            parameter[P_L1], parameter[P_L2], parameter[P_M1], parameter[P_M2], parameter[P_G]);
    value energy3 = DoublePendulumEnergy(model3[V_THETA1], model3[V_THETA2], model3[V_OMEGA1], model3[V_OMEGA2],
            parameter[P_L1], parameter[P_L2], parameter[P_M1], parameter[P_M2], parameter[P_G]);

    EXPECT_NEAR(energy1, energy0, 0.00000001);
    EXPECT_NEAR(energy2, energy0, 0.00000001);
    EXPECT_NEAR(energy3, energy0, 0.00000001);

    EXPECT_NEAR(energy1, energy2, 0.00000001);
    EXPECT_NEAR(energy1, energy3, 0.00000001);
    EXPECT_NEAR(energy2, energy3, 0.00000001);
}

TEST(DoublePendulum, EnergyConservation) {
     container variable(4),parameter(5);
     variable[V_THETA1] = M_PI / 2;
     variable[V_THETA2] = M_PI / 2;
     variable[V_OMEGA1] = 0.0;
     variable[V_OMEGA2] = 0.0;

    labels_values parameters;
    parameters["l1"]= 0.30;
    parameters["l2"]= 0.30;
    parameters["m1"]= 0.1;
    parameters["m2"]= 0.1;
    parameters["g"]= 9.8;

     parameter[P_L1]= 0.3;
     parameter[P_L2]= 0.3;
     parameter[P_M1]= 0.1;
     parameter[P_M2]= 0.1;
     parameter[P_G]= 9.8;
     value delta_t = pow(10, -5);
     
    DoublePendulumFunction f(parameters);

    value energy_0 = DoublePendulumEnergy(variable[V_THETA1],variable[V_THETA2], variable[V_OMEGA1], variable[V_OMEGA2],
            parameter[P_L1], parameter[P_L2], parameter[P_M1], parameter[P_M2], parameter[P_G]);


    RungeKutta model1(f,variable, delta_t); 
    AdamsBashforth model2(f, variable, delta_t); 
    AdamsMoulton model3(f, variable, delta_t); 

    value sum1 = 0, sum2 = 0, sum3 = 0;

    for(size_t i = 0; i < pow(10, 6); ++i){
        ++model1; ++model2; ++model3;
    value energy1 = DoublePendulumEnergy(model1[V_THETA1], model1[V_THETA2], model1[V_OMEGA1], model1[V_OMEGA2],
            parameter[P_L1], parameter[P_L2], parameter[P_M1], parameter[P_M2], parameter[P_G]);
    value energy2 = DoublePendulumEnergy(model2[V_THETA1], model2[V_THETA2], model2[V_OMEGA1], model2[V_OMEGA2],
            parameter[P_L1], parameter[P_L2], parameter[P_M1], parameter[P_M2], parameter[P_G]);
    value energy3 = DoublePendulumEnergy(model3[V_THETA1], model3[V_THETA2], model3[V_OMEGA1], model3[V_OMEGA2],
            parameter[P_L1], parameter[P_L2], parameter[P_M1], parameter[P_M2], parameter[P_G]);
        sum1 += fabs(energy1 - energy_0);
        sum2 += fabs(energy2 - energy_0);
        sum3 += fabs(energy3 - energy_0);
    }
    EXPECT_NEAR(sum1, 0.0, 0.0000001);
    EXPECT_NEAR(sum2, 0.0, 0.0000001);
    EXPECT_NEAR(sum3, 0.0, 0.0000001);
}

/*
TEST(ODE, Pendulum_Bifurcation_Diagram) {
        container variable(4),parameter(5);
        variable[V_THETA1] = 0.0;
        variable[V_THETA2] = 0.0;
        variable[V_OMEGA1] = 0.0;
        variable[V_OMEGA2] = 0.0;
        parameter[P_L1]= 0.30;
        parameter[P_L2]= 0.30;
        parameter[P_M1]= 0.10;
        parameter[P_M2]= 0.10;
        parameter[P_G]= 9.8;
     
        value dt=0.0001;
        value coordinate_value=0.0;
        int quadrant=1;
        int coordinate_x=2, coordinate_y=3;
        int time=pow(10,5),transiente=pow(10,3);
        
        std::ofstream file;
        file.open("DP_bifurcation_diagram_simetric.out");
	
        for(value theta = 0.10; theta < M_PI; theta += 0.0005){
            variable[V_THETA1] = theta;
            variable[V_THETA2] = theta;
            AdamsBashforth<DoublePendulumFunction> model(variable, parameter, dt); 
            for(size_t i(0); i < pow(10,5); i++)
                model.next();
            container zero = PhasePlaneSection(model, coordinate_x, coordinate_y, coordinate_value, 8);        
            for(auto item: zero)
                file << (180.0 * theta) / M_PI << " " << item << std::endl;
        }
        file.close();
}
//*/


/*
TEST(ODE, Pendulum_Phase_Plane) {
        container variable(4),parameter(5);
        variable[V_THETA1] = M_PI/5.0;
        variable[V_THETA2] = M_PI/5.0;
        variable[V_OMEGA1] = 0.0;
        variable[V_OMEGA2] = 0.0;
        parameter[P_L1]= 0.30;
        parameter[P_L2]= 0.30;
        parameter[P_M1]= 0.10;
        parameter[P_M2]= 0.10;
        parameter[P_G]= 9.8;
     
        value dt=0.0001;
        value coordinate_value=0.75;
        int quadrant=1;
        int coordinate_x=2, coordinate_y=3;
        value init=0.1,end=M_PI;int n_points=10000;
        int time=4*pow(10,5),transiente=pow(10,3);
        
        DoublePendulumFunction f;

        std::ofstream file;
        file.open("DP_phase_plane.out");
	
        RungeKutta model(f, variable, parameter, 0.0001); 
        container zero = PhasePlaneSection(model, coordinate_x, coordinate_y, coordinate_value, 10);        
        for(auto item: zero)
            file << item << " " << coordinate_value << std::endl;
        file << std::endl;
        file << std::endl;

        RungeKutta orbit(f, variable, parameter, 0.0001); 
        for(size_t i = 0; i <  transiente; ++i)
            orbit.next();
        for(size_t i = 0; i <  (time / 10); ++i){
            file << orbit << std::endl;
            for(size_t j = 0; j <  10; ++j){
                orbit.next();
            }
        }
}
*/

#endif /* value_PENDULUM_H */
