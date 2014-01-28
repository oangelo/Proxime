#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include "gtest/gtest.h"

//#include "../src/functions/simple_pendulum.h"
#include "../src/functions/double_pendulum.h"
//#include "../src/functions/rossler.h"
//#include "../src/functions/lorenz.h"
#include "../src/numerical_integration/adams_moulton.h"
#include "../src/numerical_integration/adams_moulton.h"
#include "../src/bifurcation_diagram.h"
//#include "../src/lyapunov.h"

//#include "rossler.h"
#include "functions.h"
#include "double_pendulum.h"
//#include "simple_pendulum.h"
//#include "t_map.h"
//#include "lorenz.h"

TEST(ODE, wrong_access) {
    container variable(4);
    variable[V_THETA1] = 5 * M_PI / 10;
    variable[V_THETA2] = 5 * M_PI / 10;
    variable[V_OMEGA1] = 0.0;
    variable[V_OMEGA2] = 0.0;

    labels_and_values parameters;
    parameters["l1"] = 0.30;
    parameters["l2"] = 0.30;
    parameters["m1"] = 0.10;
    parameters["m2"] = 0.10;
    parameters["g"] = 9.8;


    DoublePendulumFunction dp_function(parameters);
    AdamsMoulton pendulo(dp_function, variable, 0.00001);
    EXPECT_THROW(pendulo[4], Index_error);
    EXPECT_THROW(pendulo.get_variable(4), Index_error);
}

TEST(ODE, bad_integration) {
    value aux = 0;
    container variable(4), parameter(5);
    variable[V_THETA1] = 1;
    variable[V_THETA2] = 1 / aux;
    variable[V_OMEGA1] = 0.0;
    variable[V_OMEGA2] = 0.0;

    labels_and_values parameters;
    parameters["l1"] = 0.30;
    parameters["l2"] = 0.30;
    parameters["m1"] = 0.10;
    parameters["m2"] = 0.10;
    parameters["g"] = 9.8;

    DoublePendulumFunction dp_function(parameters);
    EXPECT_THROW(AdamsMoulton pendulo(dp_function, variable, 0.00001);, Value_error);

}

TEST(ODE, assigment) {
    container variable(4), parameter(5);
    variable[V_THETA1] = 1.1;
    variable[V_THETA2] = 2.2;
    variable[V_OMEGA1] = 3.3;
    variable[V_OMEGA2] = 4.4;

    labels_and_values parameters;
    parameters["l1"] = 0.30;
    parameters["l2"] = 0.30;
    parameters["m1"] = 0.10;
    parameters["m2"] = 0.10;
    parameters["g"] = 9.8;


    DoublePendulumFunction dp_function(parameters);
    AdamsMoulton pendulo(dp_function, variable, 0.00001);

    EXPECT_NEAR(pendulo[V_THETA1], 1.1, 0.1);
    EXPECT_NEAR(pendulo[V_THETA2], 2.2, 0.1);
    EXPECT_NEAR(pendulo[V_OMEGA1], 3.3, 0.1);
    EXPECT_NEAR(pendulo[V_OMEGA2], 4.4, 0.1);
}

TEST(ODE, size) {

    container variable(4), parameter(5);
    variable[V_THETA1] = 5 * M_PI / 10;
    variable[V_THETA2] = 5 * M_PI / 10;
    variable[V_OMEGA1] = 0.0;
    variable[V_OMEGA2] = 0.0;

    labels_and_values parameters;
    parameters["l1"] = 0.30;
    parameters["l2"] = 0.30;
    parameters["m1"] = 0.10;
    parameters["m2"] = 0.10;
    parameters["g"] = 9.8;

    DoublePendulumFunction dp_function(parameters);
    AdamsMoulton pendulo(dp_function, variable, 0.00001);

    EXPECT_EQ(pendulo.size_variable(), 4);

}

int main( int argc , char * argv[]) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
