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
#include "../src/lyapunov.h"

#include "rossler.h"
#include "functions.h"
#include "double_pendulum.h"
//#include "simple_pendulum.h"
//#include "t_map.h"
//#include "lorenz.h"

TEST(ODE, wrong_access) {
    labels_values variable;
    variable["theta1"] = M_PI / 10;
    variable["theta2"] = M_PI / 10;
    variable["omega1"] = 0.0;
    variable["omega2"] = 0.0;


    labels_values parameters;
    parameters["l1"] = 0.30;
    parameters["l2"] = 0.30;
    parameters["m1"] = 0.10;
    parameters["m2"] = 0.10;
    parameters["g"] = 9.8;


    DoublePendulumFunction dp_function(parameters);
    AdamsMoulton4Th pendulo(dp_function, variable, 0.00001);
    EXPECT_THROW(pendulo[4], Index_error);
}

TEST(ODE, bad_integration) {
    labels_values variable;
    variable["theta1"] = M_PI / 0.0;
    variable["theta2"] = M_PI / 0.0;
    variable["omega1"] = 0.0;
    variable["omega2"] = 0.0;

    labels_values parameters;
    parameters["l1"] = 0.30;
    parameters["l2"] = 0.30;
    parameters["m1"] = 0.10;
    parameters["m2"] = 0.10;
    parameters["g"] = 9.8;

    DoublePendulumFunction dp_function(parameters);
    EXPECT_THROW(AdamsMoulton4Th pendulo(dp_function, variable, 0.00001);, Value_error);

}

TEST(ODE, assigment) {

    labels_values variable;
    variable["theta1"] = 1.1;
    variable["theta2"] = 2.2;
    variable["omega1"] = 3.3;
    variable["omega2"] = 4.4;

    labels_values parameters;
    parameters["l1"] = 0.30;
    parameters["l2"] = 0.30;
    parameters["m1"] = 0.10;
    parameters["m2"] = 0.10;
    parameters["g"] = 9.8;


    DoublePendulumFunction dp_function(parameters);
    AdamsMoulton4Th pendulo(dp_function, variable, 0.00001);

    EXPECT_NEAR(pendulo["theta1"], 1.1, 0.1);
    EXPECT_NEAR(pendulo["theta2"], 2.2, 0.1);
    EXPECT_NEAR(pendulo["omega1"], 3.3, 0.1);
    EXPECT_NEAR(pendulo["omega2"], 4.4, 0.1);
}

TEST(ODE, ObjectFactory) {

    labels_values variable;
    variable["theta1"] = M_PI / 10;
    variable["theta2"] = M_PI / 10;
    variable["omega1"] = 0.0;
    variable["omega2"] = 0.0;

    labels_values parameters;
    parameters["l1"] = 0.30;
    parameters["l2"] = 0.30;
    parameters["m1"] = 0.10;
    parameters["m2"] = 0.10;
    parameters["g"] = 9.8;

    DoublePendulumFunction dp_function(parameters);
    AdamsMoulton4Th pendula(dp_function, variable, 0.00001);

    AdamsMoulton4Th* ptr_clone(pendula.Clone()); 
    AdamsMoulton4Th& clone(*ptr_clone);

    AdamsMoulton4Th* ptr_new(pendula.Create(dp_function, variable, 0.00001)); 
    AdamsMoulton4Th& other(*ptr_new);

    for(size_t i(0); i < 100; ++i){
        ++pendula;
        ++clone;
        ++other;
    }

    for(size_t i(0); i < pendula.size(); ++i){
        EXPECT_EQ(pendula[i], clone[i]);
        EXPECT_EQ(pendula[i], other[i]);
    }

}

TEST(ODE, size) {

    labels_values variable;
    variable["theta1"] = M_PI / 10;
    variable["theta2"] = M_PI / 10;
    variable["omega1"] = 0.0;
    variable["omega2"] = 0.0;

    labels_values parameters;
    parameters["l1"] = 0.30;
    parameters["l2"] = 0.30;
    parameters["m1"] = 0.10;
    parameters["m2"] = 0.10;
    parameters["g"] = 9.8;

    DoublePendulumFunction dp_function(parameters);
    AdamsMoulton4Th pendulo(dp_function, variable, 0.00001);

    EXPECT_EQ(pendulo.size(), 4);

}

int main( int argc , char * argv[]) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
