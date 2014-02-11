#ifndef LORENZ_test_H
#define LORENZ_test_H 

#include "../src/functions/lorenz.h"


TEST(Lorenz, LyapunovSpectrum) {
    labels_values variable;
    variable["x"] = 2.61622;
    variable["y"] = 0.32533;
    variable["z"] = 2.0335135;

    labels_values parameter;
    parameter["sigma"]= 10;
    parameter["gamma"]= 28;
    parameter["beta"] = 8.0 / 3.0;
    LorenzFunction function(parameter);

    AdamsBashforth4Th  attractor(function, variable, 0.001);
    Jacobian_LorenzFunction jacobian(parameter);
    for (size_t i = 0; i < 10000; ++i)
        ++attractor;
    LyapunovSpectrum exponent(attractor, jacobian, parameter, pow(10, 4));
    container lambda(exponent(pow(10, 6)));

    EXPECT_NEAR(lambda[0], 0.907, 0.01);
    EXPECT_NEAR(lambda[1], 0.0000, 0.001);
    EXPECT_NEAR(lambda[2], -14.574, 0.01);

}

TEST(Lorenz, MaxLyapunov) {
    labels_values variable;
    variable["x"] = 2.61622;
    variable["y"] = 0.32533;
    variable["z"] = 2.0335135;

    labels_values parameter;
    parameter["sigma"]= 10;
    parameter["gamma"]= 28;
    parameter["beta"] = 8.0 / 3.0;
    LorenzFunction function(parameter);

    AdamsBashforth4Th  attractor(function, variable, 0.001);
    Jacobian_LorenzFunction jacobian(parameter);
    for (size_t i = 0; i < 10000; ++i)
        ++attractor;
    MaxLyapunov exponent(attractor, jacobian, parameter, pow(10, 4));
    EXPECT_NEAR(exponent(pow(10, 7)), 0.907, 0.01);
}

#endif /* LORENZ_H */
