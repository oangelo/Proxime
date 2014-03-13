#ifndef ROSSLER_teste_H
#define ROSSLER_teste_H 

#include "../src/functions/rossler.h"


TEST(Rossler, LyapunovSpectrum) {
    labels_values variable;
    variable["x"] = 2.61622;
    variable["y"] = 0.32533;
    variable["z"] = 2.0335135;

    labels_values parameter;
    parameter["a"]= 0.15;
    parameter["b"]= 0.2;
    parameter["c"]= 10.0;
    Rossler function(parameter);

    AdamsBashforth4Th  attractor(function, variable, 0.001);
    Jacobian_Rossler jacobian(parameter);
    for (size_t i = 0; i < 10000; ++i)
        ++attractor;
    LyapunovSpectrum exponent(attractor, jacobian, parameter, pow(10, 4));
    container lambda(exponent(pow(10, 7)));

    EXPECT_NEAR(lambda[0], 0.0890, 0.001);
    EXPECT_NEAR(lambda[1], 0.0000, 0.001);
    EXPECT_NEAR(lambda[2], -9.802, 0.001);
}

TEST(Rossler, MaxLyapunov) {
    labels_values variable;
    variable["x"] = 2.61622;
    variable["y"] = 0.32533;
    variable["z"] = 2.0335135;

    labels_values parameter;
    parameter["a"]= 0.15;
    parameter["b"]= 0.2;
    parameter["c"]= 10.0;
    Rossler function(parameter);

    AdamsBashforth4Th  attractor(function, variable, 0.001);
    Jacobian_Rossler jacobian(parameter);
    for (size_t i = 0; i < 10000; ++i)
        ++attractor;
    MaxLyapunov exponent(attractor, jacobian, parameter, pow(10, 4));
    EXPECT_NEAR(exponent(pow(10, 7)), 0.0890, 0.0005);
}


TEST(Rossler, PhasePlanePoincareSection) {
    labels_values variable;
    variable["x"] = 2.61622;
    variable["y"] = 2.32533;
    variable["z"] = 2.0335135;

    labels_values parameter;
    parameter["a"]= 0.1;
    parameter["b"]= 0.1;
    parameter["c"]= 4.0;

    size_t coordinate_x=0, coordinate_y=1;
    unsigned transiente=pow(10,5);
    value coordinate_value = 0;

    //Period one ###############################################################
    parameter["c"]= 4.0;
    Rossler function(parameter);
    AdamsBashforth4Th  model(function, variable,0.001);
    for(size_t i = 0; i <  transiente; ++i){
        ++model;
    }
    container zero = PhasePlaneSection(model, coordinate_x, coordinate_y, coordinate_value, 4);        
    for(size_t i; i < zero.size()-1; ++i)
        EXPECT_NEAR(zero[i],zero[i + 1], 0.00001);
    //Period two ###############################################################
    parameter["c"]= 6.0;
    Rossler function2(parameter);
    AdamsBashforth4Th  model2(function2, variable, 0.001);
    for(size_t i = 0; i <  transiente; ++i){
        ++model2;
    }
    zero = PhasePlaneSection(model2, coordinate_x, coordinate_y, coordinate_value, 4);        
    for(size_t i; i < zero.size() - 2; ++i)
        EXPECT_NEAR(zero[i],zero[i + 2], 0.01);


    /* Prints the plane on a file
    std::ofstream file;
    file.open("rossler_phase_plane.out");
    for(auto item: zero)
        file << item << " " << coordinate_value << std::endl;
    file << std::endl;
    file << std::endl;
    AdamsBashforth4Th<Rossler>  orbit(variable,parameter,0.001);
    for(size_t i = 0; i <  transiente; ++i)
        orbit.next();
    for(size_t i = 0; i <  (time / 10); ++i){
        file << orbit << std::endl;
        for(size_t j = 0; j <  10; ++j){
            orbit.next();
        }
    }
    */
}

#endif /* ROSSLER_H */
