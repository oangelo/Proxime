#ifndef LORENZ_H
#define LORENZ_H 

TEST(Lorenz, Lyapunov) {
       
       vector<double> variable(3), parameter(3);
       variable[0] = -3.16;
       variable[1] = -5.31;
       variable[2] = 13.31;
       //Parameters from wikipedia
       parameter[P_SIGMA_LORENZ] = 10;
       parameter[P_GAMMA_LORENZ] = 28;
       parameter[P_BETA_LORENZ] = 8.0 / 3.0;
       
       AdamsMoulton<lorenz_func> attractor(variable, parameter, 0.0001);
       vector<double> lambda = lyapunov<Jacobian_lorenz_func> (attractor, pow(10, 5), pow(10, 5), 1, "teste");
       EXPECT_NEAR(lambda[0], 0.9, 0.1);
       EXPECT_NEAR(lambda[1], 0.0, 0.1);
       EXPECT_NEAR(lambda[2], -14.57, 0.1);
         
}

#endif /* LORENZ_H */
