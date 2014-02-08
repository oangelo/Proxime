#ifndef ROSSLER_teste_H
#define ROSSLER_teste_H 

#include "../src/functions/rossler.h"


TEST(rossler, LyapunovSpectrum) {
    labels_values variable;
    variable["x"] = 2.61622;
    variable["y"] = 0.32533;
    variable["z"] = 2.0335135;

    labels_values parameter;
    parameter["a"]= 0.15;
    parameter["b"]= 0.2;
    parameter["c"]= 10.0;
    RosslerFunction function(parameter);

    AdamsBashforth4Th  attractor(function, variable, 0.001);
    Jacobian_RosslerFunction jacobian(parameter);
    for (size_t i = 0; i < 10000; ++i)
        ++attractor;
    LyapunovSpectrum exponent(attractor, jacobian, parameter, pow(10, 4));
    container lambda(exponent(pow(10, 7)));

    EXPECT_NEAR(lambda[0], 0.0890, 0.001);
    EXPECT_NEAR(lambda[1], 0.0000, 0.006);
    EXPECT_NEAR(lambda[2], -9.802, 0.02);
}

TEST(rossler, Jacobian) {
    labels_values variable;
    variable["x"] = 2.61622;
    variable["y"] = 0.32533;
    variable["z"] = 2.0335135;

    labels_values parameter;
    parameter["a"]= 0.15;
    parameter["b"]= 0.2;
    parameter["c"]= 10.0;
    parameter["x"] = 2.61622;
    parameter["y"] = 0.32533;
    parameter["z"] = 2.0335135;


    Jacobian_RosslerFunction jacobian(parameter);
    AdamsBashforth4Th  attractor(jacobian, variable, 0.001);
    for (size_t i = 0; i < 1000; ++i)
        ++attractor;
}

TEST(rossler, PhasePlanePoincareSection) {
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
    RosslerFunction function(parameter);
    AdamsBashforth4Th  model(function, variable,0.001);
    for(size_t i = 0; i <  transiente; ++i){
        ++model;
    }
    container zero = PhasePlaneSection(model, coordinate_x, coordinate_y, coordinate_value, 4);        
    for(size_t i; i < zero.size()-1; ++i)
        EXPECT_NEAR(zero[i],zero[i + 1], 0.00001);
    //Period two ###############################################################
    parameter["c"]= 6.0;
    RosslerFunction function2(parameter);
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
    AdamsBashforth4Th<RosslerFunction>  orbit(variable,parameter,0.001);
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
///*
TEST(rossler, MaxLyapunov) {
    labels_values variable;
    variable["x"] = 2.61622;
    variable["y"] = 0.32533;
    variable["z"] = 2.0335135;

    labels_values parameter;
    parameter["a"]= 0.15;
    parameter["b"]= 0.2;
    parameter["c"]= 10.0;
    RosslerFunction function(parameter);

    AdamsBashforth4Th  attractor(function, variable, 0.001);
    Jacobian_RosslerFunction jacobian(parameter);
    for (size_t i = 0; i < 10000; ++i)
        ++attractor;
    MaxLyapunov exponent(attractor, jacobian, parameter, pow(10, 4));
    EXPECT_NEAR(exponent(pow(10, 7)), 0.0890, 0.0005);
}


/*
TEST(ODE, rossler_bifurcation) {

        container variable(3),parameter(3);
        variable[0] = 2.61622;
        variable[1] = -6.32533;
        variable[2] = 0.0335135;
        parameter[0]= 0.2;
        parameter[1]= 0.1;
        parameter[2]= 5.7;
        value dt=0.001;
        int parameter_index=1;
        value coordinate_value=0;
        int quadrant=-1;
        int coordinate_x=0, coordinate_y=1;
        value init=1.0,end=2.0;int n_points=500;
        int time=100000, transiente=1000000;

        bifurcation<rossler_func>(variable,parameter,dt,
                             parameter_index,coordinate_value,quadrant,
                             coordinate_x,coordinate_y,
                             init,end,n_points,
                             time,transiente,10);

}
//*/
/*
TEST(ODE, rossler_series) {

        container variable(3),parameter(3);
        variable[0] = 2.61622;
        variable[1] = -6.32533;
        variable[2] = 0.0335135;
        parameter[0]= 0.15;
        parameter[1]= 0.20;
        parameter[2]= 10;
	//value c[9]={4,6,8.5,8.7,9,12,12.6,13,18};

        //for(unsigned cont=0;cont<9;cont++){
	
        std::stringstream file;
        file <<  "rosler_series_f3_norbet_bstep"; //<< c[cont];
	//parameter[2]= c[cont];
	AdamsBashforth4Th<rossler_func> attractor(variable,parameter,0.001,"");
	std::ofstream dp_data;
	dp_data.open((file.str()).c_str());
	for(int j=0;j<2*pow(10,6);j++) 
	      attractor.next();
	  for(int i=0;i<pow(10,3);i++){
	    for(int j=0;j<2*pow(10,2);j++) 
	      attractor.next();
	    dp_data << attractor  << endl; 
           //old_var = attractor.get_variable();	  
	  }
	  dp_data.close();
	//}
}
//*/
void rossler_bif_inc_func(container &variable,container &parameter,value increment){
  parameter[1]=increment;
}
/*
TEST(ODE, rossler_bifurcation_MPI) {
        container variable(3),parameter(3);
        variable[0] = 2.61622;
        variable[1] = -6.32533;
        variable[2] = 0.0335135;
        parameter[0]= 0.2;
        parameter[1]= 0.1;
        parameter[2]= 5.7;
        value dt=0.001;
        value coordinate_value=0;
        int quadrant=-1;
        int coordinate_x=0, coordinate_y=1;
        value init=0.0,end=1.0;int n_points=100;
        int time=100000, transiente=1000000;

        MPI_BIFURCATIONS<rossler_func>(variable,parameter,dt,
                             coordinate_value,quadrant,
                             coordinate_x,coordinate_y,
                             init,end,n_points,
                             time,transiente,&rossler_bif_inc_func,g_argc,g_argv);

}
//*/
#endif /* ROSSLER_H */
