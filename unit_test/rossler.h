#ifndef ROSSLER_H
#define ROSSLER_H 

TEST(rossler, Lyapunov) {
    type_container variable(3),parameter(3);
    variable[0] = 2.61622;
    variable[1] = -6.32533;
    variable[2] = 0.0335135;
    parameter[0]= 0.15;
    parameter[1]= 0.2;
    parameter[2]= 10.0;
    AdamsBashforth<rossler_func>  attractor(variable,parameter,0.0001);
    type_container lambda = lyapunov<Jacobian_rossler_func>(attractor, pow(10, 6), pow(10, 6), 1, "teste");
    EXPECT_NEAR(lambda[0], 0.09, 0.02);
    EXPECT_NEAR(lambda[1], 0.0, 0.01);
    EXPECT_NEAR(lambda[2], -9.8, 0.4);
    for(auto iten: lambda)
        std::cout << iten << std::endl;
}
/*
TEST(ODE, rossler_lyapunov) {
      
  type_container variable(3),parameter(3);
  variable[0] = 2.61622;
  variable[1] = -6.32533;
  variable[2] = 0.0335135;
  parameter[0]= 0.1;
  parameter[1]= 0.1;
  parameter[2]= 9;
  AdamsBashforth<rossler_func>  attractor(variable,parameter,0.001);
  lyapunov<Jacobian_rossler_func>(attractor,2*pow(10,7),2*pow(10,6),1000,"rossler9");  

}
//*/
/*
TEST(ODE, rossler_bifurcation) {

        type_container variable(3),parameter(3);
        variable[0] = 2.61622;
        variable[1] = -6.32533;
        variable[2] = 0.0335135;
        parameter[0]= 0.2;
        parameter[1]= 0.1;
        parameter[2]= 5.7;
        type_data dt=0.001;
        int parameter_index=1;
        type_data coordinate_value=0;
        int quadrant=-1;
        int coordinate_x=0, coordinate_y=1;
        type_data init=1.0,end=2.0;int n_points=500;
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

        type_container variable(3),parameter(3);
        variable[0] = 2.61622;
        variable[1] = -6.32533;
        variable[2] = 0.0335135;
        parameter[0]= 0.15;
        parameter[1]= 0.20;
        parameter[2]= 10;
	//type_data c[9]={4,6,8.5,8.7,9,12,12.6,13,18};

        //for(unsigned cont=0;cont<9;cont++){
	
        std::stringstream file;
        file <<  "rosler_series_f3_norbet_bstep"; //<< c[cont];
	//parameter[2]= c[cont];
	AdamsBashforth<rossler_func> attractor(variable,parameter,0.001,"");
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
void rossler_bif_inc_func(type_container &variable,type_container &parameter,type_data increment){
  parameter[1]=increment;
}

/*
TEST(ODE, rossler_bifurcation_MPI) {
        type_container variable(3),parameter(3);
        variable[0] = 2.61622;
        variable[1] = -6.32533;
        variable[2] = 0.0335135;
        parameter[0]= 0.2;
        parameter[1]= 0.1;
        parameter[2]= 5.7;
        type_data dt=0.001;
        type_data coordinate_value=0;
        int quadrant=-1;
        int coordinate_x=0, coordinate_y=1;
        type_data init=0.0,end=1.0;int n_points=100;
        int time=100000, transiente=1000000;

        MPI_BIFURCATIONS<rossler_func>(variable,parameter,dt,
                             coordinate_value,quadrant,
                             coordinate_x,coordinate_y,
                             init,end,n_points,
                             time,transiente,&rossler_bif_inc_func,g_argc,g_argv);

}
//*/

#endif /* ROSSLER_H */
