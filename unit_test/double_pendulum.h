#ifndef type_data_PENDULUM_H
#define type_data_PENDULUM_H 

#include "../src/functions/double_pendulum.h"

TEST(DoublePendulum, Energy) {
     type_container variable(4),parameter(5);
     variable[V_THETA1] = 4*M_PI/10;
     variable[V_THETA2] = 4*M_PI/10;
     variable[V_OMEGA1] = 0.0;
     variable[V_OMEGA2] = 0.0;
     parameter[P_L1]= 0.30;
     parameter[P_L2]= 0.30;
     parameter[P_M1]= 0.10;
     parameter[P_M2]= 0.10;
     parameter[P_G]= 9.8;

     parameter[P_L1]= 0.30;
     parameter[P_L2]= 0.30;
     parameter[P_M1]= 0.10;
     parameter[P_M2]= 0.10;
     parameter[P_G]= 9.8;

     double delta_t = 0.0001;
     

    RungeKutta<DoublePendulumFunction> model(variable, parameter, delta_t); 
    double energy_0 = DoublePendulumEnergy(variable[V_THETA1],variable[V_THETA2], variable[V_OMEGA1], variable[V_OMEGA1],
            parameter[P_L1], parameter[P_L2], parameter[P_M1], parameter[P_M2], parameter[P_G]);
    for(size_t i = 0; i < 1000000; ++i){
        model.next();
        double energy = DoublePendulumEnergy(model[V_THETA1], model[V_THETA2], model[V_OMEGA1], model[V_OMEGA2],
        parameter[P_L1], parameter[P_L2], parameter[P_M1], parameter[P_M2], parameter[P_G]);
        EXPECT_NEAR((energy - energy_0) / energy_0, 0.0, 0.0001);
    }
}

TEST(DoublePendulum, Lyapunov) {
    std::vector<double> variable(4),parameter(5);
    variable[V_THETA1] = 5*M_PI/10.0;
    variable[V_THETA2] = 5*M_PI/10.0;
    variable[V_OMEGA1] = 0.0;
    variable[V_OMEGA2] = 0.0;
    parameter[P_L1]= 0.30;
    parameter[P_L2]= 0.30;
    parameter[P_M1]= 0.10;
    parameter[P_M2]= 0.10;
    parameter[P_G]= 9.8;

    RungeKutta<DoublePendulumFunction> model(variable, parameter, 0.0001); 
    double max_lyapunov = MaxLyapunov<Jacobian_DoublePendulumFunction>(model, 2000000, 20000);
    EXPECT_NEAR(max_lyapunov, 1.7, 0.2);
}


/*
void dp_bif_inc_func(std::type_container &variable,std::type_container &parameter,type_data increment){
      variable[V_THETA1] = increment;
      variable[V_THETA2] = 0;    
}
*/

TEST(ODE, type_data_pendulum_bifurcation) {
/*    
        type_container variable(4),parameter(5);
        variable[V_THETA1] = M_PI/2.0;
        variable[V_THETA2] = M_PI/2.0;
        variable[V_OMEGA1] = 0.0;
        variable[V_OMEGA2] = 0.0;
        parameter[P_L1]= 0.30;
        parameter[P_L2]= 0.30;
        parameter[P_M1]= 0.10;
        parameter[P_M2]= 0.10;
        parameter[P_G]= 9.8;
     
        type_data dt=0.0001;
        type_data coordinate_value=0.75;
        int quadrant=1;
        int coordinate_x=2, coordinate_y=3;
        type_data init=0.1,end=M_PI;int n_points=10000;
        int time=pow(10,5),transiente=pow(10,3);
        
	
        AttractorCrossAxis(variable,parameter,dt,
                             coordinate_value,quadrant,
                             coordinate_x,coordinate_y,
                             init,end,n_points,
                             time,transiente,&dp_bif_inc_func,g_argc,g_argv);        

        for(int angle = 1;angle < 6;angle++){
            variable[V_THETA1] = angle * M_PI / 10.0;
            variable[V_THETA2] = 0;
            std::stringstream file;
            file <<  "dp_ns_series_" << angle <<"_.out";

            AdamsBashforth<type_data_pendulum_func> attractor(variable,parameter,0.0001);
            std::ofstream dp_data;
            std::string file_name=file.str();
            dp_data.open(file_name.c_str());
            for(int i=0;i<pow(10,4);i++)
                for(int j=0;j<pow(10,3);j++) 
                    attractor.next();
 
            for(int i=0;i<pow(10,4);i++){
                for(int j=0;j<pow(10,3);j++) 
                    attractor.next();
                dp_data << attractor  << endl; 
                //old_var = attractor.get_variable();	  
            }
            dp_data.close();
        }
	//*/
        //type_data_pendulum dp_attractor(variable,parameter,0.001);



        //portrail<type_data_pendulum>(dp_attractor,1000000,1000000,0,1,0,-1);
        //gen_atractor(dp_attractor,100000,1000000);
    //*/


}

//*
//*/
/*
TEST(ODE, type_data_Pendulum_Lyapunov_MAX_MPI) {
  
    int total_steps=100;
  
  int MPI_MTAG1=1,MPI_MTAG2=2;
  int MPI_MYID , MPI_NUMPROCS,MPI_ISLAVE,MPI_SLAVE_STEPS;
  
  MPI_Status MPI_STATUS ;
  MPI_Init (&argc ,&argv );
  MPI_Comm_size ( MPI_COMM_WORLD ,&MPI_NUMPROCS);
  MPI_Comm_rank ( MPI_COMM_WORLD ,&MPI_MYID);
  
  MPI_SLAVE_STEPS=total_steps/MPI_NUMPROCS;
  type_data lyapunov[MPI_SLAVE_STEPS];
  //The Slave program 
  if(MPI_MYID!=0){
    type_data angle;
    int cont=0;
    for (int step = (MPI_MYID-1)*MPI_SLAVE_STEPS; step < MPI_MYID*MPI_SLAVE_STEPS; step++){
	    angle = ((type_data)step)/100.0;
            type_container variable(4), parameter(8);
            variable[V_THETA1] = angle * M_PI;
            variable[V_THETA2] = angle*M_PI;
            variable[V_OMEGA1] = 0.0;
            variable[V_OMEGA2] = 0.0;
            parameter[P_L1] = 0.30;
            parameter[P_L2] = 0.30;
            parameter[P_M1] = 2*0.10;
            parameter[P_M2] = 0.10;
            parameter[P_G] = 9.8;
            type_data_pendulum pendulo(variable, parameter, 0.00001);       
            lyapunov[cont] = 
lyapunov_max<Jacob_type_data_pendulum>(pendulo,640000000,200000);
	    cont++;
    }
    MPI_Send(&lyapunov,MPI_SLAVE_STEPS, MPI_type_data_PRECISION , 0, MPI_MTAG1, MPI_COMM_WORLD);
  }else{
	//The Master program
	ofstream data_lyapunov;
	std::string Filename = 
"/nfs/fiscomp/angelo/lyapunov_angles_simetric_2m2.out";
	data_lyapunov.open(Filename.c_str());
	type_data angle;
	int cont;
  	for ( MPI_ISLAVE =1; MPI_ISLAVE < MPI_NUMPROCS ; MPI_ISLAVE++){
		MPI_Recv(&lyapunov,MPI_SLAVE_STEPS, MPI_type_data_PRECISION, MPI_ISLAVE , MPI_MTAG1, MPI_COMM_WORLD,&MPI_STATUS);
		if (data_lyapunov.is_open()) {
		  cont=0;
		  for (int step = (MPI_ISLAVE-1)*MPI_SLAVE_STEPS; step < MPI_ISLAVE*MPI_SLAVE_STEPS; step++){
			angle = ((type_data)step)/10.0;
			data_lyapunov << angle << " ";
			data_lyapunov << lyapunov[cont] << endl;
			cont++;
		  }
		}else{ cout << "Unable to open file"; }
      }
      data_lyapunov.close();
  }
    
  MPI_Finalize();
  
}
//*/
/*
TEST(ODE, type_data_pendulum_bifurcation) {
    
        type_container variable(4),parameter(5);
        variable[V_THETA1] = M_PI/2.0;
        variable[V_THETA2] = M_PI/2.0;
        variable[V_OMEGA1] = 0.0;
        variable[V_OMEGA2] = 0.0;
        parameter[P_L1]= 0.30;
        parameter[P_L2]= 0.30;
        parameter[P_M1]= 0.10;
        parameter[P_M2]= 0.10;
        parameter[P_G]= 9.8;
     
        type_data dt=0.0001;
        type_data coordinate_value=0.75;
        int quadrant=1;
        int coordinate_x=2, coordinate_y=3;
        type_data init=0.1,end=M_PI;int n_points=10000;
        int time=pow(10,5),transiente=pow(10,3);
        
	
        for(int angle = 1;angle < 6;angle++){
            variable[V_THETA1] = angle * M_PI / 10.0;
            variable[V_THETA2] = 0;
            std::stringstream file;
            file <<  "dp_ns_series_" << angle <<"_.out";

            AdamsBashforth<type_data_pendulum_func> attractor(variable,parameter,0.0001);
            std::ofstream dp_data;
            std::string file_name=file.str();
            dp_data.open(file_name.c_str());
            for(int i=0;i<pow(10,3);i++)
                    attractor.next();
 
            for(int i=0;i<5*pow(10,3);i++){
                for(int j=0;j<pow(10,3);j++) 
                    attractor.next();
                dp_data << attractor  << endl; 
            }
            dp_data.close();
        }
}
*/

#endif /* type_data_PENDULUM_H */
