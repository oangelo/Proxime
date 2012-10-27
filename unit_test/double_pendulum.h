#ifndef DOUBLE_PENDULUM_H
#define DOUBLE_PENDULUM_H 
void dp_bif_inc_func(std::vector<double> &variable,std::vector<double> &parameter,double increment){
      variable[V_THETA1] = increment;
      variable[V_THETA2] = 0;    
}

TEST(ODE, Double_pendulum_bifurcation) {
    
        vector<double> variable(4),parameter(5);
        variable[V_THETA1] = M_PI/2.0;
        variable[V_THETA2] = M_PI/2.0;
        variable[V_OMEGA1] = 0.0;
        variable[V_OMEGA2] = 0.0;
        parameter[P_L1]= 0.30;
        parameter[P_L2]= 0.30;
        parameter[P_M1]= 0.10;
        parameter[P_M2]= 0.10;
        parameter[P_G]= 9.8;
     
        double dt=0.0001;
        double coordinate_value=0.75;
        int quadrant=1;
        int coordinate_x=2, coordinate_y=3;
        double init=0.1,end=M_PI;int n_points=10000;
        int time=pow(10,5),transiente=pow(10,3);
        
	
       /*    
	MPI_BIFURCATIONS<double_pendulum_func>(variable,parameter,dt,
                             coordinate_value,quadrant,
                             coordinate_x,coordinate_y,
                             init,end,n_points,
                             time,transiente,&dp_bif_inc_func,g_argc,g_argv);        
	*/
	//*

        for(int angle = 1;angle < 6;angle++){
            variable[V_THETA1] = angle * M_PI / 10.0;
            variable[V_THETA2] = 0;
            std::stringstream file;
            file <<  "dp_ns_series_" << angle <<"_.out";

            AdamsBashforth<double_pendulum_func> attractor(variable,parameter,0.0001);
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
        //double_pendulum dp_attractor(variable,parameter,0.001);



        //portrail<double_pendulum>(dp_attractor,1000000,1000000,0,1,0,-1);
        //gen_atractor(dp_attractor,100000,1000000);
    //*/


}

/*
TEST(ODE, Double_Pendulum_Lyapunov) {
    vector<double> variable(4),parameter(8);
    variable[V_THETA1] = 120*M_PI/180.0;
    variable[V_THETA2] = -120*M_PI/180.0;
    variable[V_OMEGA1] = 0.0;
    variable[V_OMEGA2] = 0.0;
    parameter[P_L1]= 0.057;
    parameter[P_L2]= 0.025;
    parameter[P_M1]= 32.243767+13.05691;
    parameter[P_M2]= 87.7+12.576;
    parameter[P_G]= 9.799403;

    

    lyapunov<jacobian_double_pendulum_func>(attractor,pow(10,6),pow(10,6),1,"test");
}
//*/
/*
TEST(ODE, Double_Pendulum_Lyapunov_MAX_MPI) {
  
    int total_steps=100;
  
  int MPI_MTAG1=1,MPI_MTAG2=2;
  int MPI_MYID , MPI_NUMPROCS,MPI_ISLAVE,MPI_SLAVE_STEPS;
  
  MPI_Status MPI_STATUS ;
  MPI_Init (&argc ,&argv );
  MPI_Comm_size ( MPI_COMM_WORLD ,&MPI_NUMPROCS);
  MPI_Comm_rank ( MPI_COMM_WORLD ,&MPI_MYID);
  
  MPI_SLAVE_STEPS=total_steps/MPI_NUMPROCS;
  double lyapunov[MPI_SLAVE_STEPS];
  //The Slave program 
  if(MPI_MYID!=0){
    double angle;
    int cont=0;
    for (int step = (MPI_MYID-1)*MPI_SLAVE_STEPS; step < MPI_MYID*MPI_SLAVE_STEPS; step++){
	    angle = ((double)step)/100.0;
            vector<double> variable(4), parameter(8);
            variable[V_THETA1] = angle * M_PI;
            variable[V_THETA2] = angle*M_PI;
            variable[V_OMEGA1] = 0.0;
            variable[V_OMEGA2] = 0.0;
            parameter[P_L1] = 0.30;
            parameter[P_L2] = 0.30;
            parameter[P_M1] = 2*0.10;
            parameter[P_M2] = 0.10;
            parameter[P_G] = 9.8;
            double_pendulum pendulo(variable, parameter, 0.00001);       
            lyapunov[cont] = 
lyapunov_max<Jacob_double_pendulum>(pendulo,640000000,200000);
	    cont++;
    }
    MPI_Send(&lyapunov,MPI_SLAVE_STEPS, MPI_DOUBLE_PRECISION , 0, MPI_MTAG1, MPI_COMM_WORLD);
  }else{
	//The Master program
	ofstream data_lyapunov;
	std::string Filename = 
"/nfs/fiscomp/angelo/lyapunov_angles_simetric_2m2.out";
	data_lyapunov.open(Filename.c_str());
	double angle;
	int cont;
  	for ( MPI_ISLAVE =1; MPI_ISLAVE < MPI_NUMPROCS ; MPI_ISLAVE++){
		MPI_Recv(&lyapunov,MPI_SLAVE_STEPS, MPI_DOUBLE_PRECISION, MPI_ISLAVE , MPI_MTAG1, MPI_COMM_WORLD,&MPI_STATUS);
		if (data_lyapunov.is_open()) {
		  cont=0;
		  for (int step = (MPI_ISLAVE-1)*MPI_SLAVE_STEPS; step < MPI_ISLAVE*MPI_SLAVE_STEPS; step++){
			angle = ((double)step)/10.0;
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
TEST(ODE, Double_pendulum_bifurcation) {
    
        vector<double> variable(4),parameter(5);
        variable[V_THETA1] = M_PI/2.0;
        variable[V_THETA2] = M_PI/2.0;
        variable[V_OMEGA1] = 0.0;
        variable[V_OMEGA2] = 0.0;
        parameter[P_L1]= 0.30;
        parameter[P_L2]= 0.30;
        parameter[P_M1]= 0.10;
        parameter[P_M2]= 0.10;
        parameter[P_G]= 9.8;
     
        double dt=0.0001;
        double coordinate_value=0.75;
        int quadrant=1;
        int coordinate_x=2, coordinate_y=3;
        double init=0.1,end=M_PI;int n_points=10000;
        int time=pow(10,5),transiente=pow(10,3);
        
	
        for(int angle = 1;angle < 6;angle++){
            variable[V_THETA1] = angle * M_PI / 10.0;
            variable[V_THETA2] = 0;
            std::stringstream file;
            file <<  "dp_ns_series_" << angle <<"_.out";

            AdamsBashforth<double_pendulum_func> attractor(variable,parameter,0.0001);
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

TEST(ODE, Double_pendulum_ENERGY) {
   /*
     vector<double> variable(4),parameter(5);
     variable[V_THETA1] = 4*M_PI/10;
     variable[V_THETA2] = 4*M_PI/10;
     variable[V_OMEGA1] = 0.0;
     variable[V_OMEGA2] = 0.0;
     parameter[P_L1]= 0.30;
     parameter[P_L2]= 0.30;
     parameter[P_M1]= 0.10;
     parameter[P_M2]= 0.10;
     parameter[P_G]= 9.8;

     cout.precision(10);
     parameter[P_L1]= 0.30;
     parameter[P_L2]= 0.30;
     parameter[P_M1]= 0.10;
     parameter[P_M2]= 0.10;
     parameter[P_G]= 9.8;

     double_pendulum dp_attractor(variable,parameter,0.00001);
     std::ofstream file;
     std::string Filename ="./energy/energia_0.4PI_Rk_0.00001.out";
     file.open(Filename.c_str());
     double Eref=2*0.882;
     //dp_attractor.next();
     double E0=dp_attractor.get_energy();
     cout << "# E0 =  " << endl;
     for(int i=0;i<10000;i++){
         for(int j=0;j<10000;j++){
             dp_attractor.next();
         }
         file  << dp_attractor.get_t()  << " " << fabs((dp_attractor.get_energy()-E0)/Eref) << endl;
     }
     file.close();
     //*/

}

#endif /* DOUBLE_PENDULUM_H */
