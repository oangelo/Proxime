#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <math.h>
#include <string>
#include <stdlib.h>

#include <mpi.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>

#include "bifurcation_diagram.h"
#include "Lyapunov.h"
#include "recurrence_relation.h"

using namespace std;

int g_argc; 
char ** g_argv;
/*
void FFT(double *data, int n,std::string file_name){
  int i;
  
  double x[2*n];
  gsl_complex_packed_array coefficient = x;
  
  gsl_fft_real_wavetable * real;
  gsl_fft_halfcomplex_wavetable * hc;
  gsl_fft_real_workspace * work;
  
  
  ofstream data_file;
  ofstream fft_file;
  ofstream invert_fft_file; 
  
  std::string Filename1 = "data_FFT_"+file_name+".out";
  std::string Filename2 = "FFT_"+file_name+".out";
  std::string Filename3 = "invert_FFT_"+file_name+".out";
  
  data_file.open(Filename1.c_str());
  fft_file.open(Filename2.c_str());
  invert_fft_file.open(Filename3.c_str()); 
  
  i=0;
  
  
  for (i = 0; i < n; i++)
  {
    data_file << i << " " << data[i] << endl;
  }
  printf ("\n");
  
  work = gsl_fft_real_workspace_alloc (n);
  real = gsl_fft_real_wavetable_alloc (n);
  
  gsl_fft_real_transform (data, 1, n, 
			  real, work);
  
  gsl_fft_real_wavetable_free (real);
  
  gsl_fft_halfcomplex_unpack (data,coefficient,1,n);
  
  for (i = 0; i < n; i=i+2)
  {
    fft_file << (i+n)/2 << " " << sqrt(coefficient[i]*coefficient[i]+coefficient[i+1]*coefficient[i+1]) << endl;
  }   
  for (i = n; i < 2*n; i=i+2)
  {
    fft_file << (i-n)/2 << " " << sqrt(coefficient[i]*coefficient[i]+coefficient[i+1]*coefficient[i+1]) << endl;
  }
  
  
  hc = gsl_fft_halfcomplex_wavetable_alloc (n);
  
  gsl_fft_halfcomplex_inverse (data, 1, n, 
			       hc, work);
  gsl_fft_halfcomplex_wavetable_free (hc);
  
  for (i = 0; i < n; i++)
  {
    invert_fft_file << i << " " << data[i] << endl;
  } 
  gsl_fft_real_workspace_free (work);
  
  data_file.close();
  fft_file.close();
  invert_fft_file.close(); 
  
  
}  


void menu(int argc , char * argv[]){
  
  // cout << "The name used to start the program: " << argv[ 0 ]  << "\nArguments are:\n";
  //for (int n = 1; n < argc; n++)
  // cout << setw( 2 ) << n << ": " << argv[ n ] << '\n';
  string measure=argv[1];
  string model;
  
  if(model=="-h"){
    cout << "measure[lyapunov] model[lorenz, Rossler] method[rk, ab, am] archive_name  parameter[0] parameter[1] parameter[2] dt steps[10^(int)] trasient[10^(int)]  " << endl;
    cout << "measure[energy] model[Double_Pendulum] method[rk, ab, am] archive_namerm  angle1*PI angle2*PI dt steps*10^4" << endl;
    cout << "./chaos GL_pendulum  0.7 0.7" << endl;
}
if(measure=="lyapunov")
{
  model=argv[2];
  if(model=="rossler"){
    vector<double> variable(3),parameter(3);
    variable[0] = 2.61622;
    variable[1] = -6.32533;
    variable[2] = 0.0335135;
    parameter[0]= atof(argv[5]);
    parameter[1]= atof(argv[6]);
    parameter[2]= atof(argv[7]);
    
    rossler attractor(variable,parameter,atof(argv[8]),(string) argv[3]);
    lyapunov<Jacobian_rossler>(attractor,pow(10,atoi(argv[9])),pow(10,atoi(argv[10])),1000,(string) argv[4]);
}
if(model=="lorenz"){
  vector<double> variable(3),parameter(3);
  variable[0] = 2.61622;
  variable[1] = -6.32533;
  variable[2] = 3.0335135;
  parameter[0]= atof(argv[5]);
  parameter[1]= atof(argv[6]);
  parameter[2]= atof(argv[7]);
  //cout << parameter[0] << " "<< parameter[1]<< " " <<parameter[2] << endl;
  lorenz attractor(variable,parameter,atof(argv[8]),(string) argv[3]);
  lyapunov<Jacobian_lorenz>(attractor,pow(10,atoi(argv[9])),pow(10,atoi(argv[10])),1000,(string) argv[4]);
  }
  }
  if(measure=="energy"){
    model=argv[2];
    if(model=="Double_Pendulum"){
      string method=argv[3];
      string file_name = argv[4];
      double angle1= atof(argv[5])*M_PI;
      double angle2= atof(argv[6])*M_PI;	
      double dt=atof(argv[7]);
      int steps = atoi(argv[8]);
      
      
      vector<double> variable(4),parameter(5);
      variable[V_THETA1] = angle1;
      variable[V_THETA2] = angle2;
      variable[V_OMEGA1] = 0.0;
      variable[V_OMEGA2] = 0.0;
      parameter[P_L1]= 0.30;
      parameter[P_L2]= 0.30;
      parameter[P_M1]= 0.10;
      parameter[P_M2]= 0.10;
      parameter[P_G]= 9.8;
      
      double_pendulum dp_attractor(variable,parameter,dt,method);
      std::ofstream file;
      std::string Filename = "Energy_"+dp_attractor.get_model_name()+"_"+dp_attractor.get_method_name()+file_name+".out";
      file.open(Filename.c_str());
      double Eref=0.882;
      double E0=dp_attractor.get_energy();
      file << "#Method: =  "<< method << endl;
      file << "#angle1: =  "<< angle1 << endl;
      file << "#angle2: =  "<< angle2 << endl;
      file << "# dt =  "<< dt << endl;
      double error_sum=0;
      for(int i=0;i<steps;i++){
	for(int j=0;j<(1.0/dt);j++){
	  error_sum+=fabs((dp_attractor.get_energy()-E0)/Eref)*dt;                
	  dp_attractor.next();
	}
	file  << dp_attractor.get_t()  << " " << fabs((dp_attractor.get_energy()-E0)/Eref)<< "  " << error_sum << endl;
      }
      
      file.close();
    }
    }
    
    if(measure=="fft"){
      model=argv[2];
      if(model=="Double_Pendulum"){
	string method=argv[3];
	string file_name = argv[4];
	double angle1= atof(argv[5])*M_PI;
	double angle2= atof(argv[6])*M_PI;	
	double dt=atof(argv[7]);
	int steps = atoi(argv[8]);
	
	
	vector<double> variable(4),parameter(5);
	variable[V_THETA1] = angle1;
	variable[V_THETA2] = angle2;
	variable[V_OMEGA1] = 0.0;
	variable[V_OMEGA2] = 0.0;
	parameter[P_L1]= 0.30;
	parameter[P_L2]= 0.30;
	parameter[P_M1]= 0.10;
	parameter[P_M2]= 0.10;
	parameter[P_G]= 9.8;
	
	double_pendulum dp_attractor(variable,parameter,dt,method);
	double data[steps];
	for(int i=0;i<steps;i++){
	  data[i]=dp_attractor[0];	
	  for(int j=0;j<(1.0/dt);j++){
	    dp_attractor.next();
	  }
	}
	FFT(data, steps,file_name);
	
      }
    }
    
    if(measure=="time_series"){
      model=argv[2];
      if(model=="Double_Pendulum"){
	string method=argv[3];
	string file_name = argv[4];
	double angle1= atof(argv[5])*M_PI;
	double angle2= atof(argv[6])*M_PI;	
	double dt=atof(argv[7]);
	int steps = atoi(argv[8]);
	
	
	vector<double> variable(4),parameter(5);
	variable[V_THETA1] = angle1;
	variable[V_THETA2] = angle2;
	variable[V_OMEGA1] = 0.0;
	variable[V_OMEGA2] = 0.0;
	parameter[P_L1]= 0.30;
	parameter[P_L2]= 0.30;
	parameter[P_M1]= 0.10;
	parameter[P_M2]= 0.10;
	parameter[P_G]= 9.8;
	
	double_pendulum dp_attractor(variable,parameter,dt,method);
	std::ofstream file;
	std::string Filename = "time_series_"+dp_attractor.get_model_name()+"_"+dp_attractor.get_method_name()+file_name+".out";
	file.open(Filename.c_str());
	
	file << "#Method: =  "<< method << endl;
	file << "#angle1: =  "<< angle1 << endl;
	file << "#angle2: =  "<< angle2 << endl;
	file << "# dt =  "<< dt << endl;
	double x,x1,y,y1;        	
	for(int i=0;i<steps;i++){
	  for(int j=0;j<(10);j++){	             
	    dp_attractor.next();
	  }
	  x = 2*sin(dp_attractor[0]);
	  y = -2*cos(dp_attractor[0]);
	  x1 = x+2*sin(dp_attractor[1]);
	  y1 = y-2*cos(dp_attractor[1]);
	  file  << x << " " << y << " " << x1 << " " << y1 << endl;
	}
	
	file.close(); 
      }
    }
    if(measure=="GL_pendulum"){
      //double angle1= atof(argv[2])*M_PI;
      //double angle2= atof(argv[3])*M_PI;
      //Pendulum_GL(argc,argv,angle1,angle2);
    }
    
    }
*/    


int main( int argc , char * argv[]) {
    g_argc=argc;
    g_argv=argv;
    vector<double> variable(3),parameter(3);
    variable[0] = 2.61622;
    variable[1] = -6.32533;
    variable[2] = 0.0335135;
    parameter[0]= 0.15;
    parameter[1]= 0.20;
    parameter[2]= 10;
    double t=1;
    rossler_func func;   
    func.set(t,variable,parameter); 
    AdamsBashforth<rossler_func> attractor(variable,parameter,0.001);
//    std::cout << "x " << attractor["y"] << std::endl;

        return 0;
}

    
