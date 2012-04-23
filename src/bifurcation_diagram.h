#include <pthread.h>
#include "Numerical_Integration.h"
#include <fstream>
#include <mpi.h>

class thread_arg{
public:
    Numerical_Integration *attractor;
    std::vector< std::vector<double> > result;
    double coordinate_value;
    int coordinate_x, coordinate_y;
    int time, transiente,quadrant;
};

void* thread_func(void* arguments);


/*Be carful with vertors of a class, because the vector will copy the object!
 * Is better to use arrays that point to those objects!
 */

std::vector< std::vector<double> > attractor_cross_axi(Numerical_Integration *attractor,int time,int transiente,int coordinate_x,int coordinate_y,double coordinate_value,int quadrant);

//*
template<class function>
void bifurcation(std::vector<double> &variable,std::vector<double> &parameter,double dt,
                 unsigned parameter_index,double coordinate_value,int quadrant,
                 int coordinate_x,int coordinate_y,
                 double init,double end,unsigned n_points,
                 unsigned time,unsigned transiente,unsigned number_threads){

    std::vector<pthread_t> threads(number_threads);
    std::vector<Numerical_Integration*> attractors(n_points);
    //Generating the atractos objects//
    unsigned cont=0;
    for (double c_parameter = init; c_parameter < end; c_parameter+=((double)(end-init))/n_points) {
        parameter[parameter_index]=c_parameter;
        //memory leak!!!!!!!!!!!!!!!!!!!!!!!!!!!
	attractors[cont]=(new AdamsBashforth<function>(variable,parameter,dt));
	cont++;
    }
    //****************************Generating the threads arguments**********************************
    std::vector< thread_arg> args(attractors.size());
    for (unsigned count = 0; count < args.size(); count++){
        args[count].coordinate_value=coordinate_value;
        args[count].coordinate_x=coordinate_x;
        args[count].coordinate_y=coordinate_y;
        args[count].time=time;
        args[count].transiente=transiente;
        args[count].quadrant=quadrant;
    }
    //**********************Creating threads and starting the calculations**************************
    std::cout << "Starting calculation of "
         << (*attractors[0]).get_model_name()
         << " Bifurcation diagram:" << std::endl;
    
    for (unsigned count = 0; count < attractors.size()-number_threads+1 ; count+=number_threads) {
        for (unsigned c_threads = 0; c_threads < number_threads; c_threads++){
            args[count+c_threads].attractor=attractors[count+c_threads];
            pthread_create(&(threads[c_threads]), NULL, thread_func, (void*) &(args[count+c_threads]));
        }
        std::cout << "Steps to complete:" << attractors.size()-count << std::endl;
        for (unsigned c_threads = 0; c_threads < number_threads; c_threads++)
            pthread_join( threads[c_threads], NULL);
    }
    //******************************Writing data to file********************************************
    std::ofstream bifurcation;
    std::string Filename = "data_out/bifurcation_"+(*attractors[0]).get_model_name()+".out";
    bifurcation.open(Filename.c_str());
    for (unsigned j = 0; j < args.size(); j++)
        for (unsigned i = 0; i < args[j].result.size(); i++)
            bifurcation << (*args[j].attractor).get_parameter(parameter_index)<< "\t" << args[j].result[i][0] << "\t"<< args[j].result[i][1] <<  std::endl;
    bifurcation.close();
    //*************************free the alocated attractors objects*********************************
    for (unsigned count = 0; count < attractors.size(); count++)
        delete attractors[count];
}
//*/
/*
template <class function>
void portrail(Numerical_Integration & attractor,
              int time, int transiente,
              int coordinate_x,int coordinate_y,
              double coordinate_value,int quadrant)
{
    std::cout << "Generating phase portrail of " << attractor.get_model_name() << std::endl;
    std::vector< std::vector<double> > zeros;
    AdamsBashforth<function> attractor_aux(attractor);
    zeros = attractor_cross_axi(attractor_aux,time/2,transiente,coordinate_x,coordinate_y,coordinate_value,quadrant);
    std::cout << "Number of zeros = " << zeros.size() << std::endl;

    std::ofstream portrail;
    std::string Filename = "data_out/portrail_"+attractor.get_model_name()+".out";
    portrail.open(Filename.c_str());

    for(int t=0;t<transiente;t++)
        attractor.next();
    for(int t=0;t<time;t++){
        portrail << attractor.get_variable(coordinate_x)
                 <<"\t"
                 << attractor.get_variable(coordinate_y)
                 << std::endl;
        attractor.next();
    }
    portrail << "   "<< std::endl;
    portrail << "   "<< std::endl;

    for(int i=0;i<zeros.size();i++){
        portrail << zeros[i][0] << "\t" << zeros[i][1] << std::endl;
    }

    portrail.close();
}
//*/
//use calback functions to alter the parameter!
template<class function>
void MPI_BIFURCATIONS(std::vector<double> &variable,std::vector<double> &parameter,double dt,
                 double coordinate_value,int quadrant,
                 int coordinate_x,int coordinate_y,
                 double init,double end,int n_points,
                 int time,int transiente,void (*increment_func)(std::vector<double> &variable,std::vector<double> &parameter,double increment),int argc , char * argv[]){

    
   int total_steps=n_points;
  
   /***********************************MPI VARS****************************************************/
   int MPI_MTAG1=1,MPI_MTAG2=2,MPI_MTAG3=3,MPI_MTAG4=4;
   int MPI_MYID , MPI_NUMPROCS,MPI_ISLAVE,MPI_SLAVE_STEPS;
   /***********************************************************************************************/
   /************************************MPI INIT***************************************************/
   MPI_Status MPI_STATUS ;
   MPI_Init (&argc ,&argv );
   MPI_Comm_size ( MPI_COMM_WORLD ,&MPI_NUMPROCS);
   MPI_Comm_rank ( MPI_COMM_WORLD ,&MPI_MYID);
   /************************************************************************************************/
   MPI_SLAVE_STEPS=total_steps/(MPI_NUMPROCS-1);
   //std::cout<< "slave steps = " << MPI_SLAVE_STEPS << std::endl;
  /*######################################The Slave program########################################*/ 
  std::cout << "Slave " << MPI_MYID << std::endl;
  if(MPI_MYID!=0){
    /*Generating the atractos objects*/
    std::vector<Numerical_Integration*> attractors(MPI_SLAVE_STEPS);
    double MY_PARAMETER_INTERVAL=end-init;
    double MY_PARAMETER_STEP=(MY_PARAMETER_INTERVAL/(MPI_NUMPROCS-1))/MPI_SLAVE_STEPS;
    //double MY_PARAMETER_END=init+(MPI_MYID)*(MY_PARAMETER_INTERVAL/(MPI_NUMPROCS-1));
    double MY_PARAMETER_BEGIN=init+(MPI_MYID-1)*(MY_PARAMETER_INTERVAL/(MPI_NUMPROCS-1));

    std::vector<double> x_zeros;
    int size_vec_x;
    std::vector<double> my_parameters(MPI_SLAVE_STEPS);
    std::vector<int> size(MPI_SLAVE_STEPS);
    
    double count_parameter=MY_PARAMETER_BEGIN;
    for (int count=0; count<MPI_SLAVE_STEPS;count++) {
      increment_func(variable,parameter,count_parameter);
      attractors[count]=(new AdamsBashforth<function>(variable,parameter,dt));
      my_parameters[count]=count_parameter;
      count_parameter+=MY_PARAMETER_STEP;
    }

     //std::cout << "ID = " <<MPI_MYID << std::endl;
    for (int step = 0; step <MPI_SLAVE_STEPS; step++){
      std::vector< std::vector<double> > zeros;    
      zeros = attractor_cross_axi((attractors[step]),time,transiente,coordinate_x,coordinate_y,coordinate_value,quadrant);
      
      size[step]=zeros.size();
      for(unsigned int i=0;i<zeros.size();i++){
        x_zeros.push_back(zeros[i][0]);
	
      }
      
    }
    for (unsigned int count = 0; count < attractors.size(); count++)
    	delete attractors[count];
    
    size_vec_x=x_zeros.size();	
    MPI_Send(&size_vec_x,1, MPI_INTEGER , 0, MPI_MTAG1, MPI_COMM_WORLD);
    MPI_Send(&x_zeros[0],size_vec_x, MPI_DOUBLE_PRECISION , 0, MPI_MTAG2, MPI_COMM_WORLD);
    MPI_Send(&size[0],MPI_SLAVE_STEPS, MPI_INTEGER , 0, MPI_MTAG3, MPI_COMM_WORLD);
    MPI_Send(&my_parameters[0],MPI_SLAVE_STEPS, MPI_DOUBLE_PRECISION , 0, MPI_MTAG4, MPI_COMM_WORLD); 
    }else{
    /*################################The Master program############################################*/
    	std::cout << "Master" << std::endl;
        std::vector<double> x_zeros;
	x_zeros.reserve(1000);
	int size_vec_x;
    	std::vector<double> my_parameters(MPI_SLAVE_STEPS);
    	std::vector<int> size(MPI_SLAVE_STEPS);
	
	std::ofstream data_out;
	std::string Filename = "Bifurcation.out";
	data_out.open(Filename.c_str());
	std::cout << "file_open" << std::endl;
  	
	for ( MPI_ISLAVE =1; MPI_ISLAVE < MPI_NUMPROCS ; MPI_ISLAVE++){
		MPI_Recv(&size_vec_x,1, MPI_INTEGER, MPI_ISLAVE , MPI_MTAG1, MPI_COMM_WORLD,&MPI_STATUS);
		x_zeros.resize(size_vec_x);		
		MPI_Recv(&x_zeros[0],size_vec_x, MPI_DOUBLE_PRECISION, MPI_ISLAVE , MPI_MTAG2, MPI_COMM_WORLD,&MPI_STATUS);
		MPI_Recv(&size[0],MPI_SLAVE_STEPS, MPI_INTEGER, MPI_ISLAVE , MPI_MTAG3, MPI_COMM_WORLD,&MPI_STATUS);
		MPI_Recv(&my_parameters[0],MPI_SLAVE_STEPS, MPI_DOUBLE_PRECISION, MPI_ISLAVE , MPI_MTAG4, MPI_COMM_WORLD,&MPI_STATUS);
		
			int count=0;
			for (int step = count; step < MPI_SLAVE_STEPS; step++){
				for (int i = count; i < count+size[step]; i++){
					data_out << my_parameters[step] << " ";
					data_out << x_zeros[i] << std::endl;
				}
				count+=size[step];
			}
				
	}
     data_out.close();
    }
    MPI_Finalize();
}

