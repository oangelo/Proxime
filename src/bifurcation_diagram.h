#include <pthread.h>
#include "numerical_integration.h"
#include <fstream>
#include <mpi.h>

class thread_arg{
public:
    Numerical_Integration *attractor;
    std::vector< type_container > result;
    type_data coordinate_value;
    int coordinate_x, coordinate_y;
    int time, transiente,quadrant;
};

void* thread_func(void* arguments);


/*Be carful with vertors of a class, because the vector will copy the object!
 * Is better to use arrays that point to those objects!
 */

std::vector< type_container > attractor_cross_axi(Numerical_Integration *attractor,int time,int transiente,int coordinate_x,int coordinate_y,type_data coordinate_value,int quadrant);

//*
template<class function>
void bifurcation(type_container &variable,type_container &parameter,type_data dt,
                 unsigned parameter_index,type_data coordinate_value,int quadrant,
                 int coordinate_x,int coordinate_y,
                 type_data init,type_data end,unsigned n_points,
                 unsigned time,unsigned transiente,unsigned number_threads){

    std::vector<pthread_t> threads(number_threads);
    std::vector<Numerical_Integration*> attractors(n_points);
    //Generating the atractos objects//
    unsigned cont=0;
    for (type_data c_parameter = init; c_parameter < end; c_parameter+=((type_data)(end-init))/n_points) {
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
              type_data coordinate_value,int quadrant)
{
    std::cout << "Generating phase portrail of " << attractor.get_model_name() << std::endl;
    std::vector< type_container > zeros;
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

