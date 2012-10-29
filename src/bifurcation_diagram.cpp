/* 
 * File:   bifurcation_diagram.cpp
 * Author: alexandre
 * 
 * Created on 6 de Agosto de 2010, 05:23
 */

#include "bifurcation_diagram.h"


void* thread_func(void* arguments) {


    thread_arg &arg= *((thread_arg*)arguments);
    std::vector< type_container > &zeros=(arg.result);
    zeros = attractor_cross_axi((arg.attractor),
                               arg.time,arg.transiente,
                               arg.coordinate_x,arg.coordinate_y,
                               arg.coordinate_value,arg.quadrant);
    return(NULL);
}

std::vector< type_container > attractor_cross_axi(Numerical_Integration *attractor,int time,int transiente,int coordinate_x,int coordinate_y,type_data coordinate_value,int quadrant)
{
    type_container value((*attractor).size_variable()),next_value((*attractor).size_variable()),aux(2);
    std::vector< type_container > zeros;

    for (int i = 0; i < transiente; i++){
         try {
            (*attractor).next();
        }
        catch (Value_error) {
            std::cout << "Problem with the numerical integration, "
                    "please use a smaller step or a better method" << std::endl;
            zeros.clear();
            return(zeros);
        }
    }
    for (int i=0; i < time; i++){
        try {
            (*attractor).next();
        }
        catch (Value_error) {
            std::cout << "Problem with the numerical integration, "
                    "please use a smaller step or a better method" << std::endl;
            zeros.clear();
            return(zeros);
        }
        value = (*attractor).get_variable();

        try {
            (*attractor).next();
        }
        catch (Value_error) {
            std::cout << "Problem with the numerical integration, "
                    "please use a smaller step or a better method" << std::endl;
            zeros.clear();
            return(zeros);
        }
        next_value = (*attractor).get_variable();

        if (quadrant*value[coordinate_x] > 0) {
                if (((value[coordinate_y] > coordinate_value) && (next_value[coordinate_y] < coordinate_value))
                 || ((value[coordinate_y] < coordinate_value) && (next_value[coordinate_y] > coordinate_value)))
                {
                    aux[0]=(value[coordinate_x]);
                    aux[1]=value[coordinate_y];
                    zeros.push_back(aux);
                }
            }
    }
    return(zeros);
}

