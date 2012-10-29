#include <vector>

#include "numerical_integration.h"   


type_data Dot(type_container v1, type_container v2) {
    type_data escalar = 0;
    for (unsigned cont = 0; cont < v1.size(); cont++)
        escalar += (v1[cont]) * v2[cont];
    return (escalar);
}

void GramSchmidt(std::vector<type_container> & vec_space, type_container & modulo) {
    std::vector<type_container> aux_space = vec_space;
    modulo.resize(vec_space.size());
    modulo[0] = sqrt(Dot(aux_space[0], aux_space[0]));
    for (unsigned i = 1; i < vec_space.size(); i++) {
        for (unsigned j = 0; j < i; j++)
            for (unsigned cont = 0; cont < vec_space.size(); cont++)
                aux_space[i][cont] += -(Dot(vec_space[i], aux_space[j]) * aux_space[j][cont]) / Dot(aux_space[j], aux_space[j]);
        modulo[i] = sqrt(Dot(aux_space[i], aux_space[i]));
    }
    //normalizing###################################################
    for (unsigned i = 0; i < vec_space.size(); i++)
        for (unsigned j = 0; j < vec_space.size(); j++)
            aux_space[i][j] = aux_space[i][j] / modulo[i];
    //##############################################################
    vec_space = aux_space;
}

template <class jacobian_function>
void ortogonal_space_norm(Numerical_Integration & fiducial,
std::vector<type_container> & ortogonal_space,
type_container & Jparameters,
type_container & modulo) {
  
    RungeKutta<jacobian_function>* jacobian[fiducial.size_variable()];
    
    //Integrating the fiducial trajectory and getting the variables for the Jacobian
    fiducial.next();
    for (unsigned i = 0; i < fiducial.size_variable(); i++)
        Jparameters[fiducial.size_parameter() + i] = fiducial.get_variable(i);
    //Integrating the linear equations and saving the vector on ortogonal_space
    for (unsigned i = 0; i < fiducial.size_variable(); i++)
        jacobian[i] = new RungeKutta<jacobian_function>(ortogonal_space[i], Jparameters, fiducial.get_dt());
    for (unsigned i = 0; i < fiducial.size_variable(); i++)
        (*jacobian[i]).next();
    for (unsigned i = 0; i < ortogonal_space.size(); i++)
        ortogonal_space[i] = (*jacobian[i]).get_variable();
    //Ortonormalizing the vectors and saving the modulus
    GramSchmidt(ortogonal_space, modulo);
    //Clear the jacobian, since in the next loop is needed to integrate a new one
    for (unsigned i = 0; i < fiducial.size_variable(); i++)
        delete jacobian[i];
    
}

template <class jacobian_function>
type_container lyapunov(Numerical_Integration & fiducial, unsigned number_steps, 
        unsigned transients_steps, int steps_to_print,std::string file_name) {
    /*Be carfull with vertors of a class, because the vector will copy the object!*/
    type_container Jparameters(fiducial.size_variable() + fiducial.size_parameter());
    std::vector<type_container> ortogonal_space;
    type_container Lambda(fiducial.size_variable(), 0), modulo(fiducial.size_variable(), 0);
    std::ofstream data_lyapunov;
    std::string Filename = "lyapunov_"+fiducial.get_model_name()+"_"+fiducial.get_method_name()+file_name+".out";
    data_lyapunov.open(Filename.c_str());

    //Initial Conditions##################################
    ortogonal_space.resize(fiducial.size_variable());
    for (unsigned i = 0; i < fiducial.size_variable(); i++) {
        ortogonal_space[i].resize(fiducial.size_variable());
    }
    for (unsigned i = 0; i < fiducial.size_variable(); i++) {
        for (unsigned j = 0; j < fiducial.size_variable(); j++) {
            if (i == j) {
                ortogonal_space[i][j] = 0;
            } else {
                ortogonal_space[i][j] = 1;
            }
        }
    }
    for (unsigned i = 0; i < fiducial.size_parameter(); i++)
        Jparameters[i] = fiducial.get_parameter(i);
    /*****************************************************************************/
    /********************************Transient Loop*******************************/
    /*****************************************************************************/
    for (unsigned steps = 1; steps < transients_steps; steps++) {
        ortogonal_space_norm<jacobian_function > (fiducial, ortogonal_space, Jparameters, modulo);
    }
    /*****************************************************************************/
    /*********************************Data Loop***********************************/
    /*****************************************************************************/
    unsigned cont_print = 0;
    for (unsigned steps = 0; steps < number_steps; steps++) {
        ortogonal_space_norm<jacobian_function > (fiducial,ortogonal_space, Jparameters, modulo);
        //Print the result
        for (unsigned i = 0; i < fiducial.size_variable(); i++)
            Lambda[i] = (steps * Lambda[i] + log(modulo[i])) / (steps + 1);
        cont_print++;
        if (cont_print == (number_steps / steps_to_print)) {
            cont_print = 0;
            data_lyapunov << (steps) * fiducial.get_dt() << " ";
            for (unsigned i = 0; i < fiducial.size_variable(); i++)
                data_lyapunov << Lambda[i] / fiducial.get_dt() << " ";
            data_lyapunov << std::endl;
        }
    }
    data_lyapunov.close();
    
    for(unsigned i=0;i<Lambda.size();i++)
      Lambda[i]=Lambda[i]/fiducial.get_dt();
    return(Lambda);
}

template <class jacobian_function>
type_data lyapunov_max(Numerical_Integration & fiducial, unsigned number_steps, unsigned transients_steps) {
    /*Be carfull with vertors of a class, because the vector will copy the object!*/
    type_container Jparameters(fiducial.size_variable() + fiducial.size_parameter());
    type_container space_vec;
    type_data Lambda;
    type_data modulo;
    jacobian_function* jacobian;
    //Initial Conditions##################################
    modulo=1;
    space_vec.resize(fiducial.size_variable());
    space_vec[0] = 1;
    for (unsigned j = 1; j < space_vec.size(); j++) {
        space_vec[j] = 0;
    }
    for (unsigned i = 0; i < fiducial.size_parameter(); i++)
        Jparameters[i] = fiducial.get_parameter(i);
/*****************************************************************************/
/********************************Transient Loop*******************************/
/*****************************************************************************/
    for (unsigned steps = 0; steps < transients_steps; steps++) {
        for (unsigned i = 0; i < fiducial.size_variable(); i++)
            Jparameters[fiducial.size_parameter() + i] = fiducial.get_variable(i);
        //Integrating the linear equations and saving the vector on ortogonal_space
        jacobian = new jacobian_function(space_vec, Jparameters, fiducial.get_dt());
        (*jacobian).next();
        space_vec = (*jacobian).get_variable();
        //Ortonormalizing the vectors and saving the modulus
        modulo=sqrt(Dot(space_vec,space_vec));
        for (unsigned i = 0; i < fiducial.size_variable(); i++)
            space_vec[i]/=modulo;
        //Clear the jacobian, since in the next loop is needed to integrate a new one
        delete jacobian;
    }
/*****************************************************************************/
/*********************************Data Loop***********************************/
/*****************************************************************************/
    for (unsigned steps = 0; steps < number_steps; steps++) {
        fiducial.next();
        for (unsigned i = 0; i < fiducial.size_variable(); i++)
            Jparameters[fiducial.size_parameter() + i] = fiducial.get_variable(i);
        //Integrating the linear equations and saving the vector on ortogonal_space
        jacobian = new jacobian_function(space_vec, Jparameters, fiducial.get_dt());
        (*jacobian).next();
        space_vec = (*jacobian).get_variable();
        //Ortonormalizing the vectors and saving the modulus
        modulo=sqrt(Dot(space_vec,space_vec));
        for (unsigned i = 0; i < fiducial.size_variable(); i++)
            space_vec[i]/=modulo;
        //Clear the jacobian, since in the next loop is needed to integrate a new one
        delete jacobian;
        //Print the result
        Lambda = (steps * Lambda + log(modulo)) / (steps + 1);
    }
    return(Lambda/fiducial.get_dt());
}

