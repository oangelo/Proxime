#include <vector>
#include <memory>

#include "numerical_integration/runge_kutta.h"   
#include "functions/rossler.h"   

class MaxLyapunov{
public:
    MaxLyapunov(NumericalIntegration& model, FunctionCapsule& _jacobian, labels_values jacobian_parameters, unsigned transient);
    value operator()(unsigned iterations);
    const value get_exponent() const;
private:
    value LocalLyapunov();

    std::unique_ptr<NumericalIntegration> fiducial;
    std::unique_ptr<FunctionCapsule> jacobian;
    value exponent;
    labels_values base;
    labels_values jacobian_parameters;
};


/*
void GramSchmidt(std::vector<container> & vec_space, labels_values & modulo) {
    std::vector<container> aux_space = vec_space;
    modulo.resize(vec_space.size());
    modulo[0] = sqrt(Dot(aux_space, aux_space));
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
void OrtogonalSpaceNorm(NumericalIntegration & fiducial, labels_values Jparameters,
std::vector<container> & ortogonal_space, container & modulo) {
  
    jacobian_function function(Jparameters);
    std::vector<RungeKutta4Th*> jacobian(fiducial.size());
    
    //Integrats the fiducial trajectory and get the variables for the Jacobian
    ++fiducial;
    labels_values result(fiducial.get_labels_values()); 
    Jparameters.insert(result.begin(), result.end());
    //Integrating the linear equations and saving the vector on ortogonal_space
    for (unsigned i = 0; i < fiducial.size(); i++)
        jacobian[i] = new RungeKutta4Th(function, ortogonal_space[i], fiducial.get_dt());
    for (unsigned i = 0; i < fiducial.size(); i++)
        ++(*jacobian[i]);
    for (unsigned i = 0; i < ortogonal_space.size(); i++)
        ortogonal_space[i] = (*jacobian[i]).get_variable();
    //Ortonormalizing the vectors and saving the modulus
    GramSchmidt(ortogonal_space, modulo);
    //Clear the jacobian, since in the next loop is needed to integrate a new one
    for (unsigned i = 0; i < fiducial.size(); i++)
        delete jacobian[i];
    
}


template <class jacobian_function>
container Lyapunov(NumericalIntegration& fiducial, labels_values parameters, 
        unsigned number_steps, unsigned transients_steps) {
    //Be carfull with vertors of a class, because the vector will copy the object!
    labels_values Jparameters(parameters);
    std::vector<container> ortogonal_space(fiducial.size(), container(fiducial.size()));
    container Lambda(fiducial.size(), 0), modulo(fiducial.size(), 0);

    //Initial Conditions##################################
    for (unsigned i = 0; i < fiducial.size(); i++) {
        for (unsigned j = 0; j < fiducial.size(); j++) {
            if (i == j) {
                ortogonal_space[i][j] = 1;
            } else {
                ortogonal_space[i][j] = 0;
            }
        }
    }
    labels_values result(fiducial.get_labels_values()); 
    Jparameters.insert(result.begin(), result.end());
    //*****************************************************************************
    //********************************Transient Loop*******************************
    //*****************************************************************************
    for (unsigned steps(0); steps < transients_steps; steps++) {
        OrtogonalSpaceNorm<jacobian_function>(fiducial, Jparameters, ortogonal_space, modulo);
    }
    //*****************************************************************************
    //*********************************Data Loop***********************************
    //*****************************************************************************
    for (unsigned steps = 0; steps < number_steps; steps++) {
        OrtogonalSpaceNorm<jacobian_function>(fiducial, Jparameters, ortogonal_space, modulo);
        //calculating the mean
        for (unsigned i = 0; i < fiducial.size(); i++)
            Lambda[i] = (steps * Lambda[i] + log(modulo[i])) / (steps + 1);
    }
    
    //returning the result
    for(unsigned i = 0; i < Lambda.size(); i++)
      Lambda[i] /= fiducial.get_dt();
    return(Lambda);
}
*/


/*
template <class jacobian_function>
value MaxLyapunov(, unsigned number_steps, unsigned transients_steps) {
    labels_values Jparameter(parameters);
    labels_values space_vec;
    value Lambda(0);
    NumericalIntegration* jacobian;
    //Initial Conditions##################################
    labels_values variables(fiducial.get_labels_values());
    labels_values::iterator it(variables.begin());
    space_vec[(*it).first] = 1;
    while(it != variables.end()) {
        space_vec[(*it).first] = 0;
    }
    for (unsigned steps = 0; steps < transients_steps; steps++) {
        ++fiducial;
        labels_values result(fiducial.get_labels_values()); 
        Jparameter.insert(result.begin(), result.end());
        //Integrating the linear equations and saving the vector on ortogonal_space
        jacobian_function function(Jparameter);
        jacobian = new RungeKutta4Th(function, space_vec, fiducial.get_dt());
        ++(*jacobian);
        space_vec = jacobian->get_labels_values();
        //Ortonormalizing the vectors and saving the modulus
        value modulo(sqrt(Dot(space_vec, space_vec)));
        for (labels_values::iterator it(space_vec.begin()); it != space_vec.end(); ++it)
            space_vec[(*it).first] = space_vec[(*it).first] / modulo;
        //Clear the jacobian, since in the next loop is needed to integrate a new one
        delete jacobian;
    }
    for (unsigned steps(0); steps < number_steps; steps++) {
        ++fiducial;
        labels_values result(fiducial.get_labels_values()); 
        Jparameter.insert(result.begin(), result.end());
        jacobian_function function(Jparameter);
        //Integrating the linear equations and saving the vector on ortogonal_space
        jacobian = new RungeKutta4Th(function, space_vec, fiducial.get_dt());
        ++(*jacobian);
        space_vec = jacobian->get_labels_values();
        //Ortonormalizing the vectors and saving the modulus
        value modulo(sqrt(Dot(space_vec, space_vec)));
        for (labels_values::iterator it(space_vec.begin()); it != space_vec.end(); ++it)
            space_vec[(*it).first] = space_vec[(*it).first] / modulo;
        //Clear the jacobian, since in the next loop is needed to integrate a new one
        delete jacobian;
        //Calculating the mean
        Lambda = (steps * Lambda + log(modulo)) / (steps + 1);
    }
    return(Lambda / fiducial.get_dt());
}
*/
