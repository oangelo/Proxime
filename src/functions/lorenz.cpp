#include "lorenz.h"

functions_capsule::functions_capsule():
  variable_name_index(), parameter_name_index(), func_name(), result()

{}


LorenzFunction::LorenzFunction():
X(),Y(),Z(), sigma(), gamma(), beta()
{
    //model name
    func_name = "Lorenz System";
    //variables names
//    name_item variable_init[3] = {name_item("x",0), name_item("y",1), name_item("z",2)};
//    name_item parameter_init[3] = {name_item("sigma",0), name_item("gamma",1), name_item("beta",2)};
//    variable_name_index.insert(variable_init, variable_init + 3);
//    parameter_name_index.insert(parameter_init, parameter_init + 3);
    //init internal variables
    result.clear();
    result.resize(3);
}

inline
value LorenzFunction::dx() {
    value func =  (-sigma*X +sigma*Y);
    return(func);
}

inline
value LorenzFunction::dy() {
    value func = (gamma-Z)*X - Y;
    return(func);
}

inline
value LorenzFunction::dz() {
    value func = X*Y-beta*Z;
    return(func);
}


inline
void LorenzFunction::set(value &t,container & variables,container & parameters){
    X=variables[0];
    Y=variables[1];
    Z=variables[2];

    sigma=parameters[0];
    gamma=parameters[1];
    beta=parameters[2];
    
    result[0]=dx();
    result[1]=dy();
    result[2]=dz();

    t=t;
}


Jacobian_LorenzFunction::Jacobian_LorenzFunction():
X(), Y(), Z(), sigma(), gamma(), beta(), X_fiducial(), Y_fiducial(), Z_fiducial()
{
    result.clear();
    result.resize(3);
}

inline
value Jacobian_LorenzFunction::dx() {
    value func = + (-sigma*X +sigma*Y);
    return(func);
}

inline
value Jacobian_LorenzFunction::dy() {
    value func =+(-Z_fiducial+gamma)*X - Y-X_fiducial*Z;
    return(func);
}

inline
value Jacobian_LorenzFunction::dz() {
    value func =  +Y_fiducial*X +X_fiducial*Y -beta*Z;
    return(func);
}


inline
void Jacobian_LorenzFunction::set(value &t,container & variables,container & parameters){

    X=variables[0];
    Y=variables[1];
    Z=variables[2];

    sigma=parameters[0];
    gamma=parameters[1];
    beta=parameters[2];
    X_fiducial=parameters[3];
    Y_fiducial=parameters[4];
    Z_fiducial=parameters[5];
    
    result[0]=dx();
    result[1]=dy();
    result[2]=dz();

    t=t;
}
