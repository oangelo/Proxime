#include "lorenz.h"

functions_capsule::functions_capsule():
  variable_name_index(), parameter_name_index(), func_name(), __result()

{}


type_data functions_capsule::get_result(unsigned i) const {
  return (__result[i]);
  
}

const type_container & functions_capsule::get_result() const {
  return(__result);
}


unsigned functions_capsule::size() const{
      return(__result.size());
}

LorenzFunction::LorenzFunction():
X(),Y(),Z(), sigma(), gamma(), beta()
{
    //model name
    func_name = "Lorenz System";
    //variables names
    name_item variable_init[3] = {name_item("x",0), name_item("y",1), name_item("z",2)};
    name_item parameter_init[3] = {name_item("sigma",0), name_item("gamma",1), name_item("beta",2)};
    variable_name_index.insert(variable_init, variable_init + 3);
    parameter_name_index.insert(parameter_init, parameter_init + 3);
    //init internal variables
    __result.clear();
    __result.resize(3);
}

inline
type_data LorenzFunction::dx() {
    type_data func =  (-sigma*X +sigma*Y);
    return(func);
}

inline
type_data LorenzFunction::dy() {
    type_data func = (gamma-Z)*X - Y;
    return(func);
}

inline
type_data LorenzFunction::dz() {
    type_data func = X*Y-beta*Z;
    return(func);
}


inline
void LorenzFunction::set(type_data &t,type_container & variables,type_container & parameters){
    X=variables[0];
    Y=variables[1];
    Z=variables[2];

    sigma=parameters[0];
    gamma=parameters[1];
    beta=parameters[2];
    
    __result[0]=dx();
    __result[1]=dy();
    __result[2]=dz();

    t=t;
}


Jacobian_LorenzFunction::Jacobian_LorenzFunction():
X(), Y(), Z(), sigma(), gamma(), beta(), X_fiducial(), Y_fiducial(), Z_fiducial()
{
    __result.clear();
    __result.resize(3);
}

inline
type_data Jacobian_LorenzFunction::dx() {
    type_data func = + (-sigma*X +sigma*Y);
    return(func);
}

inline
type_data Jacobian_LorenzFunction::dy() {
    type_data func =+(-Z_fiducial+gamma)*X - Y-X_fiducial*Z;
    return(func);
}

inline
type_data Jacobian_LorenzFunction::dz() {
    type_data func =  +Y_fiducial*X +X_fiducial*Y -beta*Z;
    return(func);
}


inline
void Jacobian_LorenzFunction::set(type_data &t,type_container & variables,type_container & parameters){

    X=variables[0];
    Y=variables[1];
    Z=variables[2];

    sigma=parameters[0];
    gamma=parameters[1];
    beta=parameters[2];
    X_fiducial=parameters[3];
    Y_fiducial=parameters[4];
    Z_fiducial=parameters[5];
    
    __result[0]=dx();
    __result[1]=dy();
    __result[2]=dz();

    t=t;
}
