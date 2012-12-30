#include "rossler.h"

RosslerFunction::RosslerFunction():
X(), Y(), Z(), a(), b(), c()
{
    //model name
    func_name="Rössler System";
    //variables names
    name_item variable_init[3] = {name_item("x",0), name_item("y",1), name_item("z",2)};
    name_item parameter_init[3] = {name_item("a",0), name_item("b",1), name_item("c",2)};
    variable_name_index.insert(variable_init, variable_init + 3);
    parameter_name_index.insert(parameter_init, parameter_init + 3);
    //init internal variables
    __result.clear();
    __result.resize(3);
}

inline
void RosslerFunction::set(type_data &t, type_container & variables, type_container & parameters){
      
    X=variables[0];
    Y=variables[1];
    Z=variables[2];

    a=parameters[0];
    b=parameters[1];
    c=parameters[2];
  
    __result[0]=RosslerFunction::dx();
    __result[1]=RosslerFunction::dy();
    __result[2]=RosslerFunction::dz();
    
    t=t;
}

inline
type_data RosslerFunction::dx() {
    type_data func =  (-Y - Z);
    return(func);
}

inline
type_data RosslerFunction::dy() {
    type_data func = X +a* Y;
    return(func);
}

inline
type_data RosslerFunction::dz() {
    type_data func = b + Z*(X-c);
    return(func);  
}


Jacobian_RosslerFunction::Jacobian_RosslerFunction():
X(),Y(),Z(),a(),b(),c(),X_fiducial(),Y_fiducial(),Z_fiducial()
{
    __result.clear();
    __result.resize(3);
}

inline
void Jacobian_RosslerFunction::set(type_data &t,type_container & variables,type_container & parameters){
    X=variables[0];
    Y=variables[1];
    Z=variables[2];

    a=parameters[0];
    b=parameters[1];
    c=parameters[2];
    X_fiducial=parameters[3];
    Y_fiducial=parameters[4];
    Z_fiducial=parameters[5];

    __result[0]=Jacobian_RosslerFunction::dx();
    __result[1]=Jacobian_RosslerFunction::dy();
    __result[2]=Jacobian_RosslerFunction::dz();

    t=t;
}

inline
type_data Jacobian_RosslerFunction::dx() {
    type_data func =  -Y-Z;
    return(func);
}

inline
type_data Jacobian_RosslerFunction::dy() {
    type_data func = X+a*Y;
    return(func);
}

inline
type_data Jacobian_RosslerFunction::dz() {
    type_data func = +Z_fiducial*X +(X_fiducial-c)*Z;
    return(func);
}