#include "rossler.h"

RosslerFunction::RosslerFunction():
functions_capsule("Rössler System", 3, 
                  dictionary{{"x",0},{"y",1},{"z",2},},
                  dictionary{{"a",0},{"b",1},{"c",2},}),
X(), Y(), Z(), a(), b(), c()
{
}

inline
void RosslerFunction::set(value &t, container & variables, container & parameters){
      
    X=variables[index_var["x"]];
    Y=variables[index_var["y"]];
    Z=variables[index_var["z"]];

    a=parameters[index_par["a"]];
    b=parameters[index_par["b"]];
    c=parameters[index_par["c"]];
  
    result[index_var["x"]]=RosslerFunction::dx();
    result[index_var["y"]]=RosslerFunction::dy();
    result[index_var["z"]]=RosslerFunction::dz();
    
    t=t;
}

inline
value RosslerFunction::dx() {
    value func =  (-Y - Z);
    return(func);
}

inline
value RosslerFunction::dy() {
    value func = X +a* Y;
    return(func);
}

inline
value RosslerFunction::dz() {
    value func = b + Z*(X-c);
    return(func);  
}


Jacobian_RosslerFunction::Jacobian_RosslerFunction():
functions_capsule("Rössler Jacobian", 3, 
                  dictionary{{"x",0},{"y",1},{"z",2},}, 
                  dictionary{{"a",0},{"b",1},{"c",2},{"x",3},{"y",4},{"z",5},}),
X(),Y(),Z(),a(),b(),c(),X_fiducial(),Y_fiducial(),Z_fiducial()
{
    result.clear();
    result.resize(3);
}

inline
void Jacobian_RosslerFunction::set(value &t,container & variables,container & parameters){
    X=variables[0];
    Y=variables[1];
    Z=variables[2];

    a=parameters[0];
    b=parameters[1];
    c=parameters[2];
    X_fiducial=parameters[3];
    Y_fiducial=parameters[4];
    Z_fiducial=parameters[5];

    result[0]=Jacobian_RosslerFunction::dx();
    result[1]=Jacobian_RosslerFunction::dy();
    result[2]=Jacobian_RosslerFunction::dz();

    t=t;
}

inline
value Jacobian_RosslerFunction::dx() {
    value func =  -Y-Z;
    return(func);
}

inline
value Jacobian_RosslerFunction::dy() {
    value func = X+a*Y;
    return(func);
}

inline
value Jacobian_RosslerFunction::dz() {
    value func = +Z_fiducial*X +(X_fiducial-c)*Z;
    return(func);
}
