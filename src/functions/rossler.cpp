#include "rossler.h"

RosslerFunction::RosslerFunction(labels_values parameters):
FunctionCapsule("Rössler System",  
                  dictionary{{"x",0},{"y",1},{"z",2},},
                  dictionary{{"a",0},{"b",1},{"c",2},},
                  parameters),
X(), Y(), Z(), a(parameters["a"]), b(parameters["b"]), c(parameters["c"])
{
}

RosslerFunction* RosslerFunction::Clone() const{
    return(new RosslerFunction(*this));
}

RosslerFunction* RosslerFunction::Create(labels_values parameters) const{
    return new RosslerFunction(parameters);
}

inline
void RosslerFunction::set(value &t, container & variables){
      
    X=variables[index_var["x"]];
    Y=variables[index_var["y"]];
    Z=variables[index_var["z"]];

    result[index_var["x"]]=RosslerFunction::dx();
    result[index_var["y"]]=RosslerFunction::dy();
    result[index_var["z"]]=RosslerFunction::dz();
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


Jacobian_RosslerFunction::Jacobian_RosslerFunction(labels_values parameters):
FunctionCapsule("Rössler Jacobian",  
                  dictionary{{"x",0},{"y",1},{"z",2},}, 
                  dictionary{{"a",0},{"b",1},{"c",2},{"x",3},{"y",4},{"z",5},},
                  parameters),
X(),Y(),Z(),
a(parameters["a"]),b(parameters["b"]),c(parameters["c"]),
X_fiducial(parameters["x"]),Y_fiducial(parameters["y"]),Z_fiducial(parameters["z"])
{}

inline
Jacobian_RosslerFunction* Jacobian_RosslerFunction::Clone() const{
    return(new Jacobian_RosslerFunction(*this));
}

inline
Jacobian_RosslerFunction* Jacobian_RosslerFunction::Create(labels_values parameters) const{
    return new Jacobian_RosslerFunction(parameters);
}

inline
void Jacobian_RosslerFunction::set(value &t,container & variables){
    X=variables[index_var["x"]];
    Y=variables[index_var["y"]];
    Z=variables[index_var["z"]];

    result[index_var["x"]]=Jacobian_RosslerFunction::dx();
    result[index_var["y"]]=Jacobian_RosslerFunction::dy();
    result[index_var["z"]]=Jacobian_RosslerFunction::dz();
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

