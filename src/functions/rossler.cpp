#include "rossler.h"

Rossler::Rossler(labels_values parameters):
FunctionCapsule("Rössler System",  
                  dictionary{{"x",0},{"y",1},{"z",2},},
                  dictionary{{"a",0},{"b",1},{"c",2},},
                  parameters),
X(), Y(), Z(), a(parameters["a"]), b(parameters["b"]), c(parameters["c"])
{
}

Rossler* Rossler::Clone() const{
    return(new Rossler(*this));
}

Rossler* Rossler::Create(labels_values parameters) const{
    return new Rossler(parameters);
}

inline
void Rossler::set(value &t, container & variables){
      
    X=variables[index_var["x"]];
    Y=variables[index_var["y"]];
    Z=variables[index_var["z"]];

    result[index_var["x"]]=Rossler::dx();
    result[index_var["y"]]=Rossler::dy();
    result[index_var["z"]]=Rossler::dz();
}

inline
value Rossler::dx() {
    value func =  (-Y - Z);
    return(func);
}

inline
value Rossler::dy() {
    value func = X +a* Y;
    return(func);
}

inline
value Rossler::dz() {
    value func = b + Z*(X-c);
    return(func);  
}


Jacobian_Rossler::Jacobian_Rossler(labels_values parameters):
FunctionCapsule("Rössler Jacobian",  
                  dictionary{{"x",0},{"y",1},{"z",2},}, 
                  dictionary{{"a",0},{"b",1},{"c",2},{"x",3},{"y",4},{"z",5},},
                  parameters),
X(),Y(),Z(),
a(parameters["a"]),b(parameters["b"]),c(parameters["c"]),
X_fiducial(parameters["x"]),Y_fiducial(parameters["y"]),Z_fiducial(parameters["z"])
{}

inline
Jacobian_Rossler* Jacobian_Rossler::Clone() const{
    return(new Jacobian_Rossler(*this));
}

inline
Jacobian_Rossler* Jacobian_Rossler::Create(labels_values parameters) const{
    return new Jacobian_Rossler(parameters);
}

inline
void Jacobian_Rossler::set(value &t,container & variables){
    X=variables[index_var["x"]];
    Y=variables[index_var["y"]];
    Z=variables[index_var["z"]];

    result[index_var["x"]]=Jacobian_Rossler::dx();
    result[index_var["y"]]=Jacobian_Rossler::dy();
    result[index_var["z"]]=Jacobian_Rossler::dz();
}

inline
value Jacobian_Rossler::dx() {
    value func =  -Y-Z;
    return(func);
}

inline
value Jacobian_Rossler::dy() {
    value func = X+a*Y;
    return(func);
}

inline
value Jacobian_Rossler::dz() {
    value func = +Z_fiducial*X +(X_fiducial-c)*Z;
    return(func);
}

