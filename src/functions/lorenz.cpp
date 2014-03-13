#include "lorenz.h"

Lorenz::Lorenz(labels_values parameters):
FunctionCapsule("Lorenz System", 
                  dictionary{{"x",0},{"y",1},{"z",2},},
                  dictionary{{"sigma",0},{"gamma",1},{"beta",2},},
                  parameters),
X(),Y(),Z(), sigma(parameters["sigma"]), gamma(parameters["gamma"]), beta(parameters["beta"])
{} 

Lorenz* Lorenz::Clone() const{
    return(new Lorenz(*this));
}

Lorenz* Lorenz::Create(labels_values parameters) const{
    return new Lorenz(parameters);
}

inline
value Lorenz::dx() {
    value func =  (-sigma*X +sigma*Y);
    return(func);
}

inline
value Lorenz::dy() {
    value func = (gamma-Z)*X - Y;
    return(func);
}

inline
value Lorenz::dz() {
    value func = X*Y-beta*Z;
    return(func);
}


inline
void Lorenz::set(value &t,container & variables){
    X=variables[index_var["x"]];
    Y=variables[index_var["y"]];
    Z=variables[index_var["z"]];

    result[index_var["x"]]=dx();
    result[index_var["y"]]=dy();
    result[index_var["z"]]=dz();
}


Jacobian_Lorenz::Jacobian_Lorenz(labels_values parameters):
FunctionCapsule("Lorenz Jacobian",  
                  dictionary{{"x",0},{"y",1},{"z",2},},
                  dictionary{{"sigma",0},{"gamma",1},{"beta",2},{"x",3},{"y",4},{"z",5}},
                  parameters),

X(), Y(), Z(), 
sigma(parameters["sigma"]), gamma(parameters["gamma"]), beta(parameters["beta"]), 
X_fiducial(parameters["x"]), Y_fiducial(parameters["y"]), Z_fiducial(parameters["z"])
{
}

Jacobian_Lorenz* Jacobian_Lorenz::Clone() const{
    return(new Jacobian_Lorenz(*this));
}

Jacobian_Lorenz* Jacobian_Lorenz::Create(labels_values parameters) const{
    return new Jacobian_Lorenz(parameters);
}

inline
value Jacobian_Lorenz::dx() {
    value func = + (-sigma*X +sigma*Y);
    return(func);
}

inline
value Jacobian_Lorenz::dy() {
    value func =+(-Z_fiducial+gamma)*X - Y-X_fiducial*Z;
    return(func);
}

inline
value Jacobian_Lorenz::dz() {
    value func =  +Y_fiducial*X +X_fiducial*Y -beta*Z;
    return(func);
}


inline
void Jacobian_Lorenz::set(value &t,container & variables){
    X=variables[index_var["x"]];
    Y=variables[index_var["y"]];
    Z=variables[index_var["z"]];
   
    result[index_var["x"]]=dx();
    result[index_var["y"]]=dy();
    result[index_var["z"]]=dz();
}
