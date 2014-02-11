#include "lorenz.h"

LorenzFunction::LorenzFunction(labels_values parameters):
FunctionCapsule("Lorenz System", 
                  dictionary{{"x",0},{"y",1},{"z",2},},
                  dictionary{{"sigma",0},{"gamma",1},{"beta",2},},
                  parameters),
X(),Y(),Z(), sigma(parameters["sigma"]), gamma(parameters["gamma"]), beta(parameters["beta"])
{} 

LorenzFunction* LorenzFunction::Clone() const{
    return(new LorenzFunction(*this));
}

LorenzFunction* LorenzFunction::Create(labels_values parameters) const{
    return new LorenzFunction(parameters);
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
void LorenzFunction::set(value &t,container & variables){
    X=variables[index_var["x"]];
    Y=variables[index_var["y"]];
    Z=variables[index_var["z"]];

    result[index_var["x"]]=dx();
    result[index_var["y"]]=dy();
    result[index_var["z"]]=dz();
}


Jacobian_LorenzFunction::Jacobian_LorenzFunction(labels_values parameters):
FunctionCapsule("Lorenz Jacobian",  
                  dictionary{{"x",0},{"y",1},{"z",2},},
                  dictionary{{"sigma",0},{"gamma",1},{"beta",2},{"x",3},{"y",4},{"z",5}},
                  parameters),

X(), Y(), Z(), 
sigma(parameters["sigma"]), gamma(parameters["gamma"]), beta(parameters["beta"]), 
X_fiducial(parameters["x"]), Y_fiducial(parameters["y"]), Z_fiducial(parameters["z"])
{
}

Jacobian_LorenzFunction* Jacobian_LorenzFunction::Clone() const{
    return(new Jacobian_LorenzFunction(*this));
}

Jacobian_LorenzFunction* Jacobian_LorenzFunction::Create(labels_values parameters) const{
    return new Jacobian_LorenzFunction(parameters);
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
void Jacobian_LorenzFunction::set(value &t,container & variables){
    X=variables[index_var["x"]];
    Y=variables[index_var["y"]];
    Z=variables[index_var["z"]];
   
    result[index_var["x"]]=dx();
    result[index_var["y"]]=dy();
    result[index_var["z"]]=dz();
}
