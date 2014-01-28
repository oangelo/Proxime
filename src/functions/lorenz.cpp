#include "lorenz.h"

LorenzFunction::LorenzFunction():
FunctionCapsule("Lorenz System", 3, 
                  dictionary{{"x",0},{"y",1},{"z",2},},
                  dictionary{{"sigma",0},{"gamma",1},{"beta",2},}),
X(),Y(),Z(), sigma(), gamma(), beta()
{} 

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
FunctionCapsule("Lorenz Jacobian", 3, 
                  dictionary{{"x",0},{"y",1},{"z",2},},
                  dictionary{{"sigma",0},{"gamma",1},{"beta",2},{"x",3},{"y",4},{"z",5}}),

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
