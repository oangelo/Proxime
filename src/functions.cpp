#include "functions.h"

const double functions_capsule::get_result(unsigned i) const {
  return (__result[i]);
  
}
const std::vector<double> & functions_capsule::get_result() const {
  return(__result);
}
const unsigned functions_capsule::size() const{
      return(__result.size());
}
/******************************************************************************/
/********************************Rosler funcs**********************************/
/******************************************************************************/
void rossler_func::set(double &t,std::vector<double> & variables,std::vector<double> & parameters){
      
    X=variables[0];
    Y=variables[1];
    Z=variables[2];

    a=parameters[0];
    b=parameters[1];
    c=parameters[2];
  
    __result[0]=rossler_func::dx();
    __result[1]=rossler_func::dy();
    __result[2]=rossler_func::dz();
    
}

double rossler_func::dx() {
    double func =  (-Y - Z);
    return(func);
}

double rossler_func::dy() {
    double func = X +a* Y;
    return(func);
}

double rossler_func::dz() {
    double func = b + Z*(X-c);
    return(func);  
}


void Jacobian_rossler_func::set(double &t,std::vector<double> & variables,std::vector<double> & parameters){
    X=variables[0];
    Y=variables[1];
    Z=variables[2];

    a=parameters[0];
    b=parameters[1];
    c=parameters[2];
    X_fiducial=parameters[3];
    Y_fiducial=parameters[4];
    Z_fiducial=parameters[5];

    __result[0]=Jacobian_rossler_func::dx();
    __result[1]=Jacobian_rossler_func::dy();
    __result[2]=Jacobian_rossler_func::dz();

}

double Jacobian_rossler_func::dx() {
    double func =  -Y-Z;
    return(func);
}

double Jacobian_rossler_func::dy() {
    double func = X+a*Y;
    return(func);
}

double Jacobian_rossler_func::dz() {
    double func = +Z_fiducial*X +(X_fiducial-c)*Z;
    return(func);
}

/******************************************************************************/
/********************************lorenz funcs**********************************/
/******************************************************************************/

double lorenz_func::dx() {
    double func =  (-sigma*X +sigma*Y);
    return(func);
}

double lorenz_func::dy() {
    double func = (gamma-Z)*X - Y;
    return(func);
}

double lorenz_func::dz() {
    double func = X*Y-beta*Z;
    return(func);
}


void lorenz_func::set(double &t,std::vector<double> & variables,std::vector<double> & parameters){
    X=variables[0];
    Y=variables[1];
    Z=variables[2];

    sigma=parameters[0];
    gamma=parameters[1];
    beta=parameters[2];
    
    __result[0]=dx();
    __result[1]=dy();
    __result[2]=dz();

}


double Jacobian_lorenz_func::dx() {
    double func = + (-sigma*X +sigma*Y);
    return(func);
}

double Jacobian_lorenz_func::dy() {
    double func =+(-Z_fiducial+gamma)*X - Y-X_fiducial*Z;
    return(func);
}

double Jacobian_lorenz_func::dz() {
    double func =  +Y_fiducial*X +X_fiducial*Y -beta*Z;
    return(func);
}


void Jacobian_lorenz_func::set(double &t,std::vector<double> & variables,std::vector<double> & parameters){

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

}

/******************************************************************************/
/*******************************Double Pendulum funcs**************************/
/******************************************************************************/

void double_pendulum_func::set(double& t, std::vector<double>& variables, std::vector<double>& parameters){ 
    theta1 = variables[0];
    theta2 = variables[1];
    omega1 = variables[2];
    omega2 = variables[3];
    
    l1 = parameters[0];
    l2 = parameters[1];
    m1 = parameters[2];
    m2 = parameters[3];
    g = parameters[4];
    
    __result[0]=dTheta1();
    __result[1]=dTheta2();
    __result[2]=dOmega1();
    __result[3]=dOmega2();
   
}
double double_pendulum_func::dTheta1() {
    double func;
    func = omega1;
    return (func);
}

double double_pendulum_func::dTheta2() {
    double func;
    func = omega2;
    return (func);
}

double double_pendulum_func::dOmega1() {
    double func;
    func = -(-m2 * g * sin(theta1) - g * sin(theta1) * m1 - m2 * l2 * pow(omega2, 0.2e1) * sin(theta1 - theta2) + m2 * cos(theta1 - theta2) * g * sin(theta2) - m2 * cos(theta1 - theta2) * l1 * pow(omega1, 0.2e1) * sin(theta1 - theta2)) / l1 / (-m2 - m1 + m2 * pow(cos(theta1 - theta2), 0.2e1));
    return (func);
}

double double_pendulum_func::dOmega2() {
    double func;
    func =  -(m2 * cos(theta1 - theta2) * g * sin(theta1) + cos(theta1 - theta2) * g * sin(theta1) * m1 + cos(theta1 - theta2) * m2 * l2 * pow(omega2, 0.2e1) * sin(theta1 - theta2) - m2 * g * sin(theta2) - g * sin(theta2) * m1 + l1 * pow(omega1, 0.2e1) * sin(theta1 - theta2) * m2 + l1 * pow(omega1, 0.2e1) * sin(theta1 - theta2) * m1) / l2 / (-m2 - m1 + m2 * pow(cos(theta1 - theta2), 0.2e1));
    return (func);
}

void jacobian_double_pendulum_func::set(double& t, std::vector<double>& variables, std::vector<double>& parameters){

    theta1 = variables[0];
    theta2 = variables[1];
    omega1 = variables[2];
    omega2 = variables[3];
      
    
    l1 = parameters[0];
    l2 = parameters[1];
    m1 = parameters[2];
    m2 = parameters[3];
    g = parameters[4];
    double _theta1 = parameters[5];
    double _theta2 = parameters[6];
    double _omega1 = parameters[7];
    double _omega2 = parameters[8];

    
    Matrix_Jacob(_theta1,_theta2,_omega1,_omega2);

    __result[0]=JdTheta1();
    __result[1]=JdTheta2();
    __result[2]=JdOmega1();
    __result[3]=JdOmega2();
}

void jacobian_double_pendulum_func::Matrix_Jacob(double theta1, double theta2,double omega1, double omega2)
{

    Jacobian[0][0] = -(-m2 * g * cos(theta1) - g * cos(theta1) * m1 - cos(theta1 - theta2) * m2 * l2 * pow(omega2, 0.2e1) - m2 * sin(theta1 - theta2) * g * sin(theta2) + m2 * pow(sin(theta1 - theta2), 0.2e1) * l1 * pow(omega1, 0.2e1) - m2 * pow(cos(theta1 - theta2), 0.2e1) * l1 * pow(omega1, 0.2e1)) / l1 / (-m2 - m1 + m2 * pow(cos(theta1 - theta2), 0.2e1)) - 0.2e1 * (-m2 * g * sin(theta1) - g * sin(theta1) * m1 - m2 * l2 * pow(omega2, 0.2e1) * sin(theta1 - theta2) + m2 * cos(theta1 - theta2) * g * sin(theta2) - m2 * cos(theta1 - theta2) * l1 * pow(omega1, 0.2e1) * sin(theta1 - theta2)) / l1 * pow(-m2 - m1 + m2 * pow(cos(theta1 - theta2), 0.2e1), -0.2e1) * m2 * cos(theta1 - theta2) * sin(theta1 - theta2);
    Jacobian[0][1] = -(cos(theta1 - theta2) * m2 * l2 * pow(omega2, 0.2e1) + m2 * sin(theta1 - theta2) * g * sin(theta2) + m2 * cos(theta1 - theta2) * g * cos(theta2) - m2 * pow(sin(theta1 - theta2), 0.2e1) * l1 * pow(omega1, 0.2e1) + m2 * pow(cos(theta1 - theta2), 0.2e1) * l1 * pow(omega1, 0.2e1)) / l1 / (-m2 - m1 + m2 * pow(cos(theta1 - theta2), 0.2e1)) + 0.2e1 * (-m2 * g * sin(theta1) - g * sin(theta1) * m1 - m2 * l2 * pow(omega2, 0.2e1) * sin(theta1 - theta2) + m2 * cos(theta1 - theta2) * g * sin(theta2) - m2 * cos(theta1 - theta2) * l1 * pow(omega1, 0.2e1) * sin(theta1 - theta2)) / l1 * pow(-m2 - m1 + m2 * pow(cos(theta1 - theta2), 0.2e1), -0.2e1) * m2 * cos(theta1 - theta2) * sin(theta1 - theta2);
    Jacobian[0][2] = 0.2e1 * m2 * cos(theta1 - theta2) * omega1 * sin(theta1 - theta2) / (-m2 - m1 + m2 * pow(cos(theta1 - theta2), 0.2e1));
    Jacobian[0][3] = 0.2e1 * m2 * l2 * omega2 * sin(theta1 - theta2) / l1 / (-m2 - m1 + m2 * pow(cos(theta1 - theta2), 0.2e1));
    Jacobian[1][0] = -(-m2 * sin(theta1 - theta2) * g * sin(theta1) + m2 * cos(theta1 - theta2) * g * cos(theta1) - sin(theta1 - theta2) * g * sin(theta1) * m1 + cos(theta1 - theta2) * g * cos(theta1) * m1 - pow(sin(theta1 - theta2), 0.2e1) * m2 * l2 * pow(omega2, 0.2e1) + pow(cos(theta1 - theta2), 0.2e1) * m2 * l2 * pow(omega2, 0.2e1) + m2 * cos(theta1 - theta2) * l1 * pow(omega1, 0.2e1) + l1 * pow(omega1, 0.2e1) * cos(theta1 - theta2) * m1) / l2 / (-m2 - m1 + m2 * pow(cos(theta1 - theta2), 0.2e1)) - 0.2e1 * (m2 * cos(theta1 - theta2) * g * sin(theta1) + cos(theta1 - theta2) * g * sin(theta1) * m1 + cos(theta1 - theta2) * m2 * l2 * pow(omega2, 0.2e1) * sin(theta1 - theta2) - m2 * g * sin(theta2) - g * sin(theta2) * m1 + l1 * pow(omega1, 0.2e1) * sin(theta1 - theta2) * m2 + l1 * pow(omega1, 0.2e1) * sin(theta1 - theta2) * m1) / l2 * pow(-m2 - m1 + m2 * pow(cos(theta1 - theta2), 0.2e1), -0.2e1) * m2 * cos(theta1 - theta2) * sin(theta1 - theta2);
    Jacobian[1][1] = -(m2 * sin(theta1 - theta2) * g * sin(theta1) + sin(theta1 - theta2) * g * sin(theta1) * m1 + pow(sin(theta1 - theta2), 0.2e1) * m2 * l2 * pow(omega2, 0.2e1) - pow(cos(theta1 - theta2), 0.2e1) * m2 * l2 * pow(omega2, 0.2e1) - m2 * g * cos(theta2) - g * cos(theta2) * m1 - m2 * cos(theta1 - theta2) * l1 * pow(omega1, 0.2e1) - l1 * pow(omega1, 0.2e1) * cos(theta1 - theta2) * m1) / l2 / (-m2 - m1 + m2 * pow(cos(theta1 - theta2), 0.2e1)) + 0.2e1 * (m2 * cos(theta1 - theta2) * g * sin(theta1) + cos(theta1 - theta2) * g * sin(theta1) * m1 + cos(theta1 - theta2) * m2 * l2 * pow(omega2, 0.2e1) * sin(theta1 - theta2) - m2 * g * sin(theta2) - g * sin(theta2) * m1 + l1 * pow(omega1, 0.2e1) * sin(theta1 - theta2) * m2 + l1 * pow(omega1, 0.2e1) * sin(theta1 - theta2) * m1) / l2 * pow(-m2 - m1 + m2 * pow(cos(theta1 - theta2), 0.2e1), -0.2e1) * m2 * cos(theta1 - theta2) * sin(theta1 - theta2);
    Jacobian[1][2] = -(0.2e1 * m2 * l1 * omega1 * sin(theta1 - theta2) + 0.2e1 * l1 * omega1 * sin(theta1 - theta2) * m1) / l2 / (-m2 - m1 + m2 * pow(cos(theta1 - theta2), 0.2e1));
    Jacobian[1][3] = -0.2e1 * m2 * cos(theta1 - theta2) * omega2 * sin(theta1 - theta2) / (-m2 - m1 + m2 * pow(cos(theta1 - theta2), 0.2e1));
    Jacobian[2][0] = 0;
    Jacobian[2][1] = 0;
    Jacobian[2][2] = 1;
    Jacobian[2][3] = 0;
    Jacobian[3][0] = 0;
    Jacobian[3][1] = 0;
    Jacobian[3][2] = 0;
    Jacobian[3][3] = 1;


}

double jacobian_double_pendulum_func::JdTheta1() {
    double func;
    func =  +Jacobian[2][0] * theta1 + Jacobian[2][1] * theta2 + Jacobian[2][2] * omega1 + Jacobian[2][3] * omega2;
    return (func);
}

double jacobian_double_pendulum_func::JdTheta2() {
    double func;
    func = + Jacobian[3][0] * theta1 + Jacobian[3][1] * theta2 + Jacobian[3][2] * omega1
            + Jacobian[3][3] * omega2;
    return (func);
}

double jacobian_double_pendulum_func::JdOmega1() {
    double func;
    func =  + Jacobian[0][0] * theta1 + Jacobian[0][1] * theta2 + Jacobian[0][2] * omega1
            + Jacobian[0][3] * omega2;
    return (func);
}

double jacobian_double_pendulum_func::JdOmega2() {
    double func;
    func =  + Jacobian[1][0] * theta1 + Jacobian[1][1] * theta2 + Jacobian[1][2] * omega1
            + Jacobian[1][3] * omega2;
    return (func);
}
