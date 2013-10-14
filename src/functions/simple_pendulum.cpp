#include "simple_pendulum.h"

SimplePendulumFunction::SimplePendulumFunction()
:theta(), omega(), l(), g()
{
    __result.clear();
    __result.resize(2);
}

void SimplePendulumFunction::set(double& t, std::vector<double>& variables, std::vector<double>& parameters){ 
    theta = variables[V_THETA];
    omega = variables[V_OMEGA];
    
    l = parameters[P_L];
    g = parameters[P_G];
    
    __result[V_THETA]=dTheta();
    __result[V_OMEGA]=dOmega();
   
}

double SimplePendulumFunction::dTheta() {
    double func;
    func = omega;
    return (func);
}

double SimplePendulumFunction::dOmega() {
    double func;
    func = -(g/l) * sin(theta);
    return (func);
}

double SimplePendulumEnergy(double theta, double omega, double l, double m, double g){
    return m * pow(l * omega, 2) / 2 - m * g * (l * cos(theta));
}

SimplePendulum_H::SimplePendulum_H()
:p(), q(), l(), g(), m()
{
    __result.clear();
    __result.resize(2);
}

void SimplePendulum_H::set(double& t, std::vector<double>& variables, std::vector<double>& parameters){ 
    p = variables[V_P];
    q = variables[V_Q];
    
    l = parameters[P_L];
    g = parameters[P_G];
    m = parameters[P_M];
    
    __result[0]=Vq();
    __result[1]=Tp();
   
}

double SimplePendulum_H::Vq() {
    double func;
    func =  l * m * g * sin(q);
    return (func);
}

double SimplePendulum_H::Tp() {
    double func;
    func =  p / (m * pow(l, 2));
    return (func);
}
