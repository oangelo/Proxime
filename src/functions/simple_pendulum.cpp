#include "simple_pendulum.h"

SimplePendulumFunction::SimplePendulumFunction()
:theta(),
l(), g()
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
