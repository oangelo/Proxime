#include "simple_pendulum.h"

SimplePendulumFunction::SimplePendulumFunction()
:theta(), omega(), l(), g()
{
    __result.clear();
    __result.resize(2);
}

void SimplePendulumFunction::set(value& t, container& variables, container& parameters){ 
    theta = variables[V_THETA];
    omega = variables[V_OMEGA];
    
    l = parameters[P_L];
    g = parameters[P_G];
    
    __result[V_THETA]=dTheta();
    __result[V_OMEGA]=dOmega();
   
}

value SimplePendulumFunction::dTheta() {
    value func;
    func = omega;
    return (func);
}

value SimplePendulumFunction::dOmega() {
    value func;
    func = -(g/l) * sin(theta);
    return (func);
}

value SimplePendulumEnergy(value theta, value omega, value l, value m, value g){
    return m * pow(l * omega, 2) / 2 - m * g * (l * cos(theta));
}

SimplePendulum_H::SimplePendulum_H()
:p(), q(), l(), g(), m()
{
    __result.clear();
    __result.resize(2);
}

void SimplePendulum_H::set(value& t, container& variables, container& parameters){ 
    p = variables[V_P];
    q = variables[V_Q];
    
    l = parameters[P_L];
    g = parameters[P_G];
    m = parameters[P_M];
    
    __result[0]=Vq();
    __result[1]=Tp();
   
}

value SimplePendulum_H::Vq() {
    value func;
    func =  l * m * g * sin(q);
    return (func);
}

value SimplePendulum_H::Tp() {
    value func;
    func =  p / (m * pow(l, 2));
    return (func);
}
