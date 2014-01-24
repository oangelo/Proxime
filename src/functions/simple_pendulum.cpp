#include "simple_pendulum.h"

SimplePendulumFunction::SimplePendulumFunction()
:functions_capsule("Simple Pendulum", 2),
theta(), omega(), l(), g()
{
}

void SimplePendulumFunction::set(value& t, container& variables, container& parameters){ 
    theta = variables[V_THETA];
    omega = variables[V_OMEGA];
    
    l = parameters[P_L];
    g = parameters[P_G];
    
    result[V_THETA]=dTheta();
    result[V_OMEGA]=dOmega();
   
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
:functions_capsule("Simple Pendulum Hamiltonian", 2),
p(), q(), l(), g(), m()
{
    result.clear();
    result.resize(2);
}

void SimplePendulum_H::set(value& t, container& variables, container& parameters){ 
    p = variables[V_P];
    q = variables[V_Q];
    
    l = parameters[P_L];
    g = parameters[P_G];
    m = parameters[P_M];
    
    result[0]=Vq();
    result[1]=Tp();
   
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
