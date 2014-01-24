#include "simple_pendulum.h"

SimplePendulumFunction::SimplePendulumFunction()
:functions_capsule("Simple Pendulum", 3, 
                  dictionary{{"theta",0},{"omega",1}},
                  dictionary{{"l",0},{"g",1},}),
theta(), omega(), l(), g()
{
}

void SimplePendulumFunction::set(value& t, container& variables, container& parameters){ 
    theta = variables[index_var["theta"]];
    omega = variables[index_var["omega"]];
    
    l = parameters[index_par["l"]];
    g = parameters[index_par["g"]];
    
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
:functions_capsule("Simple Pendulum", 2, 
                  dictionary{{"q",0},{"p",1}},
                  dictionary{{"l",0},{"g",1},}),
p(), q(), l(), g(), m()
{
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
