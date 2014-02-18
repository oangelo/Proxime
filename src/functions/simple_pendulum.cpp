#include "simple_pendulum.h"

SimplePendulumFunction::SimplePendulumFunction(labels_values parameters)
:FunctionCapsule("Simple Pendulum", 
                  dictionary{{"theta",0},{"omega",1}},
                  dictionary{{"l",0},{"g",1},}, 
                  parameters),
theta(), omega(), l(), g()
{
    l = parameters["l"];
    g = parameters["g"];
}

SimplePendulumFunction* SimplePendulumFunction::Clone() const{
    return(new SimplePendulumFunction(*this));
}

SimplePendulumFunction* SimplePendulumFunction::Create(labels_values parameters) const{
    return new SimplePendulumFunction(parameters);
}

void SimplePendulumFunction::set(value& t, container& variables){ 
    theta = variables[index_var["theta"]];
    omega = variables[index_var["omega"]];
   
    result[index_var["theta"]]=dTheta();
    result[index_var["theta"]]=dOmega();
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

