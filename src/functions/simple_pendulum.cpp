#include "simple_pendulum.h"

SimplePendulum::SimplePendulum(labels_values parameters)
:FunctionCapsule("Simple Pendulum", 
                  dictionary{{"theta",0},{"omega",1}},
                  dictionary{{"l",0},{"g",1},}, 
                  parameters),
theta(), omega(), l(), g()
{
    l = parameters["l"];
    g = parameters["g"];
}

SimplePendulum* SimplePendulum::Clone() const{
    return(new SimplePendulum(*this));
}

SimplePendulum* SimplePendulum::Create(labels_values parameters) const{
    return new SimplePendulum(parameters);
}

void SimplePendulum::set(value& t, container& variables){ 
    theta = variables[index_var["theta"]];
    omega = variables[index_var["omega"]];
   
    result[index_var["theta"]]=dTheta();
    result[index_var["theta"]]=dOmega();
}

value SimplePendulum::dTheta() {
    value func;
    func = omega;
    return (func);
}

value SimplePendulum::dOmega() {
    value func;
    func = -(g/l) * sin(theta);
    return (func);
}

value SimplePendulumEnergy(value theta, value omega, value l, value m, value g){
    return m * pow(l * omega, 2) / 2 - m * g * (l * cos(theta));
}

