#include "double_pendulum.h"

DoublePendulumFunction::DoublePendulumFunction(labels_values parameters):
FunctionCapsule("Double Pendulum",  
                  dictionary{{"theta1",0},{"theta2",1},{"omega1",2},{"omega2",3},},
                  dictionary{{"l1",0},{"l2",1},{"m1",2},{"m2",3},{"g",4},}, parameters),
theta1(), theta2(), omega1(), omega2(),
l1(parameters["l1"]), l2(parameters["l2"]), 
m1(parameters["m1"]), m2(parameters["m2"]), 
g(parameters["g"])
{} 

DoublePendulumFunction* DoublePendulumFunction::Clone() const{ 
    return new DoublePendulumFunction(*this);
}

DoublePendulumFunction* DoublePendulumFunction::Create(labels_values parameters) const{ 
    return new DoublePendulumFunction(parameters);
}

void DoublePendulumFunction::set(value& t, container& variables){ 
    theta1 = variables[index_var["theta1"]];
    theta2 = variables[index_var["theta2"]];
    omega1 = variables[index_var["omega1"]];
    omega2 = variables[index_var["omega2"]];
    
    result[index_var["theta1"]] = dTheta1();
    result[index_var["theta2"]] = dTheta2();
    result[index_var["omega1"]] = dOmega1();
    result[index_var["omega2"]] = dOmega2();
}

value DoublePendulumFunction::dTheta1() {
    value func;
    func = omega1;
    return (func);
}

value DoublePendulumFunction::dTheta2() {
    value func;
    func = omega2;
    return (func);
}

value DoublePendulumFunction::dOmega1() {
    value func;
    func = -(-m2 * g * sin(theta1) - g * sin(theta1) * m1 - m2 * l2 * pow(omega2, 0.2e1) * sin(theta1 - theta2) 
            + m2 * cos(theta1 - theta2) * g * sin(theta2) - m2 * cos(theta1 - theta2) * l1 * pow(omega1, 0.2e1) 
            * sin(theta1 - theta2)) / l1 / (-m2 - m1 + m2 * pow(cos(theta1 - theta2), 0.2e1));
    return (func);
}

value DoublePendulumFunction::dOmega2() {
    value func;
    func =  -(m2 * cos(theta1 - theta2) * g * sin(theta1) + cos(theta1 - theta2) * g * sin(theta1) 
            * m1 + cos(theta1 - theta2) * m2 * l2 * pow(omega2, 0.2e1) * sin(theta1 - theta2) - m2 * g 
            * sin(theta2) - g * sin(theta2) * m1 + l1 * pow(omega1, 0.2e1) * sin(theta1 - theta2) * m2 + l1 
            * pow(omega1, 0.2e1) * sin(theta1 - theta2) * m1) / l2 / (-m2 - m1 + m2 * pow(cos(theta1 - theta2), 0.2e1));
    return (func);

}

Jacobian_DoublePendulumFunction::Jacobian_DoublePendulumFunction(labels_values parameters): 
FunctionCapsule("Double Pendulum Jacobian", 
                  dictionary{{"theta1",0},{"theta2",1},{"omega1",2},{"omega2",3},},
                  dictionary{{"l1",0},{"l2",1},{"m1",2},{"m2",3},{"g",4},{"theta1",5},{"theta2",6},{"omega1",7},{"omega2",8},},
                  parameters),
theta1(), theta2(), omega1(), omega2(),
l1(parameters["l1"]), l2(parameters["l2"]), 
m1(parameters["m1"]), m2(parameters["m2"]), 
g(parameters["g"]),
Jacobian()
{
    value _theta1 = parameters["theta1"];
    value _theta2 = parameters["theta2"];
    value _omega1 = parameters["omega1"];
    value _omega2 = parameters["omega2"];
    Matrix_Jacob(_theta1,_theta2,_omega1,_omega2);
}

void Jacobian_DoublePendulumFunction::set(value& t, container& variables){

    theta1 = variables[index_var["theta1"]];
    theta2 = variables[index_var["theta2"]];
    omega1 = variables[index_var["omega1"]];
    omega2 = variables[index_var["omega2"]];

    result[index_var["theta1"]] = JdTheta1();
    result[index_var["theta2"]] = JdTheta2();
    result[index_var["omega1"]] = JdOmega1();
    result[index_var["omega2"]] = JdOmega2();
}

Jacobian_DoublePendulumFunction* Jacobian_DoublePendulumFunction::Clone() const{ 
    return new Jacobian_DoublePendulumFunction(*this);
}

Jacobian_DoublePendulumFunction* Jacobian_DoublePendulumFunction::Create(labels_values parameters) const{ 
    return new Jacobian_DoublePendulumFunction(parameters);
}

void Jacobian_DoublePendulumFunction::Matrix_Jacob(value theta1, value theta2,value omega1, value omega2)
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

value Jacobian_DoublePendulumFunction::JdTheta1() {
    value func;
    func =  +Jacobian[2][0] * theta1 + Jacobian[2][1] * theta2 + Jacobian[2][2] * omega1 + Jacobian[2][3] * omega2;
    return (func);
}

value Jacobian_DoublePendulumFunction::JdTheta2() {
    value func;
    func = + Jacobian[3][0] * theta1 + Jacobian[3][1] * theta2 + Jacobian[3][2] * omega1
            + Jacobian[3][3] * omega2;
    return (func);
}

value Jacobian_DoublePendulumFunction::JdOmega1() {
    value func;
    func =  + Jacobian[0][0] * theta1 + Jacobian[0][1] * theta2 + Jacobian[0][2] * omega1
            + Jacobian[0][3] * omega2;
    return (func);
}

value Jacobian_DoublePendulumFunction::JdOmega2() {
    value func;
    func =  + Jacobian[1][0] * theta1 + Jacobian[1][1] * theta2 + Jacobian[1][2] * omega1
            + Jacobian[1][3] * omega2;
    return (func);
}

value DoublePendulumEnergy(value theta1, value theta2, value omega1, value omega2,
value l1, value l2, value m1, value m2, value g){
    value energy = m1 * (l1 * l1 * pow(cos(theta1), 0.2e1) * pow(omega1, 0.2e1) +
    l1 * l1 * pow(sin(theta1), 0.2e1) * pow(omega1, 0.2e1)) / 0.2e1 + m2
    * (pow(l1 * cos(theta1) * omega1 + l2 * cos(theta2) *
    omega2, 0.2e1) + pow(l1 * sin(theta1) * omega1 + l2 *
    sin(theta2) * omega2, 0.2e1)) / 0.2e1 - m1 * g * l1 *
    cos(theta1) + m2 * g * (-l1 * cos(theta1) - l2 *
    cos(theta2));
    return energy;
}

