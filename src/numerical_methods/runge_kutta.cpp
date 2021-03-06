#include"runge_kutta.h"

RungeKutta4Th::RungeKutta4Th(FunctionCapsule & function, labels_values variable,value dt)
    :NumericalIntegration(function, variable, dt) { } 

RungeKutta4Th* RungeKutta4Th::Clone() const{
    return new RungeKutta4Th(*this);
}
RungeKutta4Th* RungeKutta4Th::Create(FunctionCapsule & function, labels_values variable,value dt) const {
    return new RungeKutta4Th(function, variable, dt);
}

RungeKutta4Th& RungeKutta4Th::operator++(){
    RungeKutta4Th_method();
    return *this;
}

void RungeKutta4Th::RungeKutta4Th_method()
{
    //value k[variable.size()][4];
    std::vector<container> k(variable.size(), container(4, 0));
    value t = time;

    unsigned i;
    container aux_argument(variable.size(), 0);


    //Generating K1*************************************************************
    function->set(time, variable);
    for (i = 0; i < variable.size(); i++) {
        k[i][0] = function->get_result(i);
        //k[i][0] = ((*functiontion[i])(time, variable, __parameter));
    }

    //Generating K2*************************************************************
    t = time + 0.5 * dt;
    for (i = 0; i < variable.size(); i++) {
        aux_argument[i] = variable[i]+(0.5 * dt * k[i][0]);
    }
    function->set(t, aux_argument);
    for (i = 0; i < variable.size(); i++) {
        k[i][1] = function->get_result(i);
    }
    //Generating K3*************************************************************
    t = time + 0.5 * dt;
    for (i = 0; i < variable.size(); i++) {
        aux_argument[i] = variable[i]+(0.5 * dt * k[i][1]);
    }
    function->set(t, aux_argument);
    for (i = 0; i < variable.size(); i++) {
        k[i][2] = function->get_result(i);
    }
    //Generating K4*************************************************************
    t = time + dt;
    for (i = 0; i < variable.size(); i++) {
        aux_argument[i] = variable[i] + dt * k[i][2];
    }
    function->set(t, aux_argument);
    for (i = 0; i < variable.size(); i++) {
        k[i][3] = function->get_result(i);
    }

    //Generating the final values***********************************************
   
    time = time + dt;
    for (i = 0; i < variable.size(); i++) {
        variable[i] = variable[i] + dt * (k[i][0] + 2 * k[i][1] + 2 * k[i][2] + k[i][3]) / 6;
        if (variable[i] != variable[i]) {
            throw Value_error("Value error in the Runge-Kutta 4th order method integrating " + function->get_name());
        }
    }
   

}


