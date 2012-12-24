#include "bifurcation_diagram.h"


type_container AttractorCrossAxis(NumericalIntegration &attractor, int steps, int transiente,
        int coordinate_x, int coordinate_y, type_data y_coordinate_value, int x_side)
{
    type_container value(attractor.size_variable()), next_value(attractor.size_variable());
    type_container zeros;

    for (int i = 0; i < transiente; i++){
        try {
            attractor.next();
        }
        catch (Value_error) {
            std::cerr << "Problem with the numerical integration, "
                "please, use a smaller step or a better method" << std::endl;
            zeros.clear();
            return(zeros);
        }
    }
    next_value = attractor.get_variable();
    for (int i=0; i < steps; i++){
        value = next_value;
        try {
            attractor.next();
        }
        catch (Value_error) {
            std::cout << "Problem with the numerical integration, "
                "please use a smaller step or a better method" << std::endl;
            zeros.clear();
            return(zeros);
        }
        next_value = attractor.get_variable();
        if (x_side * value[coordinate_x] > 0) {
            if (((value[coordinate_y] > y_coordinate_value) && (next_value[coordinate_y] < y_coordinate_value))
                    || ((value[coordinate_y] < y_coordinate_value) && (next_value[coordinate_y] > y_coordinate_value)))
            {
                zeros.push_back(value[coordinate_x]);
            }
        }
    }
    return(zeros);
}
