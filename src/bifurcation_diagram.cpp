#include "bifurcation_diagram.h"

bool CrossUpDown(double before, double after, double reference){
    if(before > reference and after < reference)
        return true;
    else
        return false;
};

bool CrossDownUp(double before, double after, double reference){
    if(before < reference and after > reference)
        return true;
    else
        return false;
};

type_container PhasePlaneSection(NumericalIntegration &attractor, unsigned transient, unsigned iterations, 
        unsigned coordinate_x, unsigned coordinate_y, double y_section_value, bool (*cross)(double, double, double)){
    type_container value(attractor.size_variable()), next_value(attractor.size_variable());
    type_container zeros;

    for (int i = 0; i < transient; i++){
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
    for (int i=0; i < iterations; i++){
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
        if(cross(value[coordinate_y], next_value[coordinate_y], y_section_value))
        {
            zeros.push_back(value[coordinate_x]);
        }
    }
    return(zeros);
}
