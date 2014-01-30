#include "bifurcation_diagram.h"

bool CrossUpDown(value before, value after, value reference){
    if(before > reference and after < reference)
        return true;
    else
        return false;
};

bool CrossDownUp(value before, value after, value reference){
    if(before < reference and after > reference)
        return true;
    else
        return false;
};

container PhasePlaneSection(NumericalIntegration& attractor, unsigned coordinate_x,
        unsigned coordinate_y, value y_section_value, unsigned qt_points, bool (*cross)(value, value, value)){

    container value(attractor.size()), next_value(attractor.size());
    container zeros;
    size_t max_iterations(pow(10,7));

    next_value = attractor.get_variable();
    
    size_t counter(0);
    while(zeros.size() < qt_points){
        counter++;
        if(counter == max_iterations && zeros.size() == 0){
            std::cerr << "No Phase Plane Crosses Detected!" << std::endl;
            return zeros;
        }
        value = next_value;
        try {
            ++attractor;
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
