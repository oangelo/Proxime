#include "numerical_integration/adams_bashforth.h"

//function that teturns the x cross value
type_container AttractorCrossAxis(NumericalIntegration &attractor, int steps,
        int transiente, int coordinate_x = 0, int coordinate_y = 1, type_data y_coordinate_value = 0, int x_side = 1);
