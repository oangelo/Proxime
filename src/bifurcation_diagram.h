#include "numerical_methods/adams_bashforth.h"

bool CrossUpDown(value before, value after, value reference);
bool CrossDownUP(value before, value after, value reference);

/* This function computes the Poincaré section of a phase plane 
   at Y = y_section_values and retuns the X values. This can be
   used to build a bifurcation diagram.
*/
container PhasePlaneSection(NumericalIntegration& attractor, unsigned coordinate_x,
        unsigned coordinate_y, value y_section_value, unsigned qt_points=50, bool (*cross)(value, value, value)=&CrossUpDown);



