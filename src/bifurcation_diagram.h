#include "numerical_integration/adams_bashforth.h"

/* This function computes the Poincaré section of a phase plane 
   at Y = y_section_values and retuns the X values. This can be
   used to build a bifurcation diagram.
*/
type_container PhasePlaneSection(NumericalIntegration &attractor, unsigned transient, unsigned iterations, 
        unsigned coordinate_x, unsigned coordinate_y, double y_section_value, bool (*cross)(double, double, double));


bool CrossUpDown(double before, double after, double reference);
bool CrossUpDown(double before, double after, double reference);
