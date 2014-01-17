#include "numerical_integration/adams_bashforth.h"
bool CrossUpDown(double before, double after, double reference);
bool CrossDownUP(double before, double after, double reference);

/* This function computes the Poincaré section of a phase plane 
   at Y = y_section_values and retuns the X values. This can be
   used to build a bifurcation diagram.
*/
type_container PhasePlaneSection(NumericalIntegration& attractor, unsigned coordinate_x,
        unsigned coordinate_y, double y_section_value, unsigned qt_points=50, bool (*cross)(double, double, double)=&CrossUpDown);



