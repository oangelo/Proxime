Proxime - Ordinary Differential Equation Integration and Analysis
================================================================

This project intends to offer the time evolution of various models, 
and also some characteristics values or visualizations of the model orbits. 
However, this project will not offer time series analysis tools, 
since here the focus are models. 
Nevertheless, these model may be very useful to generate data to 
test and develop time series analysis tools.

This lib has the following:

Numerical methods:
* [Runge-Kutta](http://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods "wikipedia"), 4Th Order
* [Adams-Bashforth](http://en.wikipedia.org/wiki/Adams_Bashforth#Adams.E2.80.93Bashforth_methods "wikipedia"), 4Th Order
* [Adams-Moulton](http://en.wikipedia.org/wiki/Adams_Bashforth#Adams.E2.80.93Moulton_methods "wikipedia"), 4Th Order

These models:
* [Rössler Attractor](http://en.wikipedia.org/wiki/R%C3%B6ssler_attractor "wikipedia")
* [Lorenz Attractor](http://en.wikipedia.org/wiki/Lorenz_attractor "wikipedia")
* [Double Pendulum](http://en.wikipedia.org/wiki/Double_pendulum "wikipedia")

These Analysis:
* [Lyapunov Spectrum](http://en.wikipedia.org/wiki/Lyapunov_exponent "wikipedia")
* [Bifurcation Diagrams](http://en.wikipedia.org/wiki/Bifurcation_diagram "wikipedia")

Installation
============

This will install the header files(.h) to /usr/include and library files(.so) to /usr/lib and the binary to /usr/bin.

### Installation steps:
1. `make`
2. `sudo make install`

Code directives
===============

* The code should use [Google C++ Style](http://google-styleguide.googlecode.com/svn/trunk/cppguide.xml).
* This project uses [googletest](https://code.google.com/p/googletest/) for tests (folder unit_test).
* For consistence, the code should only use these variable types: 
value, container and labels, defined on src/functions/functions.h


