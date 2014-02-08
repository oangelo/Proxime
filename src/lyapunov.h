#include <vector>
#include <memory>

#include "numerical_integration/runge_kutta.h"   


class MaxLyapunov{
public:
    MaxLyapunov(NumericalIntegration& model, FunctionCapsule& _jacobian, labels_values jacobian_parameters, unsigned transient);
    value operator()(unsigned iterations);
    const value get_exponent() const;
private:
    value LocalLyapunov();

    std::unique_ptr<NumericalIntegration> fiducial;
    std::unique_ptr<FunctionCapsule> jacobian;
    value exponent;
    labels_values base;
    labels_values jacobian_parameters;
};

class LyapunovSpectrum{
public:
    LyapunovSpectrum(NumericalIntegration& model, FunctionCapsule& _jacobian, labels_values jacobian_parameters, unsigned transient);
    container operator()(unsigned iterations);
    const value get_exponent() const;
private:
    container LocalLyapunov();
    container GramSchmidt() ;

    std::unique_ptr<NumericalIntegration> fiducial;
    std::unique_ptr<FunctionCapsule> jacobian;
    container exponent;
    std::vector<labels_values> base;
    labels_values jacobian_parameters;
};
