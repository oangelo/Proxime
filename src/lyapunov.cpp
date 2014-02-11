#include "lyapunov.h"


value Dot(labels_values v1, labels_values v2) {
    value escalar = 0;
    for (labels_values::iterator it(v1.begin()); it != v1.end(); ++it){
        std::string cont((*it).first);
        escalar += v1[cont] * v2[cont];
    }
    return (escalar);
}

value MaxLyapunov::LocalLyapunov(){
    ++(*fiducial);
    labels_values result(fiducial->get_labels_values()); 
    for(auto i : result)
        jacobian_parameters[i.first] = i.second;
    //Integrating the linear equations and saving the vector on ortogonal_space
    std::unique_ptr<FunctionCapsule> function(jacobian->Create(jacobian_parameters));
    RungeKutta4Th linearized_sistem(*function, base, fiducial->get_dt());
    ++linearized_sistem;
    base = linearized_sistem.get_labels_values();
    //Ortonormalizing the vectors and saving the modulus
    value modulo(sqrt(Dot(base, base)));
    for (labels_values::iterator it(base.begin()); it != base.end(); ++it)
        base[(*it).first] = base[(*it).first] / modulo;
    return modulo;
}

MaxLyapunov::MaxLyapunov(NumericalIntegration& model, FunctionCapsule& _jacobian, 
                         labels_values jacobian_parameters, unsigned transient)
:fiducial(model.Clone()), jacobian(_jacobian.Clone()), exponent(0), base(), jacobian_parameters(jacobian_parameters)
{
    labels_values variables(fiducial->get_labels_values());
    labels_values::iterator it(variables.begin());
    base[(*it).first] = 1;
    ++it;
    while(it != variables.end()) {
        base[(*it).first] = 0;
        ++it;
    }
    for (unsigned steps = 0; steps < transient; steps++) {
        LocalLyapunov();
    }
}

value MaxLyapunov::operator()(unsigned iterations){
    value sum(0);
    for (unsigned steps(0); steps < iterations; steps++) {
        sum += log(LocalLyapunov());
    }
    return sum / (iterations * fiducial->get_dt());
}

//####################### Lyapunov Spectrum ############################

container LyapunovSpectrum::LocalLyapunov(){
    ++(*fiducial);
    labels_values result(fiducial->get_labels_values()); 
    for(auto i : result)
        jacobian_parameters[i.first] = i.second;
    //Integrating the linear equations and saving the vector on ortogonal_space
    std::unique_ptr<FunctionCapsule> function(jacobian->Create(jacobian_parameters));
    for(std::vector<labels_values>::iterator it(base.begin()); it != base.end(); ++it){
        RungeKutta4Th linearized_sistem(*function, *it, fiducial->get_dt());
        ++linearized_sistem;
        *it = linearized_sistem.get_labels_values();
    }
    return GramSchmidt();
}

LyapunovSpectrum::LyapunovSpectrum(NumericalIntegration& model, FunctionCapsule& _jacobian, 
                         labels_values jacobian_parameters, unsigned transient)
:fiducial(model.Clone()), jacobian(_jacobian.Clone()), exponent(model.size(), 0), base(model.size()), jacobian_parameters(jacobian_parameters)
{
    labels_values variables(fiducial->get_labels_values());
    for(labels_values::iterator _it(variables.begin()); _it != variables.end(); ++_it){
        for(labels_values::iterator it(variables.begin()); it != variables.end(); ++it){
            size_t count(std::distance(variables.begin(), _it));
            if(it == _it){ 
                base[count][(*it).first] = 1;
            }else{
                base[count][(*it).first] = 0;
            }
        }
    }
    for (unsigned steps = 0; steps < transient; steps++) {
        LocalLyapunov();
    }
}

container LyapunovSpectrum::operator()(unsigned iterations){
    container sum(fiducial->size(), 0);
    for (unsigned steps(0); steps < iterations; steps++) {
        container aux(LocalLyapunov());
            for(size_t i(0); i < sum.size(); ++i){
                sum[i] += log(aux[i]);
            }
    }
    for(size_t i(0); i < sum.size(); ++i){
        sum[i] /= (iterations * fiducial->get_dt());
    }
    return sum;
}

container LyapunovSpectrum::GramSchmidt() {
    std::vector<labels_values> aux_base = base;
    container modulo(base.size());
    modulo[0] = sqrt(Dot(aux_base[0], aux_base[0]));
    for (size_t i(1); i < base.size(); ++i) {
        for (size_t j(0); j < i; ++j)
            for (labels_values::iterator it(base[0].begin()); it != base[0].end(); ++it){
                std::string label((*it).first);
                aux_base[i][label] += -(Dot(base[i], aux_base[j]) * aux_base[j][label]) / Dot(aux_base[j], aux_base[j]);
            }
        modulo[i] = sqrt(Dot(aux_base[i], aux_base[i]));
    }
    //normalizing###################################################
    for (unsigned i = 0; i < base.size(); i++)
        for (labels_values::iterator it(base[0].begin()); it != base[0].end(); ++it){
            std::string label((*it).first);
            aux_base[i][label] = aux_base[i][label] / modulo[i];
        }
    //##############################################################
    base = aux_base;
    return modulo;
}
