#include "functions.h"

value functions_capsule::get_result(unsigned i) const{
    return result[i];
}

value functions_capsule::get_result(std::string) const{
}

const container & functions_capsule::get_result() const{
    return result;
}

unsigned functions_capsule::size() const{
    return result.size();
}

