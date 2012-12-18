#ifndef _NUMERICLA_INTEGRATION_
#define _NUMERICLA_INTEGRATION_

#include <vector>
#include <string>
#include <memory>

#include "exceptions.h"
#include "../functions.h"


class NumericalIntegration {
    public:
       NumericalIntegration(type_container variable, type_container parameter, type_data dt);
        virtual ~NumericalIntegration(){};
        type_data get_dt() const;
        type_data get_t() const;
        type_data get_variable(unsigned n) const;
        const type_container & get_variable() const {return(__variable);};
        type_data get_parameter(unsigned n) const;
        const std::string & get_model_name() const {return(__model_name);};
        const std::string & get_method_name() const {return(__method);};

        unsigned size_variable() const;
        unsigned size_parameter() const;

        virtual void next() = 0;

        //return the variables
        type_data operator[] (const unsigned nIndex);
        type_data operator[] (std::string nIndex);
        //return the variable ready to print
        friend std::ostream& operator<< (std::ostream &out, NumericalIntegration &object);

    protected:

        std::auto_ptr<functions_capsule> __func;
        type_container __variable;
        type_container __parameter;
        type_data __h, __t;
        std::string __model_name;
        std::string __method;

};

#endif //_NUMERICLA_INTEGRATION_
