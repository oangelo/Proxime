#ifndef _NUMERICLA_INTEGRATION_
#define _NUMERICLA_INTEGRATION_

#include <vector>
#include <string>
#include <memory>
#include <iomanip>

#include "exceptions.h"
#include "../functions/functions.h"


class NumericalIntegration {
    public:
       NumericalIntegration(container variable, container parameter, value dt);
        virtual ~NumericalIntegration(){};
        value get_dt() const;
        value get_t() const;
        value get_variable(unsigned n) const;
        const container & get_variable() const {return(__variable);};
        value get_parameter(unsigned n) const;
        const std::string & get_model_name() const {return(__model_name);};
        const std::string & get_method_name() const {return(__method);};

        unsigned size_variable() const;
        unsigned size_parameter() const;

        virtual void next() = 0;

        //return the variables
        value operator[] (const unsigned nIndex);
        value operator[] (std::string nIndex);
        //return the variable ready to print
        friend std::ostream& operator<< (std::ostream &out, NumericalIntegration &object);

    protected:

        std::auto_ptr<functions_capsule> __func;
        container __variable;
        container __parameter;
        value __h, __t;
        std::string __model_name;
        std::string __method;
};

#endif //_NUMERICLA_INTEGRATION_
