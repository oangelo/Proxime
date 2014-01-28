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
       NumericalIntegration(FunctionCapsule& function, container variable, value dt);
        virtual ~NumericalIntegration(){};

        value get_dt() const;
        value get_t() const;
        const container & get_variable() const {return(variable);};
        const std::string & get_model_name() const {return(model_name);};
        const std::string & get_method_name() const {return(method);};
        labels_and_values GetLabelsValues();

        unsigned size_variable() const;

        virtual void next() = 0;

        //return the variables
        value operator[] (std::string nIndex);
        //return the variable ready to print
        friend std::ostream& operator<< (std::ostream &out, NumericalIntegration &object);

    protected:

        FunctionCapsule* function;
        container variable;
        value dt, time;
        std::string model_name;
        std::string method;
};

#endif //_NUMERICLA_INTEGRATION_
