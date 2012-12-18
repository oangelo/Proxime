/* 
 * File:   Exceptions.cpp
 * Author: angelo
 * 
 * Created on 4 de Junho de 2010, 15:55
 */

#include "exceptions.h"


Index_error::Index_error(std::string message) throw()
:message(message)
{
}

Index_error::~Index_error() throw(){}

const char* Index_error::what() const throw()
{
    return (message.c_str());
}

    Value_error::Value_error(std::string message) throw()
:message(message)
{
}

Value_error::~Value_error() throw(){}

const char* Value_error::what() const throw()
{
    return (message.c_str());
}

