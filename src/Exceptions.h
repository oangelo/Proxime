/* 
 * File:   Exceptions.h
 * Author: angelo
 *
 * Created on 4 de Junho de 2010, 15:55
 */

#ifndef _EXCEPTIONS_H
#define	_EXCEPTIONS_H

#include <iostream>
#include <string>
/*
 * Should happens when the user try to access out of range index
 */
class Index_error: public std::exception
{
    public:

        Index_error(std::string message) throw() ;
        ~Index_error() throw ();

        virtual const char* what() const throw();

    protected:
        std::string __message;

};
/*
 * Should happens when the user pass to a function/method an invalid parameter
 */

class  Value_error: public std::exception
{
    public:

         Value_error(std::string message) throw() ;
        ~ Value_error() throw ();

        virtual const char* what() const throw();

    protected:
        std::string __message;

};


#endif	/* _EXCEPTIONS_H */

