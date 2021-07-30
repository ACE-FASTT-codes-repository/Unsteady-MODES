/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
------------ e-mail: carlos.ortega@strath.ac.uk ----------------------
--------- Author: Annalisa Riccardi and Carlos Ortega Absil ----------
*/

#ifndef SMARTUQ_EXCEPTIONS_H
#define SMARTUQ_EXCEPTIONS_H

#include <exception>
#include <cassert>
#include <iostream>
#include <string>

namespace smartuq{

#define _SMARTUQ_EXCEPTION_QUOTEME(x) #x
#define SMARTUQ_EXCEPTION_QUOTEME(x) _SMARTUQ_EXCEPTION_QUOTEME(x)
#define SMARTUQ_EXCEPTION_EXCTOR(s) ((std::string(__FILE__ "," SMARTUQ_EXCEPTION_QUOTEME(__LINE__) ": ") + s) + ".")
#define SMARTUQ_EX_THROW(s) (throw smartuq::smartuq_exception(SMARTUQ_EXCEPTION_EXCTOR(s)))

#define smartuq_throw(s) SMARTUQ_EX_THROW(s)


class smartuq_exception: public std::exception {
	public:
		smartuq_exception(const std::string &s):m_what(std::string("SMART-UQ:")+s) {}
		virtual const char *what() const throw() {
			return m_what.c_str();
		}
		virtual ~smartuq_exception() throw() {}
	protected:
		std::string m_what;
};

}
#endif
