/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------ Copyright (C) 2017 University of Strathclyde and Authors ------
-------------------- e-mail: c.greco@strath.ac.uk --------------------
----------------------- Author: Cristian Greco -----------------------
*/

#include "../include/Surrogates/base_surrogate.h"

using namespace smartuq;
using namespace surrogate;


/**
 * base_filter default constructor
 *
 * @brief vectors m_X and m_Y are initialized as empty vector by default
 */
base_surrogate::base_surrogate() {}


/**
 * base_filter constructor
 *
 * @brief vectors m_X and m_Y are initialized as the input vectors
 */
base_surrogate::base_surrogate( const std::vector<std::vector<double>> &X,
                                const std::vector<double> &Y ) :
    m_X(X), m_Y(Y)
{
    if (X.size()!=Y.size())
        smartuq_throw("Input set of points X and responses Y have mismatching dimension");
}


/**
 * base_filter destructor
 */
//base_surrogate::~base_surrogate() {}


/**
 * evaluate function that call evaluate for a single point
 *
 * @brief Function that evaluate a new set of points with the constructed surrogate
 */
void base_surrogate::evaluate(const std::vector<std::vector<double>> &X, std::vector<double> &Y) const {


    // Time to execute worhp solver
    clock_t begin, end;
    begin=clock();

    unsigned int n_points = X.size();

    Y.clear();
    Y.reserve(n_points);

    double y ;
    for ( unsigned int i = 0 ; i < n_points ; i++ ) {
        std::vector<double> x(X[i]);
        evaluate(x,y);
        Y.push_back(y);
    }


    //timing
    end=clock();
    double time_elapsed = (double (end-begin))/CLOCKS_PER_SEC;
    std::cout << "time elapsed to evaluate " << n_points << " points with the surrogate: " << time_elapsed << std::endl;

}

/**
 * Function that uses evaluate but return y as output
 */
double base_surrogate::evaluateReturn ( const std::vector<double> &x ) const {

    double y;
    this->evaluate(x,y);

    return y;
}



/**
 * Get training points
 */
std::vector<std::vector<double>> base_surrogate::get_X() const {

    return m_X;
}

/**
 * Get responses
 */
std::vector<double> base_surrogate::get_Y() const {

    return m_Y;
}
