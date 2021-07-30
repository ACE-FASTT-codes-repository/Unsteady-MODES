/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------ Copyright (C) 2017 University of Strathclyde and Authors ------
-------------------- e-mail: c.greco@strath.ac.uk --------------------
----------------------- Author: Cristian Greco -----------------------
*/

#ifndef SMART_UQ_BASE_SURROGATE_H
#define SMART_UQ_BASE_SURROGATE_H

#include <vector>
#include "exception.h"

namespace smartuq {
namespace surrogate {


    /**
     * The base_surrogate is a template abstract class.
     *
     * @brief Any surrogate added to the toolbox needs to inherit from it and implement the virtual functions.
     */
    class base_surrogate {

    protected:

        // Set of grid points
        std::vector<std::vector<double>> m_X;

        // Set of responses
        std::vector<double> m_Y;

    public:


        /**
         * base_filter default constructor
         *
         * @brief vectors m_X and m_Y are initialized as empty vector by default
         */
        base_surrogate();

        /**
         * base_filter constructor
         *
         * @brief vectors m_X and m_Y are initialized as the input vectors
         */
        base_surrogate( const std::vector<std::vector<double>> &X,
                        const std::vector<double> &Y );


        /**
         * base_filter destructor
         */
        virtual ~base_surrogate() = default;


        /**
         * build pure virtual function to be implemented in derived class
         *
         * @brief Function that build the surrogate once a set of points and responses has been given
         */
        virtual void build () = 0;

        /**
         * evaluate pure virtual function to be implemented in derived class
         *
         * @brief Function that evaluate a new single point x with the constructed surrogate to compute the response y
         */
        virtual void evaluate (const std::vector<double> &x, double &y) const = 0;

        /**
         * evaluate function that call evaluate for a single point
         *
         * @brief Function that evaluate a new set of points X with the constructed surrogate to compute the responses Y
         */
        void evaluate (const std::vector<std::vector<double>> &X, std::vector<double> &Y) const;

        /**
         * Function that uses evaluate but return y as output
         */
        double evaluateReturn ( const std::vector<double> &x ) const ;


    public:

        /**
         * Get training points
         */
        std::vector<std::vector<double>> get_X() const;

        /**
         * Get responses
         */
        std::vector<double> get_Y() const;

    }; // class base_surrogate

} // surrogate
} // smartuq

#endif //SMART_UQ_BASE_SURROGATE_H
