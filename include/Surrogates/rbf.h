/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------ Copyright (C) 2017 University of Strathclyde and Authors ------
-------------------- e-mail: c.greco@strath.ac.uk --------------------
----------------------- Author: Cristian Greco -----------------------
*/

#ifndef SMART_UQ_RBF_H
#define SMART_UQ_RBF_H

#include <cmath>
#include <iostream>

#include "base_surrogate.h"
#include "../Eigen/Core"
#include "../Eigen/Dense"


namespace smartuq {
namespace surrogate {

    /**
     * Basis Function for RBF surrogate
     */
    enum RBF_FUNCTION {
        LINEAR,
        CUBIC,
        GAUSSIAN,
        THIN_PLATE,
        MULTIQUADRATICS
    };

    /**
     * Constants required by RBF
     */
    struct RBF_CONSTANTS {
        RBF_CONSTANTS ( double avgd,
                        double smth ) :
            avgdist(avgd), smooth(smth)
        {}

        double avgdist;
        double smooth;
    };

    /**
     * Constraints required when RBF THIN_PLATE is used for volume delimitation
     * and scalar check for in/out point condition
     *
     * interior: a negative value constraint placed at the cloud barycenter - suggested for planned model construction
     * exterior: positive value constraints placed outside the cloud volume - suggested for interactive modelling
     *    -- short idea: Take maximum distance r_max border point, multiply by percentage r = r_max*exteriorPercRadius, construct constraints on hyper-circle with radius r
     * normal  : negative value constraints placed inside the cloud, normal to the border - suggested for conversion from polygons
     *    -- short idea: Take minimum distance r_min border point, multiply by percentage r = r_min*normalPercRadius, construct constraints on hyper-circle with radius r on direction of border
     *    -- Disclaimer: It's not properly the NORMAL interpretation because it follows the direction defined by barycenter-border
     *
     * WARNING: INTERIOR WILL NOT WORK FOR NON-CONVEX VOLUMES
     * WARNING: NORMAL COULD NOT WORK FOR NON-CONVEX NICE VOLUMES
     */
    struct RBF_CONSTRAINTS {

        // Flags to tell which shall be used
        bool interior = false;
        bool exterior = true;
        bool normal = false;

        // Percentage distance with respect to maximum/minimum border distance
        std::vector<double> exteriorPercRadius = std::vector<double> (1,0.2);
        double normalPercRadius = 0.2;
    };


    /**
     * Routine to ADD normal/interior/exterior constraints to input borderX according to struct RBF_CONSTRAINTS
     * AND to compute responses ( 0.0 at border, !=0.0 at constraints )
     *
     * X : input  = border of delimited surface
     * X : output = border + constraints points
     * Y : output = border responses (0.0) + constraints responses
     */
    void rbf_addConstraints ( const RBF_CONSTRAINTS &rbfConstraints,
                              std::vector<std::vector<double>> &X,
                              std::vector<double> &Y ) ;


    /**
     * rbf is a publicly derived class from base_surrogate
     *
     * @brief Class that implement surrogate based on radial basis network
     *        The network is constituted by the weighted sum of a radial basis function centered in point x
     *        with a given response y. Therefore, there are as many basis functions as number of points. The
     *        basis functions are radial as their value depends only on the distance from the center.
     */
    class rbf : public base_surrogate {

    protected:

        // Class scope constant [for readibility, declaration at the end]
        static const RBF_CONSTANTS rbf_std_constants;

        // Basis function for surrogate
        RBF_FUNCTION m_rbf_function;

        // Coefficients/Weights of different radial functions
        std::vector<double> m_coeff;

        // Constants for radial basis
        RBF_CONSTANTS m_constants;

    public:

        /**
         * Constructor for rbf class
         *
         * @param rbf_function : Type of radial function
         */
        explicit rbf ( const RBF_FUNCTION &rbf_function,
                const RBF_CONSTANTS &constants = rbf_std_constants );

        /**
         * Constructor for rbf class
         *
         * @param X :            Set of points in which independent rbf are centered
         * @param Y :            Respo;
         * nse of independent rbf in points X
         * @param rbf_function : Type of radial function
         */
        rbf ( const std::vector<std::vector<double>> &X,
              const std::vector<double> &Y,
              const RBF_FUNCTION &rbf_function,
              const RBF_CONSTANTS &constants = rbf_std_constants );

//        virtual ~rbf();

        /**
         * Function to set basis points and responses
         *
         * @param X : Set of points in which independent rbf are centered
         * @param Y : Response of independent rbf in points X
         */
        void set_points( const std::vector<std::vector<double>> &X,
                         const std::vector<double> &Y ) ;



        /**
         * build virtual function implementation
         *
         * @brief Function that build the surrogate once a set of points and responses has been given
         */
        void build () override;

        /**
         * evaluate virtual function implementation
         *
         * @brief Function that evaluate a new single point x with the constructed surrogate to compute the response y
         */
        void evaluate (const std::vector<double> &x, double &y) const override;


        /**
         * Get rbf coefficient vector
         */
        std::vector<double> get_coeff() const ;


    protected:

        void evaluate_basis( std::vector<std::vector<double>> &R ) const ;

        void evaluate_basis( std::vector<double> &r ) const ;

    }; // class rbf



} // surrogate
} // smartuq

#endif //SMART_UQ_RBF_H
