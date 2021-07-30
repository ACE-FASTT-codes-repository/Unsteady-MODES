/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------ Copyright (C) 2017 University of Strathclyde and Authors ------
-------------------- e-mail: c.greco@strath.ac.uk --------------------
----------------------- Author: Cristian Greco -----------------------
*/

#include "../include/Surrogates/rbf.h"

using namespace smartuq;
using namespace surrogate;



/**
 * Routine to ADD normal/interior/exterior constraints to input borderX according to struct RBF_CONSTRAINTS
 * AND to compute responses ( 0.0 at border, !=0.0 at constraints )
 *
 * X : input  = border of delimited surface
 * X : output = border + constraints points
 * Y : output = border responses (0.0) + constraints responses
 */
void smartuq::surrogate::rbf_addConstraints ( const RBF_CONSTRAINTS &rbfConstraints,
                                              std::vector<std::vector<double>> &X,
                                              std::vector<double> &Y ) {

    // number of variables and border samples
    unsigned int nBorder = X.size();
    unsigned int nVar = X[0].size();

    // Clear response vector
    Y.clear();

    // Border Barycenter computation
    std::vector<double> barycenter (nVar,0.0);
    for ( unsigned int i = 0 ; i < nBorder ; i++ ) {

        // Add response of border point
        Y.push_back(0.0);

        for (unsigned int j = 0; j < nVar; j++)
            barycenter[j] += X[i][j];
    }
    for ( unsigned int j = 0 ; j < nVar ; j++ )
        barycenter[j] /= (double)nBorder;

    std::cout << "\nBarycenter is ";
    for ( unsigned int j = 0 ; j < nVar ; j++ )
        std::cout << barycenter[j] << " ";
    std::cout << std::endl<<std::endl;

    if (rbfConstraints.interior) {
        X.push_back(barycenter);
        Y.push_back(-1.0);
    }

    if (rbfConstraints.exterior) {

        // Vector containint distances from barycenter
        std::vector<double> barDist;

        // Compute maximum distance from barycenter
        double maxBarDist = 0.0;
        for (unsigned int i = 0; i < nBorder; i++) {
            double dist = 0.0;
            for (unsigned int j = 0; j < nVar; j++)
                dist += std::pow(X[i][j], 2.0);
            dist = std::sqrt(dist);

            if (dist > maxBarDist)
                maxBarDist = dist;

            // Update distance vector
            barDist.push_back(dist);
        }

        unsigned int exterSize = rbfConstraints.exteriorPercRadius.size();

        for ( unsigned int k = 0 ; k < exterSize ; k++) {

            if (rbfConstraints.exteriorPercRadius[k] <= 0.05)
                smartuq_throw("exteriorPercRadius shall respect the condition 0.05<exteriorPercRadius");

            // Insert constraint at scaled distance
            double scaledDist = maxBarDist * (1.0 + rbfConstraints.exteriorPercRadius[k]);

            for (unsigned int i = 0; i < nBorder; i++) {

                // New constraint and response
                std::vector<double> constraint(X[i]);
                double constraintResponse;

                // Scale to new distance
                for (unsigned int j = 0; j < nVar; j++)
                    constraint[j] *= scaledDist / barDist[i];

                // Scale response according to distance of constraint position from border point
                constraintResponse = (1.0/(double)(exterSize-k)) * (scaledDist - barDist[i]) / (scaledDist - maxBarDist);
                if (constraintResponse <= 0.0)
                    smartuq_throw("constraintResponse for exterior constraint is smaller or equal than 0");

                // Add new constraint point and response
                X.push_back(constraint);
                Y.push_back(constraintResponse);
            }
        }

    }

    if (rbfConstraints.normal) {

        if ( rbfConstraints.normalPercRadius <=0.05 || rbfConstraints.normalPercRadius >0.95 )
            smartuq_throw("normalPercRadius shall respect the condition 0.05<normalPercRadius<=0.95");

        // Vector containint distances from barycenter
        std::vector<double> barDist ;

        // Compute maximum distance from barycenter
        double minBarDist = 1.0e10;
        for ( unsigned int i = 0 ; i < nBorder ; i++ ) {
            double dist = 0.0;
            for ( unsigned int j = 0 ; j < nVar ; j++ )
                dist += std::pow(X[i][j],2.0);
            dist = std::sqrt(dist);

            if ( dist<minBarDist )
                minBarDist = dist;

            // Update distance vector
            barDist.push_back(dist);
        }

        if (minBarDist>=9.99e9)
            smartuq_throw("Error minBarDist not changed from pre-set value");

        // Insert constraint at scaled distance
        if (rbfConstraints.exteriorPercRadius.size()!=1)
            smartuq_throw("size different from 1");

        double scaledDist = minBarDist*(1.0-rbfConstraints.exteriorPercRadius[0]);

        for ( unsigned int i = 0 ; i < nBorder ; i++ ) {

            // New constraint and response
            std::vector<double> constraint (X[i]);
            double constraintResponse ;

            // Scale to new distance
            for ( unsigned int j = 0 ; j < nVar ; j++ )
                constraint[j] *= scaledDist/barDist[i];

            // Scale response according to distance of constraint position from border point
            constraintResponse = -0.1 * (scaledDist-barDist[i])/(scaledDist-minBarDist);
            if (constraintResponse>=0.0)
                smartuq_throw("constraintResponse for normal constraint is bigger or equal than 0");

            // Add new constraint point and response
            X.push_back(constraint);
            Y.push_back(constraintResponse);
        }

    }

}



/**
 * Constructor for rbf class
 *
 * @param rbf_function : Type of radial function
 */
rbf::rbf ( const RBF_FUNCTION &rbf_function,
           const RBF_CONSTANTS &constants ) :
        base_surrogate(), m_rbf_function(rbf_function), m_constants(constants)
{
    if (m_rbf_function==MULTIQUADRATICS||m_rbf_function==GAUSSIAN)
        if (m_constants.avgdist==0.0)
            smartuq_throw("Average distance is explicitly required for the requested function basis");

    switch ( m_rbf_function ) {
        case LINEAR:
        case MULTIQUADRATICS :
            m_constants.smooth = +std::abs(m_constants.smooth);
            break;
        case GAUSSIAN :
        case CUBIC :
        case THIN_PLATE :
            m_constants.smooth = -std::abs(m_constants.smooth);
            break;
    }
}


/**
 * Constructor for rbf class
 *
 * @param X :            Set of points in which independent rbf are centered
 * @param Y :            Response of independent rbf in points X
 * @param rbf_function : Type of radial function
 */
rbf::rbf ( const std::vector<std::vector<double>> &X,
           const std::vector<double> &Y,
           const RBF_FUNCTION &rbf_function,
           const RBF_CONSTANTS &constants ) :
        base_surrogate(X,Y), m_rbf_function(rbf_function), m_constants(constants)
{
    if (m_rbf_function==MULTIQUADRATICS||m_rbf_function==GAUSSIAN)
        if (m_constants.avgdist==0.0)
            smartuq_throw("Average distance is explicitly required for the requested function basis");
}


//rbf::~rbf () {}

/**
 * Function to set basis points and responses
 *
 * @param X : Set of points in which independent rbf are centered
 * @param Y : Response of independent rbf in points X
 */
void rbf::set_points( const std::vector<std::vector<double>> &X,
                      const std::vector<double> &Y) {

    if (X.size()!=Y.size())
        smartuq_throw("Input set of points X and responses Y have mismatching dimension");

    m_X.clear();
    m_Y.clear();

    m_X = X;
    m_Y = Y;
}



/**
 * build virtual function implementation
 *
 * @brief Function that build the surrogate once a set of points and responses has been given
 */
void rbf::build () {

    unsigned int np = m_X.size();
    unsigned int n  = m_X[0].size();
    unsigned int ns = np+1+n;

    if ( np==0 )
        smartuq_throw("Basis center points not yet set");

    // Construct lower triangular part of symmetric distance matrix
    std::vector<std::vector<double>> r ;
    for ( unsigned int i = 0 ; i < np ; i++ ) {
        std::vector<double> r_xi;
        for (unsigned int j = 0; j < i; j++ ) {
            double norm = 0.0;
            for ( unsigned int k = 0 ; k < n ; k++ )
                norm += std::pow(m_X[i][k]-m_X[j][k],2.0);
            r_xi.push_back(std::sqrt(norm));
        }
        r_xi.push_back(0.0);
        r.push_back(r_xi);
    }


    // Evaluate lower triangular matrix accorgingly to selected basis
    evaluate_basis(r);

    // If smoothing is requested
    for ( unsigned int i = 0 ; i < np ; i++ )
        r[i][i] -= m_constants.smooth;

    // Construct evaluation matrix
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(ns,ns);
    for ( unsigned int i = 0 ; i < np ; i++ ) {
        for ( unsigned int j = 0 ; j < i ; j++ ) {
            A(i,j) = r[i][j];
            A(j,i) = A(i,j);
        }
        A(i,i) = r[i][i];
        A(i,np) = 1.0;
        A(np,i) = 1.0;

        unsigned int col = np+1;
        for ( unsigned int j = 0 ; j < n ; j++ ) {
            A(i,col)   = m_X[i][j];
            A(col++,i) = m_X[i][j];
        }
    }

    // Clear data
    r.clear();

    // Construct known term
    Eigen::VectorXd b = Eigen::VectorXd::Zero(ns);
    for ( unsigned int i = 0 ; i < np ; i++ )
        b(i) = m_Y[i];

    // Compute coefficients
    Eigen::VectorXd coe = A.fullPivLu().solve(b) ;

    // Insert coefficients
    m_coeff.clear();
    m_coeff.reserve(ns);

    for ( unsigned int i = 0 ; i < ns ; i++ )
        m_coeff.push_back(coe(i));

}


void rbf::evaluate_basis( std::vector<std::vector<double>> &R ) const{

    unsigned int np = R.size();

    for ( unsigned int i = 0 ; i < np ; i++ ) {
        evaluate_basis(R[i]);
    }

}

void rbf::evaluate_basis( std::vector<double> &r ) const {

    unsigned int nc = r.size();
    for ( unsigned int j = 0 ; j < nc ; j++ ) {
        //// TODO ////
        // This shall be implemented properly at a later stage
        switch (m_rbf_function) {
            case THIN_PLATE :
                r[j] = r[j] * r[j] * std::log(r[j] + 1.0);
                break;
            case LINEAR :
                break;
            case CUBIC :
                r[j] = std::pow(r[j], 3.0);
                break;
            case MULTIQUADRATICS :
                r[j] = std::sqrt(1 + r[j] * r[j] / (m_constants.avgdist * m_constants.avgdist));
                break;
            case GAUSSIAN :
                r[j] = std::exp(-0.5 * r[j] * r[j] / (m_constants.avgdist * m_constants.avgdist));
                break;
        }
    }

}


/**
 * evaluate virtual function implementation
 *
 * @brief Function that evaluate a new single point x with the constructed surrogate to compute the response y
 */
void rbf::evaluate (const std::vector<double> &x, double &y) const {


    if ( m_coeff.empty() )
        smartuq_throw("Coefficients not computed yet");

    unsigned int np = m_X.size();
    unsigned int n  = m_X[0].size();

    if ( x.size()!=n )
        smartuq_throw("Input point has not correct dimension");


    std::vector<double> r;
    for ( unsigned int i = 0 ; i < np ; i++ ) {
        double norm = 0.0;
        for (unsigned int j = 0; j < n; j++)
            norm += std::pow(x[j]-m_X[i][j],2.0);
        r.push_back(std::sqrt(norm));
    }
    evaluate_basis(r);

    y = 0.0;

    // Contribution of rbf functions
    for ( unsigned int i = 0 ; i < np ; i++ )
        y += m_coeff[i]*r[i];

    // Contribution of constant term
    y += m_coeff[np]*1.0;

    // Contribution of linear term
    unsigned int cc = np+1;
    for ( unsigned int i = 0 ; i < n ; i++, cc++ )
        y += m_coeff[cc]*x[i];



}

std::vector<double> rbf::get_coeff() const {

    return m_coeff;
}



// Standard constants for rbf
const RBF_CONSTANTS rbf::rbf_std_constants{0.0,0.0};
