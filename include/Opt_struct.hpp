#ifndef OPT_STRUCT_HPP
#define OPT_STRUCT_HPP

#ifdef __USE_PAGMO

#include <initializer_list>
#include <utility>

#include "read_Inputs.hpp"
#include "Extract_Basis.hpp"
#include "pagmo.hpp"

// Define the problem PaGMO-style
using namespace Eigen;


//POD on absolute positions
struct SPOD_Adapt_Samp {

    // Empty constructor
    SPOD_Adapt_Samp( ){ }

    SPOD_Adapt_Samp( std::vector< std::vector< double > > &bounds, Eigen::MatrixXd &sn_set, prob_settings &settings,
            int &Nf ) : problemBounds_(bounds), m_sn_set(sn_set), m_settings(settings), m_Nf(Nf) {
        }
    // Fitness:
    std::vector<double> fitness(const std::vector<double> &variables) const;
    // Boundaries of the problem
    std::pair<std::vector<double>, std::vector<double>> get_bounds() const;
    //! Serialization function for Pagmo compatibility

private:

    Eigen::MatrixXd m_sn_set;
    prob_settings m_settings;
    const std::vector< std::vector< double > > problemBounds_;
    int m_Nf;
};


//POD on relative positions
struct SPOD_Adapt_Samp_ {

    // Empty constructor
    SPOD_Adapt_Samp_( ){ }

    SPOD_Adapt_Samp_( std::vector< std::vector< double > > &bounds, Eigen::MatrixXd &sn_set, prob_settings &settings,
                     int &Nf ) : problemBounds_(bounds), m_sn_set(sn_set), m_settings(settings), m_Nf(Nf) {
    }
    // Fitness:
    std::vector<double> fitness(const std::vector<double> &variables) const;
    // Boundaries of the problem
    std::pair<std::vector<double>, std::vector<double>> get_bounds() const;

    //Number of equality and inequality constraints
    pagmo::vector_double::size_type get_nec() const {
        return 0;
    }

    pagmo::vector_double::size_type get_nic() const {
        return 1;
    }

private:

    Eigen::MatrixXd m_sn_set;
    prob_settings m_settings;
    const std::vector< std::vector< double > > problemBounds_;
    int m_Nf;
};


//DMD on absolute positions
struct DMD_Adapt_Samp {

    // Empty constructor
    DMD_Adapt_Samp( ){ }

    DMD_Adapt_Samp( std::vector< std::vector< double > > &bounds, Eigen::MatrixXd &sn_set, prob_settings &settings) :
                            problemBounds_(bounds), m_sn_set(sn_set), m_settings(settings) {
    }
    // Fitness:
    std::vector<double> fitness(const std::vector<double> &variables) const;
    // Boundaries of the problem
    std::pair<std::vector<double>, std::vector<double>> get_bounds() const;
    //! Serialization function for Pagmo compatibility

private:

    Eigen::MatrixXd m_sn_set;
    prob_settings m_settings;
    const std::vector< std::vector< double > > problemBounds_;

};


// Pagmo problem with integers for POD

struct SPOD_Adapt_Samp_Int {

    // Empty constructor
    SPOD_Adapt_Samp_Int( ){ }

    SPOD_Adapt_Samp_Int( std::vector< std::vector< double > > &bounds, Eigen::MatrixXd &sn_set, prob_settings &settings,
                      int &Nf ) : problemBounds_(bounds), m_sn_set(sn_set), m_settings(settings), m_Nf(Nf) {
    }

    std::vector<double> fitness(const std::vector<double> &variables) const; // Fitness
    std::pair<std::vector<double>, std::vector<double>> get_bounds() const; // Boundaries of the problem

    //Number of continuous and integer variables
    pagmo::vector_double::size_type get_ncx() const {
        return 0;
    }
    pagmo::vector_double::size_type get_nix() const {
        return problemBounds_[0].size();
    }

private:

    Eigen::MatrixXd m_sn_set;
    prob_settings m_settings;
    const std::vector< std::vector< double > > problemBounds_;
    int m_Nf;
};

struct SPOD_Adapt_Samp_Int_ {

    // Empty constructor
    SPOD_Adapt_Samp_Int_( ){ }

    SPOD_Adapt_Samp_Int_( std::vector< std::vector< double > > &bounds, Eigen::MatrixXd &sn_set, prob_settings &settings,
                         int &Nf ) : problemBounds_(bounds), m_sn_set(sn_set), m_settings(settings), m_Nf(Nf) {
    }

    std::vector<double> fitness(const std::vector<double> &variables) const; // Fitness
    std::pair<std::vector<double>, std::vector<double>> get_bounds() const; // Boundaries of the problem

    //Number of equality and inequality constraints
    pagmo::vector_double::size_type get_nec() const {
        return 0;
    }
    pagmo::vector_double::size_type get_nic() const {
        return 1;
    }
    //Number of continuous and integer variables
    pagmo::vector_double::size_type get_ncx() const {
        return 0;
    }
    pagmo::vector_double::size_type get_nix() const {
        return problemBounds_[0].size();
    }

private:

    Eigen::MatrixXd m_sn_set;
    prob_settings m_settings;
    const std::vector< std::vector< double > > problemBounds_;
    int m_Nf;
};


struct DMD_Best_Modes {

    // Empty constructor
    DMD_Best_Modes( ){ }

    DMD_Best_Modes( std::vector< std::vector< double > > &bounds, Eigen::MatrixXd &sn_set, Eigen::MatrixXcd &Phi,
            Eigen::MatrixXcd &Coefs, prob_settings &settings ) : problemBounds_(bounds), m_sn_set(sn_set),
            m_Phi(Phi), m_Coefs(Coefs), m_settings(settings) {
    }
    // Fitness:
    std::vector<double> fitness(const std::vector<double> &variables) const;
    // Boundaries of the problem
    std::pair<std::vector<double>, std::vector<double>> get_bounds() const;
    //! Serialization function for Pagmo compatibility

    //Number of equality and inequality constraints
    pagmo::vector_double::size_type get_nec() const {return 0;}
    pagmo::vector_double::size_type get_nic() const {return 0;}
    //Number of continuous and integer variables
    pagmo::vector_double::size_type get_ncx() const {return 0;}
    pagmo::vector_double::size_type get_nix() const {return problemBounds_[0].size();}

private:

    Eigen::MatrixXd m_sn_set;
    Eigen::MatrixXcd m_Phi;
    Eigen::MatrixXcd m_Coefs;
    prob_settings m_settings;
    const std::vector< std::vector< double > > problemBounds_;
};




struct Best_Modes {

    // Empty constructor
    Best_Modes( ){ }

    Best_Modes( std::vector< std::vector< double > > &bounds, Eigen::MatrixXd &sn_set, Eigen::MatrixXd &Phi,
                Eigen::MatrixXd &Coefs, prob_settings &settings ) : problemBounds_(bounds), m_sn_set(sn_set),
                                                                         m_Phi(Phi), m_Coefs(Coefs), m_settings(settings) {
    }
    // Fitness:
    std::vector<double> fitness(const std::vector<double> &variables) const;
    // Boundaries of the problem
    std::pair<std::vector<double>, std::vector<double>> get_bounds() const;
    //! Serialization function for Pagmo compatibility

    //Number of equality and inequality constraints
    pagmo::vector_double::size_type get_nec() const {return 0;}
    pagmo::vector_double::size_type get_nic() const {return 0;}
    //Number of continuous and integer variables
    pagmo::vector_double::size_type get_ncx() const {return 0;}
    pagmo::vector_double::size_type get_nix() const {return problemBounds_[0].size();}

private:

    Eigen::MatrixXd m_sn_set;
    Eigen::MatrixXd m_Phi;
    Eigen::MatrixXd m_Coefs;
    prob_settings m_settings;
    const std::vector< std::vector< double > > problemBounds_;
};



#endif

#endif //OPT_STRUCT_HPP