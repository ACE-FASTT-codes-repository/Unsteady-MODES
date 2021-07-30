#ifndef RECONSTRUCTION_HPP
#define RECONSTRUCTION_HPP

#include "Extract_Basis.hpp"
#include "Surrogates/rbf.h"

using namespace smartuq::surrogate;


smartuq::surrogate::RBF_FUNCTION get_key_rbf ( const std::string &key_string ); 


Eigen::MatrixXd Reconstruction_S_POD ( const std::vector<double> &t_vec,
                                const Eigen::VectorXd &K_pc,
                                const Eigen::VectorXd &lam,
                                const Eigen::MatrixXd &Coeffs,
                                const Eigen::MatrixXd &phi,
                                const double time,
                                // const double En,
                                int Nrec,
                                std::string flag_prob = "VELOCITY-2D",
                                std::string flag_interp = "LINEAR");


Eigen::MatrixXcd Reconstruction_DMD ( const double time, const double dt, 
                                    const Eigen::VectorXcd &alfa,
                                    const Eigen::MatrixXcd &Phi,
                                    const Eigen::VectorXcd &lam,
                                    const std::string flag_prob );


Eigen::MatrixXd TimeEvo_SPOD ( const std::vector<double> &t_vec,    //Vector of sampling times
                            const Eigen::VectorXd &time,            //Vector of reconstructed times
                            const Eigen::MatrixXd &Coefs,           //Matrix of POD coefficients (a(t) is along a column), dimension Ns x Nmodes
                            const Eigen::MatrixXd &Phi,             //Matrix of modes, dimension Npoints x Nmodes
                            const Eigen::VectorXd &lam,             //Lambda POD                        
                            const std::string &flag_interp);         //Flag to define interpolation methods for POD coefficients


Eigen::MatrixXcd TimeEvo_DMD ( Eigen::VectorXd &time,
                            double dt,
                            const Eigen::VectorXcd &alfa,
                            const Eigen::MatrixXcd &Phi,
                            const Eigen::VectorXcd &lam );
 

 Eigen::MatrixXd TimeEvo_RDMD ( const std::vector<double> &t_vec, //definitions are almost the same as in TimeEvo_SPOD
                            const Eigen::VectorXd &time,
                            const Eigen::MatrixXd &Coefs,
                            const Eigen::MatrixXd &Phi,
                            const std::string &flag_interp);

void nodes_mrDMD_sort( std::vector<node_mrDMD> &nodes ); //Sorting nodes in ordered levels


Eigen::MatrixXcd Reconstruction_mrDMD ( const double time,                      //time desired for reconstruction                                                                                                                                                                                                                                                                                                                              
                                    const double dts,                           //time between initial set of snapshots
                                    const std::vector<node_mrDMD> &nodes,       //nodes
                                    const std::string flag_prob );              


std::vector<rbf> getSurrCoefs ( const std::vector<double> &t_vec,
                            const Eigen::MatrixXd &Coeffs,
                            std::string flag_interp = "LINEAR" );


Eigen::MatrixXcd Reconstruction_Hybrid_DMD ( const double time,
                                        const std::vector<double> t_vec,
                                        const Eigen::MatrixXcd &alfa,
                                        const Eigen::MatrixXcd &Phi,
                                        const Eigen::VectorXcd &omega,
                                        const std::string flag_prob = "VELOCITY-2D",
                                        const std::string flag_interp = "LINEAR" );


Eigen::MatrixXd Reconstruction_RDMD ( const double time,
                                    const std::vector<double> t_vec,
                                    const Eigen::MatrixXd &alfa,
                                    const Eigen::MatrixXd &Phi,
                                    const std::string flag_prob,
                                    const std::string flag_interp );


Eigen::MatrixXd Reconstruction_DMD_Interp ( const double time,
                                            const std::vector<double> t_vec,
                                            const Eigen::MatrixXcd &alfa,
                                            const Eigen::MatrixXcd &Phi,
                                            const std::string flag_prob = "VELOCITY-2D",
                                            const std::string flag_interp = "LINEAR" );


#endif // RECONSTRUCTION_HPP

