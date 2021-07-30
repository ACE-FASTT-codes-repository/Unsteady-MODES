#ifndef EXTRACT_BASIS_HPP
#define EXTRACT_BASIS_HPP


#include "ctime"
#include <complex> 
#include <cmath>
#include "read_Inputs.hpp"
#include "Pre-Process.hpp"

int Nmod ( double En, Eigen::VectorXd K_pc );


void eig_sort( Eigen::VectorXd &lam, Eigen::MatrixXd &eig_vec );


void dmd_sort( Eigen::VectorXd &En, 
                Eigen::MatrixXcd &Phi, 
                Eigen::VectorXcd &lam, 
                Eigen::VectorXcd &alfa );


void dmd_tenvelope_sort( Eigen::MatrixXcd &Phi,
                        Eigen::VectorXcd &omega,
                        Eigen::VectorXcd &alfa,
                        std::vector<double> t_vec);


int SVHT ( Eigen::VectorXd lam, int m, int n );


int Not_zero( Eigen::VectorXd lam );


struct node_mrDMD
{
    
    int l;                              //level number
    int bin_num;                        //time bin number
    int bin_size;                       //time bin size
    int step;                           //Delta between snaps based on nyquist frequency
    int start, stop;                    //starting and stopping index for snapshots
    double rho;                         //cut-off frequency
    double t_begin, t_end;              //initial and final instant of each window with respect to the whole interval
    double dt;
    int r;                              //rank reduction
    int n;                              //number of slow modes
    Eigen::MatrixXcd Modes;             //Matrix of modes
    Eigen::VectorXcd Coefs;             //Vector of optimized coefs
    Eigen::VectorXcd lam;               //Vector of eigenvalues
    Eigen::MatrixXcd Psi;               //Time evolution matrix

};

//Calculate the nullspace of a given set of snapshots
Eigen::MatrixXd nullspace( Eigen::VectorXd s, Eigen::MatrixXd vh, const double atol = 1e-13, const double rtol = 0 );


//Check if linear consistency is verified ( null(X) included in null(Y) with X = [x_0, x_1, ..., x_(N-1)], Y = [x_1, x_2, ..., x_N])
bool check_linear_consistency( Eigen::MatrixXd X, Eigen::MatrixXd Y, Eigen::MatrixXd Nullspace, const double atol = 1e-13, const double rtol = 0 );


std::string method_selected ( int n, int &Nf_SPOD, std::vector<int> Nf );

//SPOD/POD basis extraction (calculates also coefficients)
Eigen::MatrixXd SPOD_basis( const Eigen::MatrixXd &snap_set,
                                Eigen::VectorXd &lam,
                                Eigen::VectorXd &K_pc,
                                Eigen::MatrixXd &eig_vec,
                                const int Nf = 0,
                                std::string bc_flag = "ZERO", 
                                std::string filter_flag = "BOX",  
                                double sigma = 1.0);

//DMD basis extraction
Eigen::MatrixXcd DMD_basis( const Eigen::MatrixXd &snap_set,
                            Eigen::VectorXcd &lam,
                            Eigen::MatrixXcd &eig_vec,
                            Eigen::VectorXd &lam_POD,
                            Eigen::MatrixXd &eig_vec_POD,
                            const int r = 0 );

//Calculating optimized Coefficients according to Schmid and Jovanovic
Eigen::VectorXcd Calculate_Coefs_DMD ( const Eigen::MatrixXcd &eig_vec,
                                    const Eigen::MatrixXcd &eig_vec_POD,
                                    const Eigen::VectorXcd &lam,
                                    const Eigen::VectorXcd &lam_POD,
                                    const int Ns );

//Calculating optimized Coefficients following the implementation on the Weblog of DMD
Eigen::VectorXcd Calculate_Coefs_DMD_exact ( const Eigen::MatrixXd &sn_set,  //matrix of first Ns-1 snaps 
                                            const Eigen::VectorXcd &lam,  //slow eigenvalues
                                            const Eigen::MatrixXcd &Phi ); //slow exact DMD modes

//Calculate DMD coefficients solving a least square problem for each time step
Eigen::MatrixXcd Calculate_Coefs_Matrix_DMD ( const Eigen::MatrixXd &sn_set,
                                            const Eigen::MatrixXcd &Phi,
                                            const Eigen::VectorXcd &omega,
                                            const double t_0,
                                            const double dt_dmd);

//Multi-resolution DMD
std::vector<node_mrDMD> mrDMD_basis( Eigen::MatrixXd &snap_set,       //Initial set of snapshots
                                    std::vector<node_mrDMD> &nodes,          //Initialize as empty vector
                                    const int r,                                  //DMD rank
                                    double dts,                             //dt between initial snapshots
                                    double t_0 = 0.0,                       //time instant of the first snapshot of the series
                                    int level = 0,                          
                                    int bin_num = 0,
                                    int offset = 0,
                                    int max_levels = 7,
                                    int max_cycles = 2,
                                    std::string flag_coefs = "OPT" );


//Recursive DMD according to Noack (calculates also the matrix of coefficients)
Eigen::MatrixXd RDMD_modes_coefs ( const Eigen::MatrixXd &sn_set,
                                    Eigen::MatrixXd &Coefs,     //Define N_mod RDMD through the dimension of matrix Coefs
                                    Eigen::VectorXd &lambda,
                                    Eigen::VectorXd &K_pc,
                                    const int r,                //rank of pure DMD at each step
                                    int &rdmd,
                                    double En );


//Recursive DMD with least square method to compute coefficients
Eigen::MatrixXd RDMD_lsq_basis ( const Eigen::MatrixXd &sn_set,
                                 Eigen::MatrixXd &Coefs,
                                 Eigen::VectorXd &lambda,
                                 Eigen::VectorXd &K_pc,
                                 const int r,
                                 int &rdmd,
                                 double En );


//Forward-Backward DMD
Eigen::MatrixXcd fbDMD_basis ( const Eigen::MatrixXd &snap_set,
                            Eigen::VectorXcd &lam,
                            Eigen::MatrixXcd &eig_vec,
                            const int r );


//High order DMD
Eigen::MatrixXcd HODMD_basis( const Eigen::MatrixXd &snap_set,  //Initial set of snapshots
                            Eigen::VectorXcd &lam,               //High order Eigenvalues 
                            Eigen::MatrixXcd &eig_vec,           //High order eigenvectors
                            Eigen::VectorXcd &Coefs,
                            const double tol,                   //tolerance for svd (see Soledad for definition)
                            const int d );                       //levels of high order DMD



//SPOD/POD basis extraction (calculates also coefficients)
Eigen::MatrixXd GPOD_basis( const double Dt,
                                const Eigen::MatrixXd &snap_set,
                                Eigen::VectorXd &lam,
                                Eigen::MatrixXd &Coeffs,
                                const int &r);


//DMD basis extraction with non-uniform sampling
Eigen::MatrixXcd DMD_Adaptive_basis ( const Eigen::MatrixXd &snap_set,
                                      Eigen::VectorXcd &lam,
                                      Eigen::MatrixXcd &eig_vec,
                                      Eigen::VectorXd &lam_POD,
                                      Eigen::MatrixXd &eig_vec_POD,
                                      Eigen::VectorXi &tpos );


Eigen::MatrixXd RDMD_Adaptive_basis ( const Eigen::MatrixXd &sn_set,
                                   Eigen::MatrixXd &Coefs,     //Define N_mod RDMD through the dimension of matrix Coefs
                                   Eigen::VectorXd &K_pc,
                                   Eigen::VectorXi &tpos);

#endif //EXTRACT_BASIS_HPP