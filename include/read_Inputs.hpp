#ifndef READ_INPUTS_HPP
#define READ_INPUTS_HPP

#include <iostream>
#include <sstream>
#include <stdio.h> 
#include <algorithm>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"

#define EPS = 1e-16;

using namespace Eigen;

//defining class and function for matrix indexing
//this is to mimic the matlab capability to select 
//specific columns and rows from a matrix
template<class ArgType, class RowIndexType, class ColIndexType>
class indexing_functor {
  const ArgType &m_arg;
  const RowIndexType &m_rowIndices;
  const ColIndexType &m_colIndices;
public:
  typedef Eigen::Matrix<typename ArgType::Scalar,
                 RowIndexType::SizeAtCompileTime,
                 ColIndexType::SizeAtCompileTime,
                 ArgType::Flags&RowMajorBit?RowMajor:ColMajor,
                 RowIndexType::MaxSizeAtCompileTime,
                 ColIndexType::MaxSizeAtCompileTime> MatrixType;
  indexing_functor(const ArgType& arg, const RowIndexType& row_indices, const ColIndexType& col_indices)
    : m_arg(arg), m_rowIndices(row_indices), m_colIndices(col_indices)
  {}
  const typename ArgType::Scalar& operator() (Eigen::Index row, Eigen::Index col) const {
    return m_arg(m_rowIndices[row], m_colIndices[col]);
  }
};

template <class ArgType, class RowIndexType, class ColIndexType>
CwiseNullaryOp<indexing_functor<ArgType,RowIndexType,ColIndexType>, typename indexing_functor<ArgType,RowIndexType,ColIndexType>::MatrixType>
indexing(const MatrixBase<ArgType>& arg, const RowIndexType& row_indices, const ColIndexType& col_indices)
{
  typedef indexing_functor<ArgType,RowIndexType,ColIndexType> Func;
  typedef typename Func::MatrixType MatrixType;
  return MatrixType::NullaryExpr(row_indices.size(), col_indices.size(), Func(arg.derived(), row_indices, col_indices));
}

//offset for reading fortran binary files
const int RECORD_DELIMITER_LENGTH = 4;

// Structure to be filled with information from cfg file
struct prob_settings 
{

    //-----Definition of general variables--------
    int Ns;                         //Number of snapshots
    int Ds;                         //Delta between snapshots
    int nstart;                     //starting snapshot number
    int init_tres;
    int init_imode;
    int ndim;
    bool direct_error;
    bool surf_res;
    std::string solver;
    std::string in_file;            //Input filename
    std::string out_file;           //Output filename (for reconstruction)
    std::string flag_dim;           //Problem dimension (2D/3D)
    std::string flag_prob;          //Type of problem (Vector, scalar)
    std::vector<int> Cols;          //Fields columns to porcess
    std::vector<double> pointsource;
    std::vector<int> Cols_coords;   //Columns with coordinates 
    std::vector<std::string> flag_method;        //Method to use
    std::string flag_wdb_be;        //flag write database basis ectraction( modes and coefficients)
    double Dt_cfd;                  //delta t used in CFD simulation
    std::vector<double> Dt_res;     //delta t used for residual evaluation
    std::vector<int> t_pos;          //Position vector for adaptive sampling
    double alpha;                   //angle of attack of the problem
    double beta;                    //angle of sideslip
    double Mach;                    //Mach number
    double Re;                      //Reynolds number
    double mu;                      //dynamic viscosity
    double T;                       //Temperature
    double P;                       //Pressure
    double Rho;                     //Density

    //------Parameters for POD-SPOD---------
    std::vector<int> Nf;                     //filter size for SPOD
    double En;                  //Energetic content desired for reconstruction
    
    double sigma;               //sigma value in case of Gaussian filter
    std::string flag_filter;    //SPOD filter type
    std::string flag_mean;      //Subtract mean (ON/OFF), flag
    std::string flag_bc;        /*Type of boundary condition
                                    for correlation matrix, flag*/

    int r;                      //user defined rank (can be used for POD/SPOD/DMD)
    //------Parameters for DMD-----------
    
                                //if r=0 SVHT is used
    std::string dmd_coef_flag;  //method for coefficients calculation
    
    //only for mrDMD
    int max_cycles;
    int max_levels;

    //only for hoDMD
    int d;                      //level of high order approximation

    //only for RDMD
    int r_RDMD;                 //Number of modes for recursive DMD

    //----Parameter for Reconstruction---------
    std::string flag_rec;           //Activate-deactivate field reconstruction 
    std::string flag_interp;        //Interpolation technique for rbf
    std::vector<double> t_rec;                   //times desired for reconstruction
    std::vector<double> t_res;                   //times desired for computing residuals
    double tol;                     //Error tolerance for adaptive reconstruction

};


struct plot3d_info
{
    int nblocks; //number of blocks
    std::vector<int> ni, nj, nk; // number of elements along the three spacial directions for each block
    float M, alpha, Re, T; //Mach, angle of attack, Reynolds, Time non-dimensional
};



// List of Keywords in config file
enum keywords 
            {
                ADAPT_POS,
                SOLVER, 
                NS, DS, EN, NF, SIGMA,
                RHO_FREE, P_FREE,
                ALPHA, BETA, MACH, REYNOLDS, TEMPERATURE, VISCOSITY,
                NSTART,
                NDIM,
                DT_CFD,
                DT_RES,
                FLAG_DIM,
                FLAG_PROB,
                FLAG_METHOD, 
                INPUT_FILE,
                OUTPUT_FILE, 
                COLS_COORDS,
                COLS_FIELDS,
                FLAG_MEAN,
                FLAG_BC,
                FLAG_FILTER,
                FLAG_WDB_BE,
                FLAG_REC,
                FLAG_INTERP,
                T_REC,
                T_RES,
                RANK,
                DMD_COEF_FLAG,
                RANK_RDMD,
                HO_D,
                MAX_LEVELS,
                MAX_CYCLES,
                TOL,
                DIRECT_ERROR,
                SURF_RESEVAL,
                INIT_TRES,
                INIT_MODE,
                POINT_SOURCE
            };


// Compare keywords with input string
keywords read_keyword_type( const std::string &key_string );


// Read config file and store info in prob_settings
void Read_cfg ( std::string filename, prob_settings &settings );

// Read and change su2 file
void Modify_su2_cfg ( std::string file_in, std::string file_out, prob_settings settings, int it1, int it2, double U_inf = 0.0 );

//Get Number of grid points
int N_gridpoints ( const std::string file_in );


// Get fields on the specified columns 
Eigen::MatrixXd read_col( std::string filename, int Nr, std::vector<int> Cols );

Eigen::MatrixXd read_colnew( std::string filename, int Nr, std::vector<int> Cols );


// read modes for Recursive DMD
Eigen::MatrixXd read_modes( std::string filename, int Nr, int r_RDMD );

// read coeffs for Recursive DMD
Eigen::MatrixXd read_coefs( std::string filename, int Ns, int r_RDMD );

// read Errors and Jaccard index for adaptive reconstruction
Eigen::MatrixXd read_err_j ( std::string filename, int Ns );

//read .q (plot3d) file for visualization info
plot3d_info read_plot3d_info (std::string filename);

// read .q (plot3d) file, fortran binary for flow field information
std::vector<Eigen::VectorXd> read_plot3d (std::string filename, plot3d_info Info);

// read modes
void ModeDB_Read ( std::string root_file_m, std::string root_file_c, std::vector<Eigen::MatrixXd> &Phi, std::vector<Eigen::MatrixXd> &Coefs, prob_settings settings );



#endif //READ_INPUTS_HPP
