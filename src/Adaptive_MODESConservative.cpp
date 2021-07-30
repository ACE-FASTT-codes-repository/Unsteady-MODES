/*
CODE FOR ERROR EVALUATION OF CONSERVATIVE QUANTITIES USING DIFFERENT RBM TECHNIQUES
!!!!!!(restart files need to be numbered without gaps)!!!!
INPUT ARGUMENTS 
Config File RBM
*/

#include "Extract_Basis.hpp"
#include "read_Inputs.hpp"
#include "Generate_snset.hpp"
#include "Reconstruction.hpp"
#include "write_Outputs.hpp"


int main( int argc, char *argv[] )
{

    std::cout << "-----------MODES Adaptive Residual Evaluation starts-------------" << std::endl << std::endl;

    std::cout << "Initializing common variables ... " << std::endl << std::endl;
    prob_settings settings;
    std::string filecfg = argv[1];

    std::vector<Eigen::VectorXd> Err_RBM_Nm_time;
    std::vector<Eigen::VectorXd> ErrP_RBM_Nm_time;
    std::vector<Eigen::VectorXd> EN;

    Read_cfg( filecfg, settings );

    std::string root_outputfile;
    root_outputfile.assign ( settings.out_file, 0, settings.out_file.size() - 4);
    std::string root_inputfile;
    root_inputfile.assign ( settings.in_file, 0, settings.in_file.size() - 4);
    std::string input_format;
    input_format.assign ( settings.in_file, settings.in_file.size() - 3, 3);

    int s_Nf = 5;   //Number of values for the SPOD filter (POD included)
    std::vector<int> Nf(s_Nf);
    Nf[0] = 0;
    Nf[1] = std::ceil(settings.Ns/10.0);
    Nf[2] = std::ceil(settings.Ns/2.0);
    Nf[3] = std::ceil(2.0*settings.Ns/3.0);
    Nf[4] = settings.Ns;

    int nC = settings.Cols.size();
    double alpha = settings.alpha;
    double beta = settings.beta;
    if (settings.ndim == 2) beta = 0.0;

    Eigen::MatrixXd norm_sn_set = Eigen::MatrixXd::Zero(nC, settings.Ns*settings.Ds);
    
    std::stringstream buffer;
    buffer << std::setfill('0') << std::setw(5) << std::to_string(settings.nstart);
    std::string file_1 = root_inputfile + "_" + buffer.str() + "." + input_format;

    // Calculate number of grid points
    int Nr = N_gridpoints ( file_1 );
    std::cout << "Number of grid points : " << Nr << std::endl;

    std::cout << "Reading Coordinates ... \t ";
    Eigen::MatrixXd Coords = read_col( file_1, Nr, settings.Cols_coords );
    std::cout << "Done " << std::endl;

// How we upload the snapshot matrix in the most efficient way? 
// By now one big igen Matrix where each column has all the vectors of the conservative variables
    std::cout << "Storing snapshot Matrix ... \n ";
    Eigen::MatrixXd sn_set = generate_snap_matrix( Nr, settings);

    std::cout << "Initializing Vector of times ... " << std::endl; 
    std::vector<double> t_vec( settings.Ns );
    t_vec[0] = settings.nstart*settings.Dt_cfd;
    for ( int i = 1; i < settings.Ns; i++ )
        t_vec[i] = t_vec[i-1] + settings.Dt_cfd*settings.Ds;
        
    // std::vector<double> t_evaluate(2*settings.Ns-1);
    // t_evaluate[0] = settings.nstart*settings.Dt_cfd;
    // for ( int i = 1; i < t_evaluate.size(); i++)
    //     t_evaluate[i] = t_evaluate[i-1] + settings.Dt_cfd*(double)settings.Ds/2.0;
    std::vector<double> t_evaluate(settings.Ns*settings.Ds-1);

    t_evaluate[0] = settings.nstart*settings.Dt_cfd;
    for ( int i = 1; i < t_evaluate.size(); i++)
        t_evaluate[i] = t_evaluate[i-1] + settings.Dt_cfd;// + settings.Dt_cfd*(double)settings.Ds/2.0;

    std::cout << "Computing mean/Initial Condition of CFD solution ... " << std::endl;
    //Defining Initial condition
    Eigen::VectorXd mean = sn_set.rowwise().mean();
    Eigen::VectorXd Ic = Eigen::VectorXd::Zero(nC*Nr);

    std::string binary = "YES";
    Eigen::VectorXd svd_cum_sum(settings.Ns);

    if ( settings.flag_mean == "IC" ) {
        std::cout << "Subtracting Initial condition" << std::endl << std::endl;
        Ic = IC(sn_set,settings,nC,Nr);
    } else {
        std::cout << "Using data without subtracting any reference state" << std::endl << std::endl;
    }

    
    Eigen::MatrixXd sn_set_p(nC*Nr, settings.Ns);
    for ( int is = 0; is < settings.Ns; is ++)
    {
        sn_set_p.col(is) = sn_set.col(settings.Ds*is);
    }

    for ( int j = 0; j < nC; j++ )   
    {
        for ( int i = 0; i < settings.Ns*settings.Ds; i++ )
            norm_sn_set(j, i) = sn_set.middleRows(j*Nr,Nr).col(i).norm();
        
    }

    if ( settings.flag_mean == "YES" )
    {
            std::cout << "Subtracting mean from snapshots ... " << std::endl;
        for ( int it = 0; it < settings.Ns*settings.Ds; it++ )
            sn_set.col(it) -= mean;

        for ( int it = 0; it < settings.Ns; it++ )
            sn_set_p.col(it) -= mean;
    }

    if ( settings.flag_mean == "IC" )
    {
            std::cout << "Subtracting mean from snapshots ... " << std::endl;
        for ( int it = 0; it < settings.Ns*settings.Ds; it++ )
            sn_set.col(it) -= Ic;

        for ( int it = 0; it < settings.Ns; it++ )
            sn_set_p.col(it) -= Ic;
    }


//Defining common scope for POD-SPOD
    {

        Eigen::VectorXd lambda(settings.Ns);
        Eigen::VectorXd K_pc(settings.Ns);
        Eigen::MatrixXd eig_vec(settings.Ns, settings.Ns);
        int Nm;
        int N_notZero;
        //Check only for POD for now
        for ( int nfj = 0; nfj < Nf.size(); nfj++ )
        {
            std::cout << "Computing SPOD " << Nf[nfj] << " reconstruction for each conservative variable ... " << "\n";
            
            for ( int ncons = 0; ncons < nC; ncons ++ )
            {
                std::cout << "Processing conservative variable " << ncons << std::endl;
                Eigen::MatrixXd Phi = SPOD_basis( sn_set_p.middleRows(ncons*Nr,Nr),
                                        lambda, K_pc, eig_vec,
                                        Nf[nfj],
                                        settings.flag_bc, 
                                        settings.flag_filter,  
                                        settings.sigma);            
                N_notZero = Phi.cols();
                if ( settings.r == 0 ) Nm = Nmod(settings.En, K_pc);
                else Nm = std::min(settings.r, N_notZero);
                std::cout << "Number of modes used in reconstruction " << Nm << std::endl;
                std::vector<rbf> surr_coefs =  getSurrCoefs (t_vec, eig_vec, settings.flag_interp);                
                Eigen::MatrixXd coef_t(t_evaluate.size(), Nm);

                std::vector<double> tr(1);
                for ( int j = 0; j < t_evaluate.size(); j++ )
                {    
                    tr[0] = t_evaluate[j];
                    for ( int i = 0; i < Nm; i++ )
                        surr_coefs[i].evaluate(tr, coef_t(j,i));
                }
                
                Eigen::MatrixXd Sig = Eigen::MatrixXd::Zero(Nm, Nm);
                for ( int i = 0; i < Nm; i++ )
                    Sig(i,i) = std::sqrt(lambda(i));

                std::cout << "Computing error of interpolation and error from projection..." << "\t";

                Eigen::MatrixXd Err_SPOD_map = Eigen::MatrixXd::Zero( Nr, settings.Ns*settings.Ds-1 );
                Eigen::MatrixXd Err_PSPOD_map = Eigen::MatrixXd::Zero( Nr, settings.Ns*settings.Ds-1 );
                Eigen::VectorXd Err_SPOD_Nm_time = Eigen::VectorXd::Zero(t_evaluate.size());
                Eigen::VectorXd J_SPOD_Nm_time = Eigen::VectorXd::Zero(t_evaluate.size());
                
                Err_SPOD_map = sn_set.middleRows(ncons*Nr,Nr).leftCols(settings.Ns*settings.Ds-1) - Phi.leftCols(Nm)*Sig*coef_t.transpose();
                Eigen::MatrixXd PhiTPhi = Phi.leftCols(Nm).transpose()*Phi.leftCols(Nm);
                Eigen::MatrixXd dumCoefs = Phi.leftCols(Nm).transpose()*sn_set.middleRows(ncons*Nr,Nr).leftCols(settings.Ns*settings.Ds-1);

                Err_PSPOD_map = sn_set.middleRows(ncons*Nr,Nr).leftCols(settings.Ns*settings.Ds-1) - Phi.leftCols(Nm)*(PhiTPhi.inverse()*dumCoefs);

                for ( int i = 0; i < settings.Ns*settings.Ds-1; i++ )
                {

                    int count = 0;
                    for ( int j = 0; j < Nr; j++ )
                    {
                        Err_SPOD_Nm_time(i) += Err_SPOD_map(j,i)*Err_SPOD_map(j,i);
                        J_SPOD_Nm_time(i) += Err_PSPOD_map(j,i)*Err_PSPOD_map(j,i);
                    }
                    
                    // Err_SPOD_Nm_time(i) = std::sqrt(Err_SPOD_Nm_time(i))/norm_sn_set(ncons,i);
                    // J_SPOD_Nm_time(i) = std::sqrt(J_SPOD_Nm_time(i))/norm_sn_set(ncons,i);
                    Err_SPOD_Nm_time(i) = std::sqrt(Err_SPOD_Nm_time(i)/(double)Nr);
                    J_SPOD_Nm_time(i) = std::sqrt(J_SPOD_Nm_time(i)/(double)Nr);
                
                }     

                Err_RBM_Nm_time.push_back(Err_SPOD_Nm_time);
                ErrP_RBM_Nm_time.push_back(J_SPOD_Nm_time);
                EN.push_back(K_pc);

                std::cout << "Done" << std::endl;

            }

            std::ofstream errfile;
            std::string file_err_name = "Err_SPOD_" + std::to_string(nfj) + ".dat";
            errfile.open(file_err_name);

            for ( int nm = 0; nm < settings.Ns*settings.Ds-1; nm ++ )    
            {
                for( int j = 0; j < Err_RBM_Nm_time.size(); j++ ) 
                    errfile <<  std::setprecision(8) << Err_RBM_Nm_time[j](nm) << "\t";

                errfile << std::endl;

            }

            errfile.close();

            std::ofstream errp;
            std::string file_errp_name = "ErrP_SPOD_" + std::to_string(nfj) + ".dat";
            errp.open(file_errp_name);

            for ( int nm = 0; nm < settings.Ns*settings.Ds-1; nm ++ )    
            {
                for( int j = 0; j < Err_RBM_Nm_time.size(); j++ ) 
                    errp <<  std::setprecision(8) << ErrP_RBM_Nm_time[j](nm) << "\t";

                errp << std::endl;

            }

            errp.close();
            //Add mean or Initial condition if it is subtracted
            Err_RBM_Nm_time.clear();
            ErrP_RBM_Nm_time.clear();
            EN.clear();

        }
    
    }

//Take care only of this first part for now
//---------------------------------------------------------------------------

// //     for ( int nt = 0; nt < settings.Ns; nt++ )
// //         sn_set.col(nt) += mean;
    
// //     for ( int nt = 0; nt < settings.Ns-1; nt++ )
// //         sn_set_check.col(nt) += mean;

// Defining scope for DMD ( Rank=-1 preferable, Coeffs = OPT )
    {
        Eigen::VectorXd lambda_POD;
        Eigen::MatrixXd eig_vec_POD;
        Eigen::VectorXcd lambda_DMD;
        Eigen::MatrixXcd eig_vec_DMD;

        for ( int ncons = 0; ncons < nC; ncons ++ )
        {
            std::cout << "Processing conservative variable " << ncons << std::endl;     

            std::cout << "Extracting basis DMD using rank " << std::min(settings.r, settings.Ns-1) << "\t";        
            Eigen::MatrixXcd Phi;
            Eigen::VectorXcd alfa;

            if ( settings.r == 0 )
            {
                Phi = DMD_basis( sn_set_p.middleRows(ncons*Nr,Nr),
                                lambda_DMD,
                                eig_vec_DMD,
                                lambda_POD,
                                eig_vec_POD,
                                -1 );
            }
            else
            {
                Phi = DMD_basis( sn_set_p.middleRows(ncons*Nr,Nr),
                                lambda_DMD,
                                eig_vec_DMD,
                                lambda_POD,
                                eig_vec_POD,
                                settings.r );
            }

    //         int Nm = Phi.cols();
    //         std::cout << "Number of modes extracted : " << Nm << std::endl;

            Eigen::VectorXcd omega(Phi.cols());
            for ( int i = 0; i < Phi.cols(); i++ )
                omega(i) = std::log(lambda_DMD(i))/(settings.Dt_cfd*settings.Ds);

            // std::cout << "Calculating coefficients DMD ... " << "\t";            
            alfa = Calculate_Coefs_DMD_exact ( sn_set.middleRows(ncons*Nr,Nr).leftCols(settings.Ns-1),  
                                                                lambda_DMD, 
                                                                Phi );
            // std::cout << " Done! " << std::endl;
            
            // std::cout << "Reordering modes DMD ... " << "\t";
            Eigen::VectorXd En = Eigen::VectorXd::Zero(Phi.cols());
            double T = t_vec[t_vec.size()-1];

            for ( int i = 0 ; i < Phi.cols(); i ++ )
            {

                double alfa_i = alfa(i).imag();
                double alfa_r = alfa(i).real();
                double sigma = omega(i).real();
                En(i) = (alfa_r*alfa_r + alfa_i*alfa_i)*(std::exp(2.0*sigma*T) - 1.0)/(2.0*sigma);

            }

            dmd_sort( En, Phi, lambda_DMD, alfa);
            // std::cout << "Done" << std::endl;

            double sum = 0;
            Eigen::VectorXd K_pc = Eigen::VectorXd::Zero(settings.Ns);
            for (int i = 0; i < Phi.cols(); i++)
            {
                sum += En(i)/En.sum();
                K_pc(i) = sum;
            }
            int Nm;
            if ( settings.r == 0)
            {
                Nm = Nmod(settings.En, K_pc);
                std::cout << "Number of modes for the desired energetic content : " << Nm << std::endl;   
            }
            else
            {
                Nm = std::min(settings.r,settings.Ns-1);
                std::cout << "Number of modes (fixed) : " << Nm << std::endl;
            }
            

            Eigen::MatrixXcd V_and(lambda_DMD.size(), t_evaluate.size());      
            for ( int i = 0; i < lambda_DMD.size(); i++ )
            {
                for ( int j = 0; j < t_evaluate.size(); j++ )
                    V_and(i,j) = std::pow(lambda_DMD(i), (double)j/(double)settings.Ds);                                                                                         
            }        
            Eigen::MatrixXcd Psi = Eigen::MatrixXcd::Zero(alfa.size(), t_evaluate.size());
            for ( int i = 0; i < t_evaluate.size(); i++ )
                Psi.col(i) = alfa.cwiseProduct(V_and.col(i));

            Eigen::MatrixXcd D_dmd = Phi.leftCols(Nm)*Psi.topRows(Nm);
            Eigen::MatrixXd Err_DMD_map = Eigen::MatrixXd::Zero( Nr, settings.Ns*settings.Ds-1 );
            Eigen::MatrixXd Err_PDMD_map = Eigen::MatrixXd::Zero( Nr, settings.Ns*settings.Ds-1 );
            Eigen::VectorXd Err_DMD_Nm_time = Eigen::VectorXd::Zero(t_evaluate.size());
            Eigen::VectorXd J_DMD_Nm_time = Eigen::VectorXd::Zero(t_evaluate.size());
            
            Err_DMD_map = sn_set.middleRows(ncons*Nr,Nr).leftCols(settings.Ns*settings.Ds-1) - D_dmd.real();

            Eigen::MatrixXcd PhiTPhi = Phi.leftCols(Nm).transpose()*Phi.leftCols(Nm);
            Eigen::MatrixXcd dumCoefs = Phi.leftCols(Nm).transpose()*sn_set.middleRows(ncons*Nr,Nr).leftCols(settings.Ns*settings.Ds-1);
            Eigen::MatrixXcd P_u = Phi.leftCols(Nm)*(PhiTPhi.inverse()*dumCoefs);
            Err_PDMD_map = sn_set.middleRows(ncons*Nr,Nr).leftCols(settings.Ns*settings.Ds-1) - P_u.real();

            // Err_SPOD_map = sn_set.middleRows(ncons*Nr,Nr).leftCols(settings.Ns*settings.Ds-1) - Phi.leftCols(Nm)*Sig*coef_t.transpose();
            // Eigen::MatrixXd PhiTPhi = Phi.leftCols(Nm).transpose()*Phi.leftCols(Nm);
            // Eigen::MatrixXd dumCoefs = Phi.leftCols(Nm).transpose()*sn_set.middleRows(ncons*Nr,Nr).leftCols(settings.Ns*settings.Ds-1);

            // Err_PSPOD_map = sn_set.middleRows(ncons*Nr,Nr).leftCols(settings.Ns*settings.Ds-1) - Phi.leftCols(Nm)*(PhiTPhi.inverse()*dumCoefs);

            for ( int i = 0; i < settings.Ns*settings.Ds-1; i++ )
            {

                int count = 0;
                for ( int j = 0; j < Nr; j++ )
                {
                    Err_DMD_Nm_time(i) += Err_DMD_map(j,i)*Err_DMD_map(j,i);
                    J_DMD_Nm_time(i) += Err_PDMD_map(j,i)*Err_PDMD_map(j,i);
                }
                
                // Err_DMD_Nm_time(i) = std::sqrt(Err_DMD_Nm_time(i))/norm_sn_set(ncons,i);
                // J_DMD_Nm_time(i) = std::sqrt(J_DMD_Nm_time(i))/norm_sn_set(ncons,i);
                Err_DMD_Nm_time(i) = std::sqrt(Err_DMD_Nm_time(i)/(double)Nr);
                J_DMD_Nm_time(i) = std::sqrt(J_DMD_Nm_time(i)/(double)Nr);
            
            }     

            Err_RBM_Nm_time.push_back(Err_DMD_Nm_time);
            ErrP_RBM_Nm_time.push_back(J_DMD_Nm_time);
            EN.push_back(K_pc);

        }

        std::ofstream errfile;
        std::string file_err_name = "Err_DMD.dat";
        errfile.open(file_err_name);

        for ( int nm = 0; nm < settings.Ns*settings.Ds-1; nm ++ )    
        {
            for( int j = 0; j < Err_RBM_Nm_time.size(); j++ ) 
                errfile <<  std::setprecision(8) << Err_RBM_Nm_time[j](nm) << "\t";

            errfile << std::endl;

        }

        errfile.close();

        std::ofstream errp;
        std::string file_errp_name = "ErrP_DMD.dat";
        errp.open(file_errp_name);

        for ( int nm = 0; nm < settings.Ns*settings.Ds-1; nm ++ )    
        {
            for( int j = 0; j < Err_RBM_Nm_time.size(); j++ ) 
                errp <<  std::setprecision(8) << ErrP_RBM_Nm_time[j](nm) << "\t";

            errp << std::endl;

        }

        errp.close();

        Err_RBM_Nm_time.clear();
        ErrP_RBM_Nm_time.clear();
        EN.clear();
        
        std::cout << "Done" << std::endl;

    }

// // //Defining scope for RDMD
// // //if using the function RDMD_modes_coefs for energybased select 
// // //energy level and rank rdmd to zero, for mode based just select
// //rank rdmd to the number of desired modes
    {

        Eigen::VectorXd lambda = Eigen::VectorXd::Zero(settings.Ns);
        Eigen::VectorXd K_pc = Eigen::VectorXd::Zero(settings.Ns);
        Eigen::MatrixXd Coefs = Eigen::MatrixXd::Zero(settings.Ns, settings.Ns);
        Eigen::MatrixXd Phi;

        for ( int ncons = 0; ncons < nC; ncons++ )
        {
            std::cout << "Processing conservative variable " << ncons << std::endl;
            std::cout << "Extracting basis and Coeffs RDMD ... " << "\t";        
        
            Phi = RDMD_modes_coefs ( sn_set_p.middleRows(ncons*Nr,Nr),
                                    Coefs,
                                    lambda,
                                    K_pc,     
                                    -1, //Performing singular value hard threshold for DMD reduction at each step
                                    settings.r_RDMD,
                                    settings.En );
                
//             // else
//             // {
//             //     std::cout << "Reading basis and extracting Coeffs RDMD ... " << "\t";
//             //     std::string file_modes = argv[2];
//             //     std::string file_coefs = argv[3];
//             //     std::string file_En = argv[4];
//             //     Phi = read_modes( file_modes, settings.ndim*Nr, settings.Ns );
//             //     Coefs = read_coefs( file_coefs, settings.Ns, settings.Ns );


//             //     std::ifstream En_data;
//             //     En_data.open( file_En );
//             //     if ( !En_data.is_open() )
//             //     {
//             //         std::cout << "File : " << file_En << " not found" << std::endl;    
//             //         exit (EXIT_FAILURE);
//             //     }
//             //     std::string line_flow_data ;
//             //     getline( En_data, line_flow_data );
//             //     std::istringstream iss(line_flow_data);
//             //     std::string token;

//             //     int count = 0;
//             //     while( getline( iss, token, ' ') && count < K_pc.size() )
//             //     {
//             //         K_pc(count) = std::stod(token);
//             //         count ++;
//             //     } 
//             //     En_data.close();
//             // }

            int Nm;
            if ( settings.r == 0 )
            {
                Nm = Nmod(settings.En, K_pc);
                std::cout << "number of modes for the desired energetic content " << Nm << std::endl;
            }
            else
            {
                Nm = std::min(settings.r,settings.Ns-1);
                std::cout << "number of modes (fixed) " << Nm << std::endl;
            }

            std::vector<rbf> surr_coefs =  getSurrCoefs (t_vec,
                                                        Coefs.transpose(),
                                                        settings.flag_interp);
            
            Eigen::MatrixXd coef_t(t_evaluate.size(), Nm);

            std::vector<double> tr(1);
            for ( int j = 0; j < t_evaluate.size(); j++ )
            {    
                tr[0] = t_evaluate[j];
                for ( int i = 0; i < Nm; i++ )
                    surr_coefs[i].evaluate(tr, coef_t(j,i));
            }

            std::cout << "Computing error RDMD and error from projection ... " << std::endl << std::endl;
            Eigen::MatrixXd Err_RDMD_map;
            Eigen::MatrixXd Err_PRDMD_map;
            Eigen::VectorXd Err_RDMD_Nm_time = Eigen::VectorXd::Zero(t_evaluate.size());
            Eigen::VectorXd J_RDMD_Nm_time = Eigen::VectorXd::Zero(t_evaluate.size());

    
            // Err_RDMD_map = sn_set_check - Phi.leftCols(Nm)*coef_t.transpose();
            Err_RDMD_map = sn_set.middleRows(ncons*Nr,Nr).leftCols(settings.Ns*settings.Ds-1) - Phi.leftCols(Nm)*coef_t.transpose();
            Err_PRDMD_map = sn_set.middleRows(ncons*Nr,Nr).leftCols(settings.Ns*settings.Ds-1) - Phi.leftCols(Nm)*(Phi.leftCols(Nm).transpose()*sn_set.middleRows(ncons*Nr,Nr).leftCols(settings.Ns*settings.Ds-1));

            // for ( int i = 0; i < Err_RDMD_map.cols(); i++ )
            // {
            //     Err_RDMD_map.col(i) -= mean;
            //     Err_PRDMD_map.col(i) -= mean;
            // }

            for ( int i = 0; i < settings.Ns*settings.Ds-1; i++ )
            {
                int count = 0;

                for ( int j = 0; j < Nr; j++ )
                {
                    Err_RDMD_Nm_time(i) += Err_RDMD_map(j,i)*Err_RDMD_map(j,i);
                    J_RDMD_Nm_time(i) += Err_PRDMD_map(j,i)*Err_PRDMD_map(j,i);

                }

                // Err_RDMD_Nm_time(i) = std::sqrt(Err_RDMD_Nm_time(i))/norm_sn_set(ncons,i);
                // J_RDMD_Nm_time(i) = std::sqrt(J_RDMD_Nm_time(i))/norm_sn_set(ncons,i);
                Err_RDMD_Nm_time(i) = std::sqrt(Err_RDMD_Nm_time(i)/(double)Nr);
                J_RDMD_Nm_time(i) = std::sqrt(J_RDMD_Nm_time(i)/(double)Nr);
            }
    
            Err_RBM_Nm_time.push_back(Err_RDMD_Nm_time);
            ErrP_RBM_Nm_time.push_back(J_RDMD_Nm_time);
            EN.push_back(K_pc);

        }

        std::ofstream errfile;
        std::string file_err_name = "Err_RDMD.dat";
        errfile.open(file_err_name);

        for ( int nm = 0; nm < settings.Ns*settings.Ds-1; nm ++ )    
        {
            for( int j = 0; j < Err_RBM_Nm_time.size(); j++ ) 
                errfile <<  std::setprecision(8) << Err_RBM_Nm_time[j](nm) << "\t";

            errfile << std::endl;

        }

        errfile.close();

        std::ofstream errp;
        std::string file_errp_name = "ErrP_RDMD.dat";
        errp.open(file_errp_name);

        for ( int nm = 0; nm < settings.Ns*settings.Ds-1; nm ++ )    
        {
            for( int j = 0; j < Err_RBM_Nm_time.size(); j++ ) 
                errp <<  std::setprecision(8) << ErrP_RBM_Nm_time[j](nm) << "\t";

            errp << std::endl;

        }

        errp.close();

        Err_RBM_Nm_time.clear();
        ErrP_RBM_Nm_time.clear();
        EN.clear();

        std::cout << "Done" << std::endl;

    }

    std::cout << "MODES Adaptive Residual Evaluation ends" << std::endl << std::endl;
    return 0;

}