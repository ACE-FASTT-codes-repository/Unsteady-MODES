/*
Code for adaptive reconstruction
Input config file + error file (+ Modes,Coefs and Encontent RDMD if already available)

Output reconstructed field at the desired time instants with the adaptive technique
*/

#include "Extract_Basis.hpp"
#include "read_Inputs.hpp"
#include "Generate_snset.hpp"
#include "Reconstruction.hpp"
#include "write_Outputs.hpp"


int main( int argc, char *argv[] )
{
    std::cout << "Adaptive Reconstruction RBM-Clyde starts " << std::endl;

    prob_settings settings;
    std::string filecfg = argv[1];
    
    //Reading configuration file
    Read_cfg( filecfg, settings );
    double t_0 = settings.nstart*settings.Dt_cfd;

    int s_Nf = 5;
    int Nmethods = s_Nf + 2;
    std::vector<int> Nf(s_Nf);
    Nf[0] = 0;
    Nf[1] = std::ceil(settings.Ns/10.0);
    Nf[2] = std::ceil(settings.Ns/2.0);
    Nf[3] = std::ceil(2.0*settings.Ns/3.0);
    Nf[4] = settings.Ns;
    int Nf_SPOD = 0;

    std::string root_inputfile;
    root_inputfile.assign ( settings.in_file, 0, settings.in_file.size() - 4);
    std::string input_format;
    input_format.assign ( settings.in_file, settings.in_file.size() - 3, 3);

    std::stringstream buffer;
    buffer << std::setfill('0') << std::setw(5) << std::to_string(settings.nstart);
    std::string file_1 = root_inputfile + "_" + buffer.str() + "." + input_format;

    int Nr = N_gridpoints ( file_1 );
    std::cout << "Number of grid points : " << Nr << std::endl;

    std::cout << "Reading Coordinates ... \t ";
    Eigen::MatrixXd Coords = read_col( file_1, Nr, settings.Cols_coords );
    std::cout << "Done " << std::endl;

    std::cout << "Storing snapshot Matrix ... \n ";
    Eigen::MatrixXd sn_set = generate_snap_matrix( Nr, settings);

    Eigen::VectorXd mean = sn_set.rowwise().mean();
    Eigen::MatrixXd sn_set_p = Eigen::MatrixXd::Zero(sn_set.rows(), sn_set.cols());


    if ( settings.flag_mean == "IC" ) {
        std::cout << "Adaptive Direct Error can only subtract the mean as reference solution\n Exiting..." << std::endl << std::endl;
        exit(EXIT_FAILURE);

    } else if (settings.flag_mean == "YES") {
        std::cout << "Subtracting mean from snapshots ... " << std::endl;

        for ( int nt = 0; nt < settings.Ns; nt++ )
            sn_set.col(nt) -= mean;

    }


    for ( int nt = 0; nt < settings.Ns; nt++ )
        sn_set_p.col(nt) = sn_set.col(nt) - mean;

    std::cout << "Reading error information ... " << std::endl;

    std::string filename = argv[2];

    std::ifstream file_data;
    file_data.open( filename );
        if ( !file_data.is_open() )
    {
        std::cout << "File : " << filename << " not found" << std::endl;    
        exit (EXIT_FAILURE);
    }


    std::string line_flow_data ;
    int n_row = 0, count = 0;
    Eigen::MatrixXd Err_RBM = Eigen::MatrixXd::Zero(settings.Ns - 1, Nmethods); 
    while ( getline( file_data, line_flow_data ) )
    {
        
        std::istringstream iss(line_flow_data);
        std::string token;
        double err;
        count = 0; 
        while( getline( iss, token, '\t') )
        {
            err = std::stod(token);
            Err_RBM(n_row, count) = err;
            count ++;
        } 

        n_row++;
    }

    file_data.close();

//Adaptive reconstruction on each selected time step
    int best_method_idx;

    std::cout << "Initializing Vector of time ... " << std::endl; 
    Eigen::VectorXd t_vec( settings.Ns - 1);
    t_vec(0) = settings.nstart*settings.Dt_cfd + 0.5*settings.Dt_cfd*settings.Ds;
    for ( int i = 1; i < settings.Ns-1; i++ )
        t_vec(i) = t_vec(i-1) + settings.Dt_cfd*settings.Ds;
    
    double tol = settings.tol;
    int index1, index2;
    Eigen::VectorXd Err_interp(Nmethods);

    for ( int i = 0; i < settings.t_rec.size(); i++ )
    {
        std::vector<int> pos = {};
        std::cout << " Adaptive reconstruction at time : " << settings.t_rec[i] << std::endl;

        index1 = 0;
        index2 = 0;
        for ( int nt = 0; nt < t_vec.size()-1; nt ++ )
        {
            if ( (settings.t_rec[i] >= t_vec(nt)) && (settings.t_rec[i] <= t_vec(nt+1)) )
            {
                index1 = nt;
                index2 = nt+1;
                break;
            }
        }

        if ( index1 == index2 )
        {
            std::cout << "Time for reconstruction out of interval!" << std::endl;
            continue;
        }

        int count = 0;
        for ( int k = 0; k < Err_RBM.cols(); k ++ )
        {
            Err_interp(k) = Err_RBM(index1,k) + (Err_RBM(index2,k) - Err_RBM(index1,k))/
                                (settings.Dt_cfd*settings.Ds)*(settings.t_rec[i] - t_vec[index1]);
        }

        double eps = Err_interp.minCoeff( &best_method_idx );


        std::string method = method_selected ( best_method_idx, Nf_SPOD, Nf );
        std::cout << "Best method is " << method << " and Nf ( value meaningful only for SPOD ) : " << Nf_SPOD << std::endl;
        
        std::cout << " Error : " << Err_interp(best_method_idx) << std::endl;
                        
            
        std::cout << "Computing Reconstruction using selected method " << std::endl;
        
        if ( method == "SPOD" )
        {
            Eigen::VectorXd lambda(settings.Ns);
            Eigen::VectorXd K_pc(settings.Ns);
            Eigen::MatrixXd eig_vec(settings.Ns, settings.Ns);        

            Eigen::MatrixXd Phi = SPOD_basis( sn_set,
                                    lambda, K_pc, eig_vec,
                                    Nf_SPOD,
                                    settings.flag_bc, 
                                    settings.flag_filter,  
                                    settings.sigma);

            int Nm;

            if ( settings.r == 0 )
            {
                Nm = Nmod(settings.En, K_pc);
                std::cout << "Number of modes for desired energetic content: " << Nm << std::endl;
            }
            else
            {
                Nm = settings.r;
                std::cout << "Number of modes (fixed): " << Nm << std::endl;
            }

            std::vector<double> t_v( settings.Ns );
            t_v[0] = settings.nstart*settings.Dt_cfd;

            for ( int kt = 1; kt < settings.Ns; kt++ )
                t_v[kt] = t_v[kt-1] + settings.Dt_cfd*settings.Ds;

            Eigen::MatrixXd Rec = Reconstruction_S_POD ( t_v,
                                K_pc, lambda, eig_vec.transpose(),
                                Phi, settings.t_rec[i],
                                Nm,
                                settings.flag_prob,
                                settings.flag_interp ) ;

            if (settings.flag_mean == "YES") {
                for (int kt = 0; kt < Rec.cols(); kt++)
                    Rec.col(kt) = Rec.col(kt) + mean.segment(kt * Nr, Nr);
            }
            std::cout << "Writing reconstructed field ..." << "\t";

            write_Reconstructed_fields ( Rec, Coords,
                                    settings.out_file,
                                    settings.flag_prob, i );

            std::cout << "Done" << std::endl << std::endl;
            
        }


        if ( method == "DMD" )
        {

            Eigen::VectorXd lambda_POD;
            Eigen::MatrixXd eig_vec_POD;
            Eigen::VectorXcd lambda_DMD;
            Eigen::MatrixXcd eig_vec_DMD;      
            Eigen::MatrixXcd Phi;
            Eigen::VectorXcd alfa;    

            if ( settings.r == 0 )
            {
                Phi = DMD_basis( sn_set,
                                lambda_DMD,
                                eig_vec_DMD,
                                lambda_POD,
                                eig_vec_POD,
                                -1 );
            }
            else
            {
                Phi = DMD_basis( sn_set,
                                lambda_DMD,
                                eig_vec_DMD,
                                lambda_POD,
                                eig_vec_POD,
                                settings.r );
            }

            alfa = Calculate_Coefs_DMD_exact ( sn_set.leftCols(settings.Ns-1),  
                                                                lambda_DMD,  
                                                                Phi );                    

            std::cout << "Reordering modes DMD ... " << "\t";
            Eigen::VectorXd En = Eigen::VectorXd::Zero(Phi.cols());
            double T = t_vec[t_vec.size()-1] + 0.5*settings.Dt_cfd*settings.Ds;

            Eigen::VectorXcd omega(Phi.cols());
            for ( int i = 0; i < Phi.cols(); i++ )
                omega(i) = std::log(lambda_DMD(i))/(settings.Dt_cfd*settings.Ds);


            for ( int i = 0 ; i < Phi.cols(); i ++ )
            {

                double alfa_i = alfa(i).imag();
                double alfa_r = alfa(i).real();
                double sigma = omega(i).real();
                En(i) = (alfa_r*alfa_r + alfa_i*alfa_i)*(std::exp(2.0*sigma*T) - 1.0)/(2.0*sigma);

            }

            dmd_sort( En, Phi, lambda_DMD, alfa);
            std::cout << "Done" << std::endl;

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
                Nm = settings.r;
                std::cout << "Number of modes (fixed) : " << Nm << std::endl;
            }
        
            double t_dmd = settings.t_rec[i] - t_vec[0];
            Eigen::MatrixXcd Rec = Reconstruction_DMD ( t_dmd,
                                                    settings.Dt_cfd*settings.Ds,
                                                    alfa.topRows(Nm),
                                                    Phi.leftCols(Nm),
                                                    lambda_DMD.head(Nm),
                                                    settings.flag_prob );

            if (settings.flag_mean == "YES") {
                for (int kt = 0; kt < Rec.cols(); kt++)
                    Rec.real().col(kt) = Rec.real().col(kt) + mean.segment(kt * Nr, Nr);
            }
            std::cout << "Writing reconstructed field ..." << "\t";

            write_Reconstructed_fields ( Rec.real(), Coords,
                                    settings.out_file,
                                    settings.flag_prob, i );

            std::cout << "Done" << std::endl << std::endl;
        
        }


        if ( method == "RDMD" )
        {
        
            Eigen::VectorXd lambda = Eigen::VectorXd::Zero(settings.Ns);
            Eigen::VectorXd K_pc = Eigen::VectorXd::Zero(settings.Ns);
            Eigen::MatrixXd Coefs = Eigen::MatrixXd::Zero(settings.Ns, settings.Ns);
            Eigen::MatrixXd Phi;

            if ( argc == 3 )
            {
                std::cout << "Extracting basis and Coeffs RDMD ... " << "\t";        
                //You can define rank DMD at each time step from the config file ( use -1 for the adaptive study adviced)
                Phi = RDMD_modes_coefs ( sn_set_p,
                                        Coefs,
                                        lambda,
                                        K_pc,     
                                        settings.r,
                                        settings.r_RDMD,
                                        settings.En );
            }
            else
            {
                std::cout << "Reading basis and extracting Coeffs RDMD ... " << "\t";
                std::string file_modes = argv[3];
                std::string file_coefs = argv[4];
                std::string file_En = argv[5];
                Phi = read_modes( file_modes, settings.ndim*Nr, settings.Ns );
                Coefs = read_coefs( file_coefs, settings.Ns, settings.Ns );

                std::ifstream En_data;
                En_data.open( file_En );
                if ( !En_data.is_open() )
                {
                    std::cout << "File : " << file_En << " not found" << std::endl;    
                    exit (EXIT_FAILURE);
                }
                std::string line_flow_data ;
                getline( En_data, line_flow_data );
                std::istringstream iss(line_flow_data);
                std::string token;

                int count = 0;
                while( getline( iss, token, ' ') && count < K_pc.size() )
                {
                    K_pc(count) = std::stod(token);
                    count ++;
                } 
                En_data.close();

            }


            int Nm;
            if ( settings.r == 0 )
            {
                Nm = Nmod(settings.En, K_pc);
                std::cout << "number of modes for the desired energetic content " << Nm << std::endl;
            }
            else
            {
                Nm = settings.r;
                std::cout << "number of modes (fixed) " << Nm << std::endl;
            }


            std::vector<double> t_st_vec(settings.Ns);
            t_st_vec[0] = t_0;

            for ( int i = 1; i < settings.Ns; i++ )
                t_st_vec[i] = t_st_vec[i-1] + settings.Dt_cfd*settings.Ds;


            Eigen::MatrixXd Rec = Reconstruction_RDMD ( settings.t_rec[i],
                                                        t_st_vec,
                                                        Coefs.topRows(Nm),
                                                        Phi.leftCols(Nm),
                                                        settings.flag_prob,
                                                        settings.flag_interp );
            if (settings.flag_mean == "YES") {
                for (int kt = 0; kt < Rec.cols(); kt++)
                    Rec.col(kt) = Rec.col(kt) + mean.segment(kt * Nr, Nr);
            }
            std::cout << "Writing reconstructed field ..." << "\t";

            write_Reconstructed_fields ( Rec, Coords,
                                    settings.out_file,
                                    settings.flag_prob, i );

            std::cout << "Done" << std::endl << std::endl << std::endl;            

        }

    }

    std::cout << "Adaptive Reconstruction RBM-Clyde ends " << std::endl;

    return 0;
}