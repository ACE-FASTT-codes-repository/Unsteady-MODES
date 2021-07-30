#include "Extract_Basis.hpp"
#include "read_Inputs.hpp"
#include "Generate_snset.hpp"
#include "Reconstruction.hpp"
#include "write_Outputs.hpp"

int main(int argc, char *argv[]) {

    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "-----------Single-MODES start-------------" << std::endl << std::endl;

    std::string filecfg = argv[1];
    prob_settings settings;

    //Reading configuration file
    Read_cfg( filecfg, settings );
    Config_stream ( settings );

    int Nr;
    Eigen::MatrixXd Coords;
    std::vector<double> t_vec(settings.Ns);
    common_vars( Nr, Coords, t_vec, settings);

    // Create matrix of snapshots
    std::cout << "Storing snapshot Matrix ... \n ";
    Eigen::MatrixXd sn_set = generate_snap_matrix( Nr, settings);

    Eigen::VectorXd mean = sn_set.rowwise().mean();
//    Eigen::VectorXd Ic = IC(settings, ,Nr)

    if ( settings.flag_mean == "YES" ) {
        std::cout << "Subtracting mean from snapshots ... " << std::endl << std::endl;
        for ( int i = 0; i < settings.Ns; i++ )
            sn_set.col(i) -= mean;
    } else if ( settings.flag_mean == "IC") {
        std::cout << "Mean is the only reference solution implemented in SingleMODES\n Exiting ... " << std::endl << std::endl;
        exit(EXIT_FAILURE);
    }

    Eigen::VectorXd K_pc(settings.Ns);

    auto methods = settings.flag_method;
    //Defining common scope for POD-SPOD
    std::vector<std::string>::iterator itSPOD;
    itSPOD = std::find (methods.begin(), methods.end(), "SPOD");
    if (itSPOD != methods.end()) {

        std::vector<double> t_vec( settings.Ns );
        t_vec[0] = 0.0;

        for ( int i = 1; i < settings.Ns; i++ )
            t_vec[i] = t_vec[i-1] + settings.Dt_cfd*settings.Ds;

        std::cout << std::endl;
        std::cout << "Initialized vector of times " << std::endl;

        Eigen::VectorXd lambda(settings.Ns);
        Eigen::MatrixXd eig_vec(settings.Ns, settings.Ns);

        std::cout << "Extracting basis ... " << "\t";        

        Eigen::MatrixXd Phi = SPOD_basis( sn_set,
                                lambda, K_pc, eig_vec,
                                settings.Nf[0],
                                settings.flag_bc, 
                                settings.flag_filter,  
                                settings.sigma);

        std::cout << " Done! " << std::endl << std::endl;

        std::cout << "Number of non-zero modes : " << Phi.cols() << std::endl;

        int Nrec; 
        if ( settings.r == 0 ) {
            Nrec = Nmod( settings.En, K_pc);
            std::cout << "Number of modes for the desired energy content : " << Nrec << std::endl;
        } else {
            int Nmax = Phi.cols();
            Nrec = std::min(settings.r,Nmax);
            std::cout << " Number of modes fixed : " << Nrec << std::endl;
        }


        if ( settings.flag_wdb_be == "YES" ) {
            
            std::cout << "Writing modes ..." << "\t";
            write_modes_sPOD ( Phi.leftCols(Nrec), Coords, settings.flag_prob );
            std::cout << "Complete!" << std::endl;

            std::cout << "Writing Coefficients ..." << "\t";
            write_coeffs_sPOD ( eig_vec, t_vec, lambda );
            std::cout << "Complete!" << std::endl;
            std::cout << std::endl;
        }

        if ( settings.flag_rec == "YES" ) {
            for ( int nt = 0; nt < settings.t_rec.size(); nt++ ) {
                std::cout << "Reconstructing field at time : " << settings.t_rec[nt] << "\t";
                Eigen::MatrixXd Rec = Reconstruction_S_POD ( t_vec,
                                    K_pc, lambda, eig_vec.transpose(),
                                    Phi, settings.t_rec[nt],
                                    Nrec,
                                    settings.flag_prob,
                                    settings.flag_interp ) ;

                std::cout << "Done" << std::endl;

                if ( settings.flag_mean == "YES" ) {
                    for ( int i = 0; i < Rec.cols(); i++)
                        Rec.col(i) = Rec.col(i) + mean.segment(i*Nr, Nr);

                }

                std::cout << "Writing reconstructed field ..." << "\t";

                write_Reconstructed_fields ( Rec, Coords,
                                        settings.out_file,
                                        settings.flag_prob, nt );

                std::cout << "Done" << std::endl << std::endl;
            }
        }
    }

    //Defining common scope for POD-SPOD
    std::vector<std::string>::iterator itDMD;
    itDMD = std::find (methods.begin(), methods.end(), "DMD");
    if (itDMD != methods.end()) {

        double tol = 1e-8;
        double t_0 = 0.0;
        Eigen::VectorXd t_vec( settings.Ns );
        t_vec(0) = t_0;

        for ( int i = 1; i < settings.Ns; i++ )
            t_vec(i) = t_vec(i-1) + settings.Dt_cfd*settings.Ds;

        std::cout << std::endl;
        std::cout << "Initialized vector of times " << std::endl;

        Eigen::VectorXd lambda_POD;
        Eigen::MatrixXd eig_vec_POD;
        Eigen::VectorXcd lambda_DMD;
        Eigen::MatrixXcd eig_vec_DMD;

        std::cout << "Extracting basis ... " << "\t";        
        Eigen::MatrixXcd Phi;
        Eigen::VectorXcd alfa;
        Eigen::MatrixXcd Alfas;

//        if ( settings.flag_method[0] == "DMD") {
            Phi = DMD_basis( sn_set,
                            lambda_DMD,
                            eig_vec_DMD,
                            lambda_POD,
                            eig_vec_POD,
                            settings.r );
//        }

//        if ( settings.flag_method[0] == "fbDMD") {
//            Phi = fbDMD_basis( sn_set,
//                            lambda_DMD,
//                            eig_vec_DMD,
//                            settings.r );
//        }
//
//        if ( settings.flag_method[0] == "HODMD") {
//            Phi = HODMD_basis( sn_set,
//                            lambda_DMD,
//                            eig_vec_DMD,
//                            alfa,
//                            tol,
//                            settings.d);
//        }

        int Nm = Phi.cols();
        std::cout << "Number of modes extracted : " << Nm << std::endl;

        std::cout << " Done! " << std::endl << std::endl;

        Eigen::VectorXcd omega(Nm);
        for ( int i = 0; i < Nm; i++ )
                omega(i) = std::log(lambda_DMD(i))/(settings.Dt_cfd*settings.Ds);

        // std::cout << " DMD omegas :\n " << omega << std::endl;
        // std::cout << " DMD eigen-values :\n " << lambda_DMD << std::endl;

        if ( settings.flag_method[0] == "DMD" || settings.flag_method[0] == "fbDMD" ) {
            std::cout << "Calculating coefficients DMD ... " << "\t";

        //Calculating coefficients solving optimization problem
        // Eigen::MatrixXcd alfa = Calculate_Coefs_DMD ( eig_vec_DMD,
                                            // eig_vec_POD,
                                            // lambda_DMD,
                                            // lambda_POD,
                                            // settings.Ns - 1 );

            if ( settings.dmd_coef_flag == "OPT" ) {
                alfa = Calculate_Coefs_DMD_exact ( sn_set.leftCols(settings.Ns-1),  //matrix of first Ns-1 snaps 
                                                                    lambda_DMD,  //slow eigenvalues
                                                                    Phi ); //slow exact DMD modes
            }

            //Calculating coefficients with ls
            else if ( settings.dmd_coef_flag == "LS" ) {
                Eigen::VectorXcd b = Eigen::VectorXcd::Zero(sn_set.rows()); 
                for ( int k = 0; k < sn_set.rows(); k++ )
                    b(k).real(sn_set(k,0)); 

                alfa = Phi.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);

            }
            //Calculating coefficients with Hybrid method
            else if ( settings.dmd_coef_flag == "HYBRID" ) {
            
                Alfas = Calculate_Coefs_Matrix_DMD ( sn_set,
                                                    Phi,
                                                    omega,
                                                    t_0,
                                                    settings.Dt_cfd*settings.Ds );


                std::cout << "Writing training points ..." << std::endl;
                std::ofstream train_real;
                train_real.open("train_real.dat");

                for ( int k = 0; k < settings.Ns; k++ ) {
                    for( int j = 0; j < Alfas.cols(); j++ ) 
                        train_real << Alfas(k,j).real() << " ";   

                train_real << std::endl;
                }
                train_real.close();

                std::ofstream train_imag;
                train_imag.open("train_imag.dat");

                for ( int k = 0; k < settings.Ns; k++ ) {
                    for( int j = 0; j < Alfas.cols(); j++ ) 
                        train_imag << Alfas(k,j).imag() << " ";
                train_imag << std::endl;
                }
                train_imag.close();

            } else {
                std::cout << "Method to Calculate DMD coefficients not available! " << std::endl;
                std::cout << "Exiting ... " << std::endl;
                std::exit( EXIT_FAILURE );
            }

            std::cout << " Done! " << std::endl << std::endl;
        }

        if ( settings.flag_wdb_be == "YES" && settings.dmd_coef_flag!= "HYBRID" ) {
            
            std::cout << "Writing modes ..." << "\t";
            write_modes_DMD ( Phi, Coords, settings.flag_prob );
            std::cout << "Complete!" << std::endl;



            std::cout << "Writing time dynamics ..." << "\t";
            write_TimeDynamics_DMD ( omega, alfa, t_vec );
            write_alfa_lam_DMD( alfa, lambda_DMD);
            std::cout << "Complete!" << std::endl;
            std::cout << std::endl;

        }

        if ( settings.flag_rec == "YES" ) {
            
            std::vector<double> t_st_vec(settings.Ns);
            t_st_vec[0] = t_0;

            for ( int i = 1; i < settings.Ns; i++ )
                t_st_vec[i] = t_st_vec[i-1] + settings.Dt_cfd*settings.Ds;

            for ( int nt = 0; nt < settings.t_rec.size(); nt ++) {
                Eigen::MatrixXcd Rec;
                std::cout << "Reconstructing field at time : " << settings.t_rec[nt] << "\t";

                if ( settings.dmd_coef_flag == "OPT" || settings.dmd_coef_flag == "LS" ) {
                    Rec = Reconstruction_DMD ( settings.t_rec[nt],
                                            settings.Dt_cfd*settings.Ds,
                                            alfa,
                                            Phi,
                                            lambda_DMD,
                                            settings.flag_prob );
                } else if ( settings.dmd_coef_flag == "HYBRID" ) {
                    Rec = Reconstruction_Hybrid_DMD ( settings.t_rec[nt],
                                                    t_st_vec,
                                                    Alfas,
                                                    Phi,
                                                    omega,
                                                    settings.flag_prob,
                                                    settings.flag_interp );
                } else {
                    std::cout << "Wrong method to calculate DMD coefficients! " << std::endl;
                    std::cout << "Exiting ... " << std::endl;
                    std::exit( EXIT_FAILURE );
                }

                std::cout << "Done" << std::endl;

                if ( settings.flag_mean == "YES" ) {
                    for ( int i = 0; i < Rec.cols(); i++)
                        Rec.col(i) = Rec.col(i) + mean.segment(i*Nr, Nr);

                }

                std::cout << "Writing reconstructed field ..." << "\t";

                write_Reconstructed_fields ( Rec.real(), Coords,
                                        settings.out_file,
                                        settings.flag_prob, nt );

                std::cout << "Done" << std::endl << std::endl;
            }
        }
    }

    //Defining common scope for mrDMD
    std::vector<std::string>::iterator itmrDMD;
    itmrDMD = std::find (methods.begin(), methods.end(), "mrDMD");
    if (itmrDMD != methods.end()) {

        Eigen::VectorXd t_vec( settings.Ns );
        t_vec(0) = 0.0;

        for ( int i = 1; i < settings.Ns; i++ )
            t_vec(i) = t_vec(i-1) + settings.Dt_cfd*settings.Ds;

        double dts = t_vec(1) - t_vec(0);

        std::cout << std::endl;
        std::cout << "Initialized vector of times " << std::endl;

        std::cout << "Computing nodes for Multi-Resolution Analysis... " << std::endl <<std::endl;        

        double t_0 = 0.0;
        int level = 0, bin_num = 0, offset = 0, max_levels = settings.max_levels, max_cycles = settings.max_cycles;
        std::vector<node_mrDMD> nodes = {}; 
        nodes = mrDMD_basis( sn_set,      
                            nodes,  
                            settings.r,                    
                            dts,                      
                            t_0,                
                            level,                          
                            bin_num,
                            offset,
                            max_levels,
                            max_cycles,
                            settings.dmd_coef_flag);

        std::cout << " Done! " << std::endl << std::endl;
        int sum = 0;
        for ( int n_nodes = 0; n_nodes < nodes.size(); n_nodes++ )
            sum += nodes[n_nodes].Modes.cols();

        std::cout << "Number of nodes stored : " << nodes.size() << std::endl;
        std::cout << "Number of total modes stored : " << sum << std::endl;

        if ( settings.flag_wdb_be == "YES" ) {
            std::cout << "Write coefs dynamics : " << std::endl;

            for ( int l = 0; l < max_levels; l ++ ) {
                std::cout << "Level " << l << "\t";
                write_CoefsDynamics_mrDMD( nodes, l, t_vec.size()/std::pow(2,l), max_levels);
                std::cout << "Done" << std::endl;
            }
        }

        if ( settings.flag_rec == "YES" ) {

            for ( int nt = 0; nt < settings.t_rec.size(); nt ++) {
                
                std::cout << "Reconstructing field at time : " << settings.t_rec[nt] << "\t";
                Eigen::MatrixXcd Rec = Reconstruction_mrDMD ( settings.t_rec[nt],                                                                                                                                                                                                                                                                                                                              
                                                            dts,       
                                                            nodes,     
                                                            settings.flag_prob );  

                std::cout << "Done" << std::endl;

                if ( settings.flag_mean == "YES" ) {
                    for ( int i = 0; i < Rec.cols(); i++)
                        Rec.col(i) = Rec.col(i) + mean.segment(i*Nr, Nr);

                }

                std::cout << "Writing reconstructed field ..." << "\t";

                write_Reconstructed_fields ( Rec.real(), Coords,
                                        settings.out_file,
                                        settings.flag_prob, nt );

                std::cout << "Done" << std::endl << std::endl;
            }
        }
    }

    //Defining common scope for RDMD
    std::vector<std::string>::iterator itRDMD;
    itRDMD = std::find (methods.begin(), methods.end(), "RDMD");
    if (itRDMD != methods.end())
    if ( settings.flag_method[0] == "RDMD") {

        double t_0 = 0.0;

        std::vector<double> t_st_vec(settings.Ns);
        t_st_vec[0] = t_0;

        for ( int i = 1; i < settings.Ns; i++ )
            t_st_vec[i] = t_st_vec[i-1] + settings.Dt_cfd*settings.Ds;

        std::cout << std::endl;
        std::cout << "Initialized vector of times " << std::endl;

        Eigen::VectorXd lambda = Eigen::VectorXd::Zero(settings.Ns);
        Eigen::MatrixXd Coefs = Eigen::MatrixXd::Zero(settings.Ns, settings.Ns);
        Eigen::MatrixXd Phi;


        if ( argc == 2 ) {
            std::cout << "Extracting basis and Coeffs RDMD ... " << "\t";        
        
            Phi = RDMD_modes_coefs ( sn_set,
                                    Coefs,
                                    lambda,
                                    K_pc,     
                                    settings.r,
                                    settings.r_RDMD,
                                    settings.En );
                                
        } else {

            std::cout << "Reading basis and extracting Coeffs RDMD ... " << "\t";
            std::string file_modes = argv[2];
            std::string file_coefs = argv[3];
            Phi = read_modes( file_modes, settings.ndim*Nr, settings.r_RDMD );
            Coefs = read_coefs( file_coefs, settings.Ns, settings.r_RDMD );

        }

        // std::cout << "Check mode orthogonality\n PhiT*Phi :\n " << Phi.transpose()*Phi << std::endl; 

        int Nm = Phi.cols();

        std::cout << " Done! " << std::endl << std::endl;

        if ( settings.flag_wdb_be == "YES" ) {
            std::cout << "Writing modes and coeffs..." << "\t"; //Writing one mode for all variables (ex (u,v)---> Phi1 = (Phiu1;Phiv1))
            write_modes ( Phi.leftCols(settings.r_RDMD) );
            write_coefs ( Coefs.leftCols(settings.r_RDMD));

            std::ofstream flow_data;
            flow_data.open("EnRDMD.dat"); 

            for ( int i = 0; i < settings.Ns; i++ )
                flow_data << std::setprecision(12) << std::scientific << K_pc(i) <<  " ";           
            // Close file
            flow_data.close();

            std::cout << "Complete!" << std::endl;

            // std::cout << "Writing Coefficients ..." << "\t";
            // write_coeffs_sPOD ( Coefs.transpose(), t_st_vec, lambda );
            // std::cout << "Complete!" << std::endl;
            // std::cout << std::endl;
        }


        if ( settings.flag_rec == "YES" ) {
                               
            for ( int nt = 0; nt < settings.t_rec.size(); nt ++) {
                Eigen::MatrixXd Rec;
                std::cout << "Reconstructing field at time : " << settings.t_rec[nt] << "\t";

                Rec = Reconstruction_RDMD ( settings.t_rec[nt],
                                        t_st_vec,
                                        Coefs,
                                        Phi,
                                        settings.flag_prob,
                                        settings.flag_interp );

                std::cout << "Done" << std::endl;

                if ( settings.flag_mean == "YES" ) {
                    for ( int i = 0; i < Rec.cols(); i++)
                        Rec.col(i) = Rec.col(i) + mean.segment(i*Nr, Nr);

                }

                std::cout << "Writing reconstructed field ..." << "\t";

                write_Reconstructed_fields ( Rec, Coords,
                                        settings.out_file,
                                        settings.flag_prob, nt );

                std::cout << "Done" << std::endl << std::endl;

            }
        }
    }

//    if ( settings.flag_method[0] == "GPOD") {
//
//        std::vector<double> t_vec( settings.Ns );
//        t_vec[0] = 0.0;
//        double Dt = settings.Dt_cfd*settings.Ds;
//
//        for ( int i = 1; i < settings.Ns; i++ )
//            t_vec[i] = t_vec[i-1] + settings.Dt_cfd*settings.Ds;
//
//        std::cout << std::endl;
//        std::cout << "Initialized vector of times " << std::endl;
//
//        Eigen::VectorXd lambda(settings.Ns);
//        Eigen::MatrixXd eig_vec(settings.Ns, settings.Ns);
//
//        std::cout << "Extracting basis ... " << "\t";
//        Eigen::MatrixXd Phi = GPOD_basis( Dt, sn_set,
//                                lambda, K_pc, eig_vec);
//
//        std::cout << " Done! " << std::endl << std::endl;
//
//        std::cout << "Number of non-zero modes : " << Phi.cols() << std::endl;
//
//        int Nrec;
//        if ( settings.r == 0 )
//        {
//            Nrec = Nmod( settings.En, K_pc);
//            std::cout << "Number of modes for the desired energy content (needs fixes for sPOD) : " << Nrec << std::endl;
//        }
//        else
//        {
//            Nrec = std::min(settings.r, settings.Ns -2);
//            std::cout << " Number of modes (fixed, needs fixes for sPOD) : " << Nrec << std::endl;
//        }
//
//
//        if ( settings.flag_wdb_be == "YES" )
//        {
//
//            std::cout << "Writing modes ..." << "\t";
//            write_modes_sPOD ( Phi.leftCols(Nrec), Coords, settings.flag_prob );
//            std::cout << "Complete!" << std::endl;
//
//            std::cout << "Writing Coefficients ..." << "\t";
//            write_coeffs_sPOD ( eig_vec.transpose(), t_vec, lambda );
//            std::cout << "Complete!" << std::endl;
//            std::cout << std::endl;
//
//        }
//
//        if ( settings.flag_rec == "YES" )
//        {
//
//            for ( int nt = 0; nt < settings.t_rec.size(); nt++ )
//            {
//
//                std::cout << "Reconstructing field at time : " << settings.t_rec[nt] << "\t";
//
//                Eigen::MatrixXd Rec = Reconstruction_RDMD ( settings.t_rec[nt],
//                                                        t_vec,
//                                                        eig_vec.transpose(),
//                                                        Phi,
//                                                        settings.flag_prob,
//                                                        settings.flag_interp );
//
//                std::cout << "Done" << std::endl;
//
//                if ( settings.flag_mean == "YES" )
//                {
//
//                    for ( int i = 0; i < Rec.cols(); i++)
//                        Rec.col(i) = Rec.col(i) + mean.segment(i*Nr, Nr);
//
//                }
//
//                std::cout << "Writing reconstructed field ..." << "\t";
//
//                write_Reconstructed_fields ( Rec, Coords,
//                                        settings.out_file,
//                                        settings.flag_prob, nt );
//
//                std::cout << "Done" << std::endl << std::endl;
//            }
//
//        }
//
//
//    }


    std::cout << std::endl;    
    std::cout << "-----------Single-MODES end-------------" << std::endl << std::endl;

    return 0;

}
