/*
CODE FOR PARETO FRONT ANALYSIS OF DIFFERENT RBM TECHNIQUES USING SU2 CODE
INPUT ARGUMENTS
Config File RBM + Config File SU2
*/

#include "Extract_Basis.hpp"
#include "read_Inputs.hpp"
#include "Pre-Process.hpp"
#include "System_Calls.hpp"
#include "Post-Process.hpp"
#include "Generate_snset.hpp"
#include "Reconstruction.hpp"
#include "write_Outputs.hpp"


int main( int argc, char *argv[] )
{
    
    std::cout << "-----------MODES Sensitivity analysis starts-------------" << std::endl << std::endl;

    std::cout << "Initializing common variables ... " << std::endl << std::endl;

    prob_settings settings;
    std::string filecfg = argv[1];
    std::string su2_conf = argv[2];
    Read_cfg( filecfg, settings );

    int nC = settings.Cols.size();
    std::vector<double> Dt_res = settings.Dt_res;

    int Nr;
    Eigen::MatrixXd Coords;
    std::vector<double> t_vec(settings.Ns);
    common_vars( Nr, Coords, t_vec, settings);

    std::cout << "Generating snapshot Matrix ... \n ";
    Eigen::MatrixXd sn_set = Eigen::MatrixXd::Zero(Nr,settings.Ns);
//    if ( settings.flag_wdb_be == "READ" && settings.flag_method[0] == "RDMD")
//        sn_set = generate_snap_matrix( Nr, settings);

    if ( settings.flag_wdb_be == "NO" || settings.flag_wdb_be == "WRITE"){
        std::cout << "Storing snapshot Matrix ... \n ";
        sn_set = generate_snap_matrix( Nr, settings);
    }


    std::cout << "Computing mean/Initial Condition of CFD solution ... " << std::endl;
    //Defining Initial condition
//    Eigen::VectorXd mean = sn_set.rowwise().mean();
    Eigen::VectorXd Ic = Eigen::VectorXd::Zero(nC*Nr);

    std::string binary = "YES";
    Eigen::VectorXd svd_cum_sum(settings.Ns);

    if ( settings.flag_mean == "IC" ) {
        std::cout << "Subtracting Initial condition" << std::endl << std::endl;
        Ic = IC(sn_set,settings,nC,Nr);
    } else {
        std::cout << "Using data without subtracting any reference state" << std::endl << std::endl;
    }


    auto methods = settings.flag_method;
    //Defining common scope for POD-SPOD
    std::vector<std::string>::iterator itPOD;
    itPOD = std::find (methods.begin(), methods.end(), "POD");
    if (itPOD != methods.end()) {
        std::cout << "--------------------------------------" << std::endl ;
        std::cout << "---Performin POD Nmodes sensitivity---" << std::endl ;
        std::cout << "--------------------------------------" << std::endl ;
        //Vector of MatrixXd where to store the evolution in time of conservative variables
        Eigen::MatrixXd Sn_Cons_time = Eigen::MatrixXd::Zero(nC*Nr, 3);
        std::vector<Eigen::MatrixXd> Phi(nC);
        std::vector<Eigen::VectorXd> lambda(nC);

        std::vector< std::vector<rbf> > surr_coefs(nC);
        Eigen::VectorXd K_pc(settings.Ns);
        Eigen::MatrixXd eig_vec(settings.Ns, settings.Ns);
        Eigen::VectorXi N_notZero(nC);
        //Check only for POD for now
        for (int i = 0; i < nC; i++) {
            Phi[i] = Eigen::MatrixXd::Zero(Nr,settings.Ns);
            lambda[i] = Eigen::VectorXd::Zero(settings.Ns);
        }
        std::cout << std::endl << "Extraction of the basis" << std::endl << std::endl;
        for ( int ncons = 0; ncons < nC; ncons ++ ) {
            std::cout << "Processing conservative variable " << ncons << std::endl;
            Phi[ncons] = SPOD_basis(sn_set.middleRows(ncons * Nr, Nr),
                                    lambda[ncons], K_pc, eig_vec,
                                    0,
                                    settings.flag_bc,
                                    settings.flag_filter,
                                    settings.sigma);
            N_notZero(ncons) = Phi[ncons].cols();
            surr_coefs[ncons] = getSurrCoefs(t_vec, eig_vec, settings.flag_interp);
        }
        std::cout << std::endl;

        for ( int idtr = 0; idtr < settings.Dt_res.size(); idtr++ ) {
            std::cout << " --------------DT_RES = " << settings.Dt_res[idtr] << "--------------"<< std::endl;
            for ( int in_mode = settings.init_imode; in_mode <= N_notZero.maxCoeff(); in_mode++ ) {
                std::cout << "For iMode " << in_mode << std::endl;
                for (int itr = settings.init_tres; itr < settings.t_res.size(); itr++) {
                    std::cout << "Computing residuals at time t = " << settings.t_res[itr] << std::endl;

//                    std::vector<Eigen::MatrixXd> Coef_Rec;
                    for (int ncons = 0; ncons < nC; ncons++) {

                        int Nm = in_mode;
                        if ( Nm > N_notZero[ncons]) Nm = N_notZero[ncons];
                        Eigen::MatrixXd coef_t(3, Nm);

                        std::vector<double> tr(1);
                        std::vector<double> t_evaluate = {settings.t_res[itr] - 2.0 * settings.Dt_res[idtr],
                                                          settings.t_res[itr] - settings.Dt_res[idtr],
                                                          settings.t_res[itr]};

                        for (int j = 0; j < 3; j++) {
                            tr[0] = t_evaluate[j];
                            for (int i = 0; i < Nm; i++)
                                surr_coefs[ncons][i].evaluate(tr, coef_t(j, i));
                        }
//                        Coef_Rec.push_back(coef_t);

                        Eigen::MatrixXd Sig = Eigen::MatrixXd::Zero(Nm, Nm);
                        for (int i = 0; i < Nm; i++)
                            Sig(i, i) = std::sqrt(lambda[ncons](i));

                        Sn_Cons_time.middleRows(ncons * Nr, Nr) = Phi[ncons].leftCols(Nm) * (Sig * coef_t.transpose());

                        }
                        if (settings.flag_mean == "IC") {
                            for (int it = 0; it < 3; it++)
                                Sn_Cons_time.col(it) += Ic;
                        }
                        settings.r = in_mode;
                        //Launching SU2_DTR and saving errors and Residuals to file
                        int iter = std::round(settings.t_res[itr] / settings.Dt_cfd);
                        Write_Restart_Cons_Time(Sn_Cons_time, Coords, settings.out_file, iter, nC, settings.alpha,
                                                settings.beta, binary);
                        SU2_DTR(settings, su2_conf, "POD", idtr, itr);

                }
            }
        }
    }

    //Defining scope for SPOD
    std::vector<std::string>::iterator itSPOD;
    itSPOD = std::find (methods.begin(), methods.end(), "SPOD");
    if (itSPOD != methods.end()) {
        std::cout << "--------------------------------------" << std::endl;
        std::cout << "---Performin SPOD Nmodes sensitivity--" << std::endl;
        std::cout << "--------------------------------------" << std::endl;
        //Vector of MatrixXd where to store the evolution in time of conservative variables
        Eigen::MatrixXd Sn_Cons_time = Eigen::MatrixXd::Zero(nC * Nr, 3);
        std::vector<Eigen::MatrixXd> Phi(nC);
        std::vector<Eigen::VectorXd> lambda(nC);

        std::vector<std::vector<rbf> > surr_coefs(nC);
        Eigen::VectorXd K_pc(settings.Ns);
        Eigen::MatrixXd eig_vec(settings.Ns, settings.Ns);
        Eigen::VectorXi N_notZero(nC);
        //Check only for POD for now
        for (int i = 0; i < nC; i++) {
            Phi[i] = Eigen::MatrixXd::Zero(Nr, settings.Ns);
            lambda[i] = Eigen::VectorXd::Zero(settings.Ns);
        }

        for (int iNf = 0; iNf < settings.Nf.size(); iNf++) {

            std::cout << std::endl << "Extraction of the basis for Nf " << settings.Nf[iNf] << std::endl << std::endl;
            std::string method = "SPOD" + std::to_string(settings.Nf[iNf]);
            for (int ncons = 0; ncons < nC; ncons++) {
                std::cout << "Processing conservative variable " << ncons << std::endl;
                Phi[ncons] = SPOD_basis(sn_set.middleRows(ncons * Nr, Nr),
                                        lambda[ncons], K_pc, eig_vec,
                                        settings.Nf[iNf],
                                        settings.flag_bc,
                                        settings.flag_filter,
                                        settings.sigma);
                N_notZero[ncons] = Phi[ncons].cols();
                surr_coefs[ncons] = getSurrCoefs(t_vec, eig_vec, settings.flag_interp);
            }
            std::cout << std::endl;

            for (int idtr = 0; idtr < settings.Dt_res.size(); idtr++) {
                std::cout << " --------------DT_RES = " << settings.Dt_res[idtr] << "--------------" << std::endl;

                for (int itr = settings.init_tres; itr < settings.t_res.size(); itr++) {
                    std::cout << "Computing residuals at time t = " << settings.t_res[itr] << std::endl;
                    std::vector<Eigen::MatrixXd> Coef_Rec;
                    for (int ncons = 0; ncons < nC; ncons++) {

                        Eigen::MatrixXd coef_t(3, N_notZero(ncons));

                        std::vector<double> tr(1);
                        std::vector<double> t_evaluate = {settings.t_res[itr] - 2.0 * settings.Dt_res[idtr],
                                                          settings.t_res[itr] - settings.Dt_res[idtr],
                                                          settings.t_res[itr]};

                        for (int j = 0; j < 3; j++) {
                            tr[0] = t_evaluate[j];
                            for (int i = 0; i < N_notZero(ncons); i++)
                                surr_coefs[ncons][i].evaluate(tr, coef_t(j, i));
                        }
                        Coef_Rec.push_back(coef_t);
                    }

                    for (int in_mode = settings.init_imode; in_mode <= N_notZero.maxCoeff(); in_mode++) {
                        std::cout << "For iMode " << in_mode << std::endl;
                        //    }
                        for (int ncons = 0; ncons < nC; ncons++) {

                            int Nm = in_mode;
                            if (Nm > N_notZero[ncons]) Nm = N_notZero[ncons];

                            Eigen::MatrixXd Sig = Eigen::MatrixXd::Zero(Nm, Nm);
                            for (int i = 0; i < Nm; i++)
                                Sig(i, i) = std::sqrt(lambda[ncons](i));
                            Sn_Cons_time.middleRows(ncons * Nr, Nr) =
                                    Phi[ncons].leftCols(Nm) * Sig * Coef_Rec[ncons].leftCols(Nm).transpose();

                        }
                        if (settings.flag_mean == "IC") {
                            for (int it = 0; it < 3; it++)
                                Sn_Cons_time.col(it) += Ic;
                        }
                        settings.r = in_mode;
                        //Launching SU2_DTR and saving errors and Residuals to file
                        int iter = std::round(settings.t_res[itr] / settings.Dt_cfd);
                        Write_Restart_Cons_Time(Sn_Cons_time, Coords, settings.out_file, iter, nC, settings.alpha,
                                                settings.beta, binary);
                        SU2_DTR(settings, su2_conf, method, idtr, itr);

                    }
                }

            }
        }
    }


//Defining scope for DMD ( Rank=-1 preferable, Coeffs = OPT )
    std::vector<std::string>::iterator itDMD;
    itDMD = std::find (methods.begin(), methods.end(), "DMD");
    if (itDMD != methods.end()) {
        std::cout << "--------------------------------------" << std::endl;
        std::cout << "---Performin DMD Nmodes sensitivity---" << std::endl;
        std::cout << "--------------------------------------" << std::endl;
        Eigen::MatrixXd Sn_Cons_time = Eigen::MatrixXd::Zero(nC * Nr, 3);
        Eigen::VectorXd lambda_POD;
        Eigen::MatrixXd eig_vec_POD;
        std::vector<Eigen::VectorXcd> lambda_DMD(nC);
        Eigen::MatrixXcd eig_vec_DMD;
        std::vector<Eigen::MatrixXcd> Phi(nC);
        std::vector<Eigen::VectorXcd> alfa(nC);
        for ( int in_mode = settings.init_imode; in_mode < settings.Ns; in_mode++ ) {

            settings.r = in_mode;
            for (int i = 0; i < nC; i++) {
                    Phi[i] = Eigen::MatrixXcd::Zero(Nr, settings.Ns);
                    alfa[i] = Eigen::VectorXcd::Zero(settings.Ns);
                    lambda_DMD[i] = Eigen::VectorXcd::Zero(settings.Ns);
                }

            for (int ncons = 0; ncons < nC; ncons++) {

                std::cout << "Processing conservative variable " << ncons << std::endl;

                Phi[ncons] = DMD_basis(sn_set.middleRows(ncons * Nr, Nr),
                                       lambda_DMD[ncons],
                                       eig_vec_DMD,
                                       lambda_POD,
                                       eig_vec_POD,
                                       in_mode);

                std::cout << "Number of DMD modes extracted : " << Phi[ncons].cols() << std::endl;

                Eigen::VectorXcd omega(Phi[ncons].cols());
                for (int i = 0; i < Phi[ncons].cols(); i++)
                    omega(i) = std::log(lambda_DMD[ncons](i)) / (settings.Dt_cfd * settings.Ds);

                //Computing alpha if exponential rec is needed
                if (settings.dmd_coef_flag == "OPT") {

                    alfa[ncons] = Calculate_Coefs_DMD_exact(sn_set.middleRows(ncons * Nr, Nr).leftCols(settings.Ns - 1),
                                                            lambda_DMD[ncons],
                                                            Phi[ncons]);
                } else if (settings.dmd_coef_flag == "LS") {

                    Eigen::VectorXcd b = Eigen::VectorXcd::Zero(sn_set.rows());
                    for (int k = 0; k < sn_set.rows(); k++) b(k).real(sn_set(k, 0));
                    alfa[ncons] = Phi[ncons].jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);

                }

            }

            for (int idtr = 0; idtr < settings.Dt_res.size(); idtr++) {

                std::cout << " --------------DT_RES = " << settings.Dt_res[idtr] << "--------------" << std::endl;
                for (int itr = settings.init_tres; itr < settings.t_res.size(); itr++) {
                    std::cout << "Computing residual at time t = " << settings.t_res[itr] << std::endl;

                    std::vector<double> t_evaluate = {settings.t_res[itr] - 2.0 * settings.Dt_res[idtr],
                                                      settings.t_res[itr] - settings.Dt_res[idtr],
                                                      settings.t_res[itr]};

                    for (int ncons = 0; ncons < nC; ncons++) {

                        for (int j = 0; j < 3; j++) {
                            double t_dmd = t_evaluate[j] - t_vec[0];
                            Eigen::MatrixXcd Rec = Reconstruction_DMD(t_dmd,
                                                                      settings.Dt_cfd * settings.Ds,
                                                                      alfa[ncons],
                                                                      Phi[ncons],
                                                                      lambda_DMD[ncons],
                                                                      "SCALAR");

                            Sn_Cons_time.middleRows(ncons * Nr, Nr).col(j) = Rec.real();

                        }

                    }

                    if (settings.flag_mean == "IC") {
                        for (int it = 0; it < 3; it++)
                            Sn_Cons_time.col(it) += Ic;
                    }

                    //Launching SU2_DTR and saving errors and Residuals to file
                    int iter = std::round(settings.t_res[itr] / settings.Dt_cfd);
                    Write_Restart_Cons_Time(Sn_Cons_time, Coords, settings.out_file, iter, nC, settings.alpha,
                                            settings.beta, binary);
                    SU2_DTR(settings, su2_conf, "DMD", idtr, itr);
    //                Write_History_ResError(settings, "DMD", idtr, itr);
                    std::cout << std::endl;
                }
            }

        }
    }

// // //Defining scope for RDMD
// // //if using the function RDMD_modes_coefs for energybased select
// // //energy level and rank rdmd to zero, for mode based just select
// //rank rdmd to the number of desired modes
    std::vector<std::string>::iterator itRDMD;
    itRDMD = std::find (methods.begin(), methods.end(), "RDMD");
    if (itRDMD != methods.end()) {

        std::cout << "--------------------------------------" << std::endl ;
        std::cout << "--Performin RDMD Nmodes sensitivity---" << std::endl ;
        std::cout << "--------------------------------------" << std::endl ;
        //Vector of MatrixXd where to store the evolution in time of conservative variables
        Eigen::MatrixXd Sn_Cons_time = Eigen::MatrixXd::Zero(nC*Nr, 3);
        std::vector<Eigen::MatrixXd> Phi(nC);
        std::vector<Eigen::MatrixXd> COEF(nC);
        std::vector< std::vector<rbf> > surr_coefs(nC);
        std::vector<Eigen::VectorXd> lambda(nC);

        Eigen::VectorXd K_pc(settings.Ns);
        Eigen::MatrixXd Coefs = Eigen::MatrixXd::Zero(settings.Ns, settings.Ns);
        int Nm;
        Eigen::VectorXi N_notZero(nC);


        for (int i = 0; i < nC; i++) {
            Phi[i] = Eigen::MatrixXd::Zero(Nr,settings.Ns);
            COEF[i] = Eigen::MatrixXd::Zero(settings.Ns,settings.Ns);
            lambda[i] = Eigen::VectorXd::Zero(settings.Ns);
        }

        if ( settings.flag_wdb_be == "WRITE" || settings.flag_wdb_be == "NO") {
            for (int ncons = 0; ncons < nC; ncons++) {

                std::cout << "Processing conservative variable " << ncons << std::endl;
                Phi[ncons] = RDMD_modes_coefs(sn_set.middleRows(ncons * Nr, Nr),
                                              Coefs,
                                              lambda[ncons],
                                              K_pc,
                                              -1, //Performing singular value hard threshold for DMD reduction at each step
                                              settings.Ns,
                                              settings.En);

                COEF[ncons] = Coefs.transpose();
                surr_coefs[ncons] = getSurrCoefs(t_vec,
                                                 Coefs.transpose(),
                                                 settings.flag_interp);

                N_notZero(ncons) = Phi[ncons].cols();

            }

            if (settings.flag_wdb_be == "WRITE") ModesDB_Write(Phi, COEF, settings);

        } else if ( settings.flag_wdb_be == "READ") {
            std::cout << "Reading RDMD basis" << std::endl;
            settings.r = settings.Ns;
            ModeDB_Read("ModesRDMD_", "CoefsRDMD_", Phi, COEF, settings);

            for ( int icons = 0; icons < nC; icons++ ) {
                surr_coefs[icons] = getSurrCoefs(t_vec,
                                                 COEF[icons],
                                                 settings.flag_interp);
                N_notZero(icons) = Phi[icons].cols();
            }

        }

        for ( int idtr = 0; idtr < settings.Dt_res.size(); idtr++ ) {
            std::cout << " --------------DT_RES = " << settings.Dt_res[idtr] << "--------------"<< std::endl;
            for ( int in_mode = settings.init_imode; in_mode <= N_notZero.maxCoeff(); in_mode++ ) {
                std::cout << "For iMode " << in_mode << std::endl;
                for (int itr = settings.init_tres; itr < settings.t_res.size(); itr++) {
//                    std::vector<Eigen::MatrixXd> Coefs_Rec;
                    std::cout << " Computing residuals at time = " << settings.t_res[itr] << std::endl;

                    for (int ncons = 0; ncons < nC; ncons++) {
                        int Nm = in_mode;
                        if (Nm > N_notZero[ncons]) Nm = N_notZero[ncons];
                        Eigen::MatrixXd coef_t(3, Nm);
                        std::vector<double> tr(1);
                        std::vector<double> t_evaluate = {settings.t_res[itr] - 2.0 * settings.Dt_res[idtr],
                                                          settings.t_res[itr] - settings.Dt_res[idtr],
                                                          settings.t_res[itr]};

                        for (int j = 0; j < 3; j++) {
                            tr[0] = t_evaluate[j];
                            for (int i = 0; i < Nm; i++)
                                surr_coefs[ncons][i].evaluate(tr, coef_t(j, i));
                        }
//                        Coefs_Rec.push_back(coef_t);


    //                        Eigen::MatrixXd Sig = Eigen::MatrixXd::Zero(Nm, Nm);
    //                        for (int i = 0; i < Nm; i++)
    //                            Sig(i, i) = std::sqrt(lambda[ncons](i));
                        Sn_Cons_time.middleRows(ncons * Nr, Nr) =
                                Phi[ncons].leftCols(Nm) * coef_t.transpose();

                    }
                    if (settings.flag_mean == "IC") {
                        for (int it = 0; it < 3; it++)
                            Sn_Cons_time.col(it) += Ic;
                    }
                    settings.r = in_mode;
                    //Launching SU2_DTR and saving errors and Residuals to file
                    int iter = std::round(settings.t_res[itr] / settings.Dt_cfd);
                    Write_Restart_Cons_Time(Sn_Cons_time, Coords, settings.out_file, iter, nC, settings.alpha,
                                            settings.beta, binary);
                    SU2_DTR(settings, su2_conf, "RDMD", idtr, itr);
                }

            }

        }

    }


    std::cout << "MODES sensitivity analysis ends" << std::endl << std::endl;
    return 0;

}
