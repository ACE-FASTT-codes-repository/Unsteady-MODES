/*
CODE FOR RESIDUAL EVALUATION OF DIFFERENT RBM TECHNIQUES USING SU2 CODE
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

    std::cout << "-----------MODES Adaptive Residual Evaluation starts-------------" << std::endl << std::endl;

    std::cout << "Initializing common variables ... " << std::endl << std::endl;

    prob_settings settings;
    std::string filecfg = argv[1];
    std::string su2_conf = argv[2];
    Read_cfg( filecfg, settings );

    // Pre-processing common vars
    int nC = settings.Cols.size();
    std::vector<double> Dt_res = settings.Dt_res;
    int Nr;
    Eigen::MatrixXd Coords;
    std::vector<double> t_vec(settings.Ns);
    common_vars( Nr, Coords, t_vec, settings);

    std::cout << "Storing snapshot Matrix ... \n ";
    Eigen::MatrixXd sn_set = generate_snap_matrix( Nr, settings);

    std::cout << std::endl;

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


    auto methods = settings.flag_method;
    //Defining common scope for POD-SPOD
    std::vector<std::string>::iterator itPOD;
    itPOD = std::find (methods.begin(), methods.end(), "POD");
    if (itPOD != methods.end()) {
        std::cout << "--------------------------------------" << std::endl ;
        std::cout << "-------Performin POD ResEval----------" << std::endl ;
        std::cout << "--------------------------------------" << std::endl ;
        //Vector of MatrixXd where to store the evolution in time of conservative variables
        Eigen::MatrixXd Sn_Cons_time = Eigen::MatrixXd::Zero(nC*Nr, 3);
        std::vector<Eigen::MatrixXd> Phi(nC);
        std::vector<Eigen::VectorXd> lambda(nC);

        std::vector< std::vector<rbf> > surr_coefs(nC);
        Eigen::VectorXd K_pc(settings.Ns);
        Eigen::MatrixXd eig_vec(settings.Ns, settings.Ns);
        std::vector<int> Nm(nC);
        int N_notZero;
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
            N_notZero = Phi[ncons].cols();
            if (settings.r == 0) Nm[ncons] = Nmod(settings.En, K_pc);
            else Nm[ncons] = std::min(settings.r, N_notZero);
            std::cout << "Number of modes used in reconstruction " << Nm[ncons] << std::endl;
            surr_coefs[ncons] = getSurrCoefs(t_vec, eig_vec, settings.flag_interp);
        }
        std::cout << std::endl;

//        std::cout << "Computing SPOD " << Nf[nfj] << " reconstruction for each conservative variable ... " << "\n";

        for ( int idtr = 0; idtr < settings.Dt_res.size(); idtr++ ) {
            std::cout << " --------------DT_RES = " << settings.Dt_res[idtr] << "--------------"<< std::endl;

            for (int itr = 0; itr < settings.t_res.size(); itr++) {
                std::cout << "Computing residuals at time t = " << settings.t_res[itr] << std::endl;

                for (int ncons = 0; ncons < nC; ncons++) {

                    Eigen::MatrixXd coef_t(3, Nm[ncons]);

                    std::vector<double> tr(1);
                    std::vector<double> t_evaluate = {settings.t_res[itr] - 2.0 * settings.Dt_res[idtr],
                                                      settings.t_res[itr] - settings.Dt_res[idtr],
                                                      settings.t_res[itr]};

                    for (int j = 0; j < 3; j++) {
                        tr[0] = t_evaluate[j];
                        for (int i = 0; i < Nm[ncons]; i++)
                            surr_coefs[ncons][i].evaluate(tr, coef_t(j, i));
                    }

                //    }
                    Eigen::MatrixXd Sig = Eigen::MatrixXd::Zero(Nm[ncons], Nm[ncons]);
                    for (int i = 0; i < Nm[ncons]; i++)
                        Sig(i, i) = std::sqrt(lambda[ncons](i));
                    Sn_Cons_time.middleRows(ncons * Nr, Nr) = Phi[ncons].leftCols(Nm[ncons]) * Sig * coef_t.transpose();
                }

                if (settings.flag_mean == "IC") {
                    for (int it = 0; it < 3; it++)
                        Sn_Cons_time.col(it) += Ic;
                }

                //Launching SU2_DTR and saving errors and Residuals to file
                int iter = std::round(settings.t_res[itr]/settings.Dt_cfd);
                Write_Restart_Cons_Time(Sn_Cons_time, Coords, settings.out_file, iter, nC, settings.alpha, settings.beta, binary);
                SU2_DTR(settings, su2_conf, "POD", idtr, itr);
//                Write_History_ResError(settings, "POD", idtr, itr);
                std::cout << std::endl;
            }
        }
    }

    //Defining scope for SPOD
    std::vector<std::string>::iterator itSPOD;
    itSPOD = std::find (methods.begin(), methods.end(), "SPOD");
    if (itSPOD != methods.end()) {
        std::cout << "--------------------------------------" << std::endl ;
        std::cout << "-------Performin SPOD ResEval----------" << std::endl ;
        std::cout << "--------------------------------------" << std::endl ;
        //Vector of MatrixXd where to store the evolution in time of conservative variables
        Eigen::MatrixXd Sn_Cons_time = Eigen::MatrixXd::Zero(nC*Nr, 3);
        std::vector<Eigen::MatrixXd> Phi(nC);
        std::vector<Eigen::VectorXd> lambda(nC);

        std::vector< std::vector<rbf> > surr_coefs(nC);
        Eigen::VectorXd K_pc(settings.Ns);
        Eigen::MatrixXd eig_vec(settings.Ns, settings.Ns);
        std::vector<int> Nm(nC);
        int N_notZero;
        //Check only for POD for now
        for (int i = 0; i < nC; i++) {
            Phi[i] = Eigen::MatrixXd::Zero(Nr,settings.Ns);
            lambda[i] = Eigen::VectorXd::Zero(settings.Ns);
        }

        for ( int iNf = 0; iNf < settings.Nf.size(); iNf++ ) {

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
                N_notZero = Phi[ncons].cols();
                if (settings.r == 0) Nm[ncons] = Nmod(settings.En, K_pc);
                else Nm[ncons] = std::min(settings.r, N_notZero);
                std::cout << "Number of modes used in reconstruction " << Nm[ncons] << std::endl;
                surr_coefs[ncons] = getSurrCoefs(t_vec, eig_vec, settings.flag_interp);
            }
            std::cout << std::endl;

//        std::cout << "Computing SPOD " << Nf[nfj] << " reconstruction for each conservative variable ... " << "\n";

            for (int idtr = 0; idtr < settings.Dt_res.size(); idtr++) {
                std::cout << " --------------DT_RES = " << settings.Dt_res[idtr] << "--------------" << std::endl;

                for (int itr = 0; itr < settings.t_res.size(); itr++) {
                    std::cout << "Computing residuals at time t = " << settings.t_res[itr] << std::endl;

                    for (int ncons = 0; ncons < nC; ncons++) {

                        Eigen::MatrixXd coef_t(3, Nm[ncons]);

                        std::vector<double> tr(1);
                        std::vector<double> t_evaluate = {settings.t_res[itr] - 2.0 * settings.Dt_res[idtr],
                                                          settings.t_res[itr] - settings.Dt_res[idtr],
                                                          settings.t_res[itr]};

                        for (int j = 0; j < 3; j++) {
                            tr[0] = t_evaluate[j];
                            for (int i = 0; i < Nm[ncons]; i++)
                                surr_coefs[ncons][i].evaluate(tr, coef_t(j, i));
                        }

                        //    }
                        Eigen::MatrixXd Sig = Eigen::MatrixXd::Zero(Nm[ncons], Nm[ncons]);
                        for (int i = 0; i < Nm[ncons]; i++)
                            Sig(i, i) = std::sqrt(lambda[ncons](i));
                        Sn_Cons_time.middleRows(ncons * Nr, Nr) =
                                Phi[ncons].leftCols(Nm[ncons]) * Sig * coef_t.transpose();
                    }

                    if (settings.flag_mean == "IC") {
                        for (int it = 0; it < 3; it++)
                            Sn_Cons_time.col(it) += Ic;
                    }

                    //Launching SU2_DTR and saving errors and Residuals to file
                    int iter = std::round(settings.t_res[itr] / settings.Dt_cfd);
                    Write_Restart_Cons_Time(Sn_Cons_time, Coords, settings.out_file, iter, nC, settings.alpha,
                                            settings.beta, binary);
                    SU2_DTR(settings, su2_conf, method, idtr, itr);
//                Write_History_ResError(settings, "POD", idtr, itr);
                    std::cout << std::endl;
                }
            }

        }
    }


//Defining scope for DMD ( Rank=-1 preferable, Coeffs = OPT )
    std::vector<std::string>::iterator itDMD;
    itDMD = std::find (methods.begin(), methods.end(), "DMD");
    if (itDMD != methods.end()) {
        std::cout << "--------------------------------------" << std::endl ;
        std::cout << "-------Performin DMD ResEval----------" << std::endl ;
        std::cout << "--------------------------------------" << std::endl ;
        Eigen::MatrixXd Sn_Cons_time = Eigen::MatrixXd::Zero(nC*Nr, 3);
        Eigen::VectorXd lambda_POD;
        Eigen::MatrixXd eig_vec_POD;
        std::vector<Eigen::VectorXcd> lambda_DMD(nC);
        Eigen::MatrixXcd eig_vec_DMD;
        std::vector<Eigen::MatrixXcd> Phi(nC);
        std::vector<Eigen::VectorXcd> alfa(nC);
        std::vector<int> Nm(nC);

        //Defining quantities needed for uniform interpolation
        std::vector<std::vector<smartuq::surrogate::rbf> > surr_coefs_DMD_r(nC);
        std::vector<std::vector<smartuq::surrogate::rbf> > surr_coefs_DMD_i(nC);


        for ( int i = 0; i < nC; i++ ) {
            Phi[i] = Eigen::MatrixXcd::Zero(Nr, settings.Ns);
            alfa[i] = Eigen::VectorXcd::Zero(settings.Ns);
            lambda_DMD[i] = Eigen::VectorXcd::Zero(settings.Ns);
        }
        std::cout << std::endl << "Extraction of the basis" << std::endl << std::endl;
        for ( int ncons = 0; ncons < nC; ncons++ ) {

            std::cout << "Processing conservative variable " << ncons << std::endl;
//            std::cout << "Extracting basis DMD using rank " << std::min(settings.r, settings.Ns - 1) << "\t";

            if (settings.r == 0) {
                Phi[ncons] = DMD_basis(sn_set.middleRows(ncons * Nr, Nr),
                                       lambda_DMD[ncons],
                                       eig_vec_DMD,
                                       lambda_POD,
                                       eig_vec_POD,
                                       -1);
            } else {
                Phi[ncons] = DMD_basis(sn_set.middleRows(ncons * Nr, Nr),
                                       lambda_DMD[ncons],
                                       eig_vec_DMD,
                                       lambda_POD,
                                       eig_vec_POD,
                                       settings.r);

                std::cout << "Number of DMD modes extracted : " << Phi[ncons].cols() << std::endl;
            }

            Eigen::VectorXcd omega(Phi[ncons].cols());
            for (int i = 0; i < Phi[ncons].cols(); i++)
                omega(i) = std::log(lambda_DMD[ncons](i)) / (settings.Dt_cfd * (double)settings.Ds);

            //Computing alpha if exponential rec is needed
            if ( settings.dmd_coef_flag == "OPT") {

                alfa[ncons] = Calculate_Coefs_DMD_exact(sn_set.middleRows(ncons * Nr, Nr).leftCols(settings.Ns - 1),
                                                        lambda_DMD[ncons],
                                                        Phi[ncons]);
            } else if ( settings.dmd_coef_flag == "LS") {

                Eigen::VectorXcd b = Eigen::VectorXcd::Zero(sn_set.rows());
                for ( int k = 0; k < sn_set.rows(); k++ ) b(k).real(sn_set(k,0));
                alfa[ncons] = Phi[ncons].jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);

            }

            Eigen::VectorXd En = Eigen::VectorXd::Zero(Phi[ncons].cols());
            double T = t_vec[t_vec.size() - 1];

            for (int i = 0; i < Phi[ncons].cols(); i++) {

                double alfa_i = alfa[ncons](i).imag();
                double alfa_r = alfa[ncons](i).real();
                double sigma = omega(i).real();
                En(i) = (alfa_r * alfa_r + alfa_i * alfa_i) * (std::exp(2.0 * sigma * T) - 1.0) / (2.0 * sigma);

            }

            dmd_sort(En, Phi[ncons], lambda_DMD[ncons], alfa[ncons]);
//            dmd_tenvelope_sort( Phi[ncons], omega, alfa[ncons], t_vec );

            double sum = 0;
            Eigen::VectorXd K_pc = Eigen::VectorXd::Zero(settings.Ns);
            for (int i = 0; i < Phi[ncons].cols(); i++) {
                sum += En(i) / En.sum();
                K_pc(i) = sum;
            }
            if (settings.r == 0) {
                Nm[ncons] = Nmod(settings.En, K_pc);
                std::cout << "Number of modes for the desired energetic content : " << Nm[ncons] << std::endl;
            } else {
                Nm[ncons] = std::min(settings.r, settings.Ns - 1);
                std::cout << "Number of modes (fixed) : " << Nm[ncons] << std::endl;
            }
            //Defining quantities to get a uniform DMD reconstruction if needed (not using exprec)
//            Eigen::MatrixXcd PhiTPhi = Phi[ncons].leftCols(Nm[ncons]).transpose()*Phi[ncons].leftCols(Nm[ncons]);
//            Eigen::MatrixXcd Coeffs = PhiTPhi.inverse()*(Phi[ncons].leftCols(Nm[ncons]).transpose()*sn_set.middleRows(ncons*Nr,Nr));
            Eigen::MatrixXcd PhiTPhi = Phi[ncons].transpose()*Phi[ncons];
            Eigen::MatrixXcd Coeffs = PhiTPhi.inverse()*(Phi[ncons].transpose()*sn_set.middleRows(ncons*Nr,Nr));
//            surr_coefs_DMD_r[ncons].resize(Nm[ncons]);
//            surr_coefs_DMD_i[ncons].resize(Nm[ncons]);
            surr_coefs_DMD_r[ncons] = getSurrCoefs(t_vec, Coeffs.real().transpose(), settings.flag_interp);
            surr_coefs_DMD_i[ncons] = getSurrCoefs(t_vec, Coeffs.imag().transpose(), settings.flag_interp);

        }

        std::cout << std::endl;



        for ( int idtr = 0; idtr < settings.Dt_res.size(); idtr++) {

            std::cout << " --------------DT_RES = " << settings.Dt_res[idtr] << "--------------"<< std::endl;
            for (int itr = 0; itr < settings.t_res.size(); itr++) {
                std::cout << "Computing residual at time t = " << settings.t_res[itr] << std::endl;

                std::vector<double> t_evaluate = {settings.t_res[itr] - 2.0 * settings.Dt_res[idtr],
                                                  settings.t_res[itr] - settings.Dt_res[idtr],
                                                  settings.t_res[itr]};

                for (int ncons = 0; ncons < nC; ncons++) {


                    if ( settings.flag_rec == "UNIFORM" ) {
                        std::cout << "Using POD reconstruction formula for DMD" << std::endl;
                        double tmp_r, tmp_i;
                        std::vector<double> tr(1);
                        Eigen::MatrixXcd coef_t(3, Nm[ncons]);
                        //DMD using same rec formula as POD
//                        --------------------------------------------------------------------------------------------
                        for (int j = 0; j < 3; j++) {
                            tr[0] = t_evaluate[j];
                            for (int i = 0; i < Nm[ncons]; i++) {
                                surr_coefs_DMD_r[ncons][i].evaluate(tr, tmp_r);
                                surr_coefs_DMD_i[ncons][i].evaluate(tr, tmp_i);
                                std::complex<double> c(tmp_r,tmp_i);
                                coef_t(j,i) = c;
                            }
                        }

                        Eigen::MatrixXcd Appo = Phi[ncons].leftCols(Nm[ncons]) * coef_t.transpose();
                        Sn_Cons_time.middleRows(ncons * Nr, Nr) = Appo.real();
//                      ----------------------------------------------------------------------------------------------
                    } else {
                        std::cout << "Using exponential DMD reconstruction" << std::endl;

                        for (int j = 0; j < 3; j++) {
                            double t_dmd = t_evaluate[j] - t_vec[0];
                            Eigen::MatrixXcd Rec = Reconstruction_DMD(t_dmd,
                                                                      settings.Dt_cfd * settings.Ds,
                                                                      alfa[ncons].head(Nm[ncons]),
                                                                      Phi[ncons].leftCols(Nm[ncons]),
                                                                      lambda_DMD[ncons].head(Nm[ncons]),
                                                                      "SCALAR");

                            Sn_Cons_time.middleRows(ncons * Nr, Nr).col(j) = Rec.real();

                        }
                    }
                }

                if (settings.flag_mean == "IC") {
                    for (int it = 0; it < 3; it++)
                        Sn_Cons_time.col(it) += Ic;
                }

                //Launching SU2_DTR and saving errors and Residuals to file
                int iter = std::round(settings.t_res[itr]/settings.Dt_cfd);
                Write_Restart_Cons_Time(Sn_Cons_time, Coords, settings.out_file, iter, nC, settings.alpha, settings.beta, binary);
                SU2_DTR(settings, su2_conf, "DMD", idtr, itr);
//                Write_History_ResError(settings, "DMD", idtr, itr);
                std::cout << std::endl;
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
        std::cout << "-------Performin RDMD ResEval---------" << std::endl ;
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
        int N_notZero;


        for (int i = 0; i < nC; i++) {
            Phi[i] = Eigen::MatrixXd::Zero(Nr,settings.Ns);
            COEF[i] = Eigen::MatrixXd::Zero(settings.Ns, settings.Ns);
            lambda[i] = Eigen::VectorXd::Zero(settings.Ns);
        }

        if ( settings.flag_wdb_be == "WRITE" || settings.flag_wdb_be == "NO"){
            std::cout << std::endl << "Extraction of the basis" << std::endl << std::endl;
            for ( int ncons = 0; ncons < nC; ncons ++ ) {

                std::cout << "Processing conservative variable " << ncons << std::endl;
                Phi[ncons] = RDMD_modes_coefs(sn_set.middleRows(ncons * Nr, Nr),
                                              Coefs,
                                              lambda[ncons],
                                              K_pc,
                                              -1, //Performing singular value hard threshold for DMD reduction at each step
                                              settings.r_RDMD,
                                              settings.En);
//            Phi[ncons] = RDMD_lsq_basis(sn_set.middleRows(ncons * Nr, Nr),
//                                              Coefs,
//                                              lambda[ncons],
//                                              K_pc,
//                                              -1, //Performing singular value hard threshold for DMD reduction at each step
//                                              settings.r_RDMD,
//                                              settings.En);
                COEF[ncons] = Coefs.transpose();
                surr_coefs[ncons] = getSurrCoefs(t_vec,
                                                 Coefs.transpose(),
                                                 settings.flag_interp);

                N_notZero = Phi[ncons].cols();
                if (settings.r == 0) Nm = Nmod(settings.En, K_pc);
                else Nm = std::min(settings.r, N_notZero);
                std::cout << "Number of modes used in reconstruction " << Nm << std::endl;

            }

            if ( settings.flag_wdb_be == "WRITE" ) ModesDB_Write( Phi, COEF, settings);

        } else if ( settings.flag_wdb_be == "READ") {
            std::cout << "Reading RDMD basis" << std::endl;
            Nm = settings.r;
            ModeDB_Read("ModesRDMD_", "CoefsRDMD_", Phi, COEF, settings);
            std::cout << "Number of Modes to use in Reconstruction: " << Nm << std::endl;

            for ( int icons = 0; icons < nC; icons++ ) {
                surr_coefs[icons] = getSurrCoefs(t_vec,
                                                 COEF[icons],
                                                 settings.flag_interp);
            }

        }


        for ( int idtr = 0; idtr < settings.Dt_res.size(); idtr++ ) {
            std::cout << " --------------DT_RES = " << settings.Dt_res[idtr] << "--------------"<< std::endl;
            for (int itr = 0; itr < settings.t_res.size(); itr++) {
                std::cout << " Computing residuals at time = " << settings.t_res[itr] << std::endl;

                for (int ncons = 0; ncons < nC; ncons++) {

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

                    Sn_Cons_time.middleRows(ncons * Nr, Nr) = Phi[ncons].leftCols(Nm) * coef_t.transpose();

                }

                if (settings.flag_mean == "IC") {
                    for (int it = 0; it < 3; it++)
                        Sn_Cons_time.col(it) += Ic;
                }

                //Launching SU2_DTR and saving errors and Residuals to file
                int iter = std::round(settings.t_res[itr]/settings.Dt_cfd);
                Write_Restart_Cons_Time(Sn_Cons_time, Coords, settings.out_file, iter, nC, settings.alpha, settings.beta, binary);
                SU2_DTR(settings, su2_conf, "RDMD", idtr, itr);
//                Write_History_ResError(settings, "RDMD", idtr, itr);
                std::cout << std::endl;

            }

        }
    }


    //Defining scope for GPOD
//    std::vector<std::string>::iterator itGPOD;
//    itGPOD = std::find (methods.begin(), methods.end(), "GPOD");
//    if (itGPOD != methods.end()) {
//        std::cout << "--------------------------------------" << std::endl ;
//        std::cout << "---Performing Gradient POD ResEval----" << std::endl ;
//        std::cout << "--------------------------------------" << std::endl ;
//        //Vector of MatrixXd where to store the evolution in time of conservative variables
//        Eigen::MatrixXd Sn_Cons_time = Eigen::MatrixXd::Zero(nC*Nr, 3);
//        std::vector<Eigen::MatrixXd> Phi(nC);
//        std::vector<Eigen::VectorXd> lambda(nC);
//
//        std::vector< std::vector<rbf> > surr_coefs(nC);
//        Eigen::VectorXd K_pc(settings.Ns);
//        Eigen::MatrixXd Coeffs = Eigen::MatrixXd::Zero(settings.Ns, settings.Ns);
//        std::vector<int> Nm(nC);
//        int N_notZero;
//        //Check only for POD for now
//        for (int i = 0; i < nC; i++) {
//            Phi[i] = Eigen::MatrixXd::Zero(Nr,settings.Ns);
//            lambda[i] = Eigen::VectorXd::Zero(settings.Ns);
//        }
//        std::cout << std::endl << "Extraction of the basis" << std::endl << std::endl;
//        for ( int ncons = 0; ncons < nC; ncons ++ ) {
//            std::cout << "Processing conservative variable " << ncons << std::endl;
//            Phi[ncons] = GPOD_basis( settings.Dt_cfd*settings.Ds,
//                                 sn_set.middleRows(ncons * Nr, Nr),
//                                      lambda[ncons], Coeffs, settings.r );
//            N_notZero = Phi[ncons].cols();
//            if (settings.r == 0) Nm[ncons] = Nmod(settings.En, K_pc);
//            else Nm[ncons] = std::min(settings.r, N_notZero);
//            std::cout << "Size of matrix Coeffs : [" << Coeffs.rows() << "; " << Coeffs.cols() << "]" << std::endl;
//            std::cout << "Number of modes used in reconstruction " << Nm[ncons] << std::endl;
//            surr_coefs[ncons] = getSurrCoefs(t_vec, Coeffs.transpose(), settings.flag_interp);
//        }
//        std::cout << std::endl;
//
////        std::cout << "Computing SPOD " << Nf[nfj] << " reconstruction for each conservative variable ... " << "\n";
//
//        for ( int idtr = 0; idtr < settings.Dt_res.size(); idtr++ ) {
//            std::cout << " --------------DT_RES = " << settings.Dt_res[idtr] << "--------------"<< std::endl;
//
//            for (int itr = 0; itr < settings.t_res.size(); itr++) {
//                std::cout << "Computing residuals at time t = " << settings.t_res[itr] << std::endl;
//                Modify_su2_cfg ( su2_conf, su2_conf_new, settings, idtr, itr, V_inf );
//                for (int ncons = 0; ncons < nC; ncons++) {
//
//                    Eigen::MatrixXd coef_t(3, Nm[ncons]);
//                    std::vector<double> tr(1);
//                    std::vector<double> t_evaluate = {settings.t_res[itr] - 2.0 * settings.Dt_res[idtr],
//                                                      settings.t_res[itr] - settings.Dt_res[idtr],
//                                                      settings.t_res[itr]};
//                    for (int j = 0; j < 3; j++) {
//                        tr[0] = t_evaluate[j];
//                        for (int i = 0; i < Nm[ncons]; i++)
//                            surr_coefs[ncons][i].evaluate(tr, coef_t(j, i));
//                    }
//
//                    //    }
//                    Sn_Cons_time.middleRows(ncons * Nr, Nr) = Phi[ncons].leftCols(Nm[ncons]) * coef_t.transpose();
//                }
//
//                if (settings.flag_mean == "IC") {
//                    for (int it = 0; it < 3; it++)
//                        Sn_Cons_time.col(it) += Ic;
//                }
//
//                std::string mv_string;
//                if (settings.Ns == settings.r)
//                    mv_string = "mv history_rbm_00002.csv history_gpod_AllModes_" +
//                                std::to_string(settings.Dt_res[idtr]) + "_" + std::to_string(itr) + ".csv";
//                else
//                    mv_string = "mv history_rbm_00002.csv history_gpod_" +
//                                std::to_string(Nm[0]) + "_" + std::to_string(settings.Dt_res[idtr]) + "_" + std::to_string(itr) +
//                                ".csv";
//
//                len_s = mv_string.length();
//                char mv_sys_call[len_s + 10];
//                strcpy(mv_sys_call, mv_string.c_str());
//                // Write_Restart_Cons_Time( Sn_Cons_time, Coords, settings.out_file, t_evaluate.size(), nC, alpha );
//                Write_Restart_Cons_Time(Sn_Cons_time, Coords, settings.out_file, 3, nC, settings.alpha, settings.beta, binary);
//                //Executing SU2, removing all useless files, renaming files with residuals
//                std::cout << "Calling SU2 for residual evaluation and writing file to history " << std::endl;
//                auto opt = std::system(su2_sys_call);
//                opt = std::system(rmf_sys_call);
//                opt = std::system(mv_sys_call);
//
//                std::cout << std::endl;
//            }
//        }
//    }
//
//    //Defining scope for DMD adaptive sampling
//    std::vector<std::string>::iterator itDMD_adaptsample;
//    itDMD_adaptsample = std::find (methods.begin(), methods.end(), "DMD-AS");
//    if (itDMD_adaptsample != methods.end()) {
//
//        int nVar = settings.t_pos.size();
//
//        if ( nVar != settings.r ) std::cout << "-----------------WARNING-----------------\n "
//                                       "number of samples different from rank, resetting to number of samples!" << std::endl << std::endl;
//
//        int Nm = nVar;
//        Eigen::VectorXi Ipos = Eigen::VectorXi::Zero(nVar);
//        for ( int ipos = 0; ipos < nVar; ipos++ ) Ipos(ipos) = settings.t_pos[ipos];
//
//
//        //Vector of MatrixXd where to store the evolution in time of conservative variables
//        Eigen::MatrixXd Sn_Cons_time = Eigen::MatrixXd::Zero(nC*Nr, 3);
//        Eigen::MatrixXd Phi_POD = Eigen::MatrixXd::Zero(Nr,nVar);
//        Eigen::VectorXd lambda_POD = Eigen::VectorXd::Zero(nVar);
//        Eigen::MatrixXcd Phi_DMD = Eigen::MatrixXd::Zero(Nr,settings.Ns);
//        Eigen::VectorXcd lambda_DMD = Eigen::VectorXd::Zero(settings.Ns);
//        Eigen::VectorXd K_pc(nVar);
//        Eigen::MatrixXd eig_vec(nVar, nVar);
//
//        std::vector<smartuq::surrogate::rbf> surr_coefs_POD ;
//        std::vector<smartuq::surrogate::rbf> surr_coefs_DMD_r;
//        std::vector<smartuq::surrogate::rbf> surr_coefs_DMD_i;
//
//        int N_notZero;
//
//        Eigen::MatrixXd sub_sn_set = indexing(sn_set, Eigen::ArrayXi::LinSpaced(nC*Nr,0,nC*Nr-1),Ipos);
//
//        Eigen::MatrixXcd eig_vec_DMD;
//        std::cout << "Computing adaptive DMD modes" << std::endl;
////            if (settings.r == 0) Nm = Nmod(settings.En, K_pc);
////            else Nm = std::min(settings.r, settings.Ns);
//        Phi_DMD = DMD_Adaptive_basis(sn_set,
//                                     lambda_DMD,
//                                     eig_vec_DMD,
//                                     lambda_POD,
//                                     eig_vec,
//                                     Ipos);
//
//        Eigen::MatrixXcd PhiTPhi = Phi_DMD.transpose()*Phi_DMD;
//        Eigen::MatrixXcd Coeffs = PhiTPhi.inverse()*(Phi_DMD.transpose()*sn_set);
//        surr_coefs_DMD_r = getSurrCoefs(t_vec, Coeffs.real().transpose(), settings.flag_interp);
//        surr_coefs_DMD_i = getSurrCoefs(t_vec, Coeffs.imag().transpose(), settings.flag_interp);
//
//        for ( int idtr = 0; idtr < settings.Dt_res.size(); idtr++ ) {
//            std::cout << " --------------DT_RES = " << settings.Dt_res[idtr] << "--------------" << std::endl;
//
//            for (int itr = 0; itr < settings.t_res.size(); itr++) {
//                std::cout << " Computing residuals at time : " << settings.t_res[itr] << std::endl;
//                Modify_su2_cfg ( su2_conf, su2_conf_new, settings, idtr, itr, V_inf );
//
//                    Eigen::MatrixXcd coef_t(3, Nm);
//
//                    std::vector<double> tr(1);
//                    std::vector<double> t_evaluate = {settings.t_res[itr] - 2.0 * settings.Dt_res[0],
//                                                      settings.t_res[itr] - settings.Dt_res[0],
//                                                      settings.t_res[itr]};
//
//                    double tmp_r, tmp_i;
//                    for (int j = 0; j < 3; j++) {
//                        tr[0] = t_evaluate[j];
//                        for (int i = 0; i < Nm; i++) {
//                            surr_coefs_DMD_r[i].evaluate(tr, tmp_r);
//                            surr_coefs_DMD_i[i].evaluate(tr, tmp_i);
//                            std::complex<double> c(tmp_r, tmp_i);
//                            coef_t(j, i) = c;
//                        }
//                    }
//
//                    Eigen::MatrixXcd Appo = Phi_DMD * coef_t.transpose();
//                    Sn_Cons_time = Appo.real();
//
//
//                //Introduce an if on the number of conservative variables
////                Sn_Cons_time.middleRows(0, Nr) =
////                        Sn_Cons_time.middleRows(0, Nr) * (rho_max - rho_min) + Eigen::MatrixXd::Ones(Nr, 3) * rho_min;
////                Sn_Cons_time.middleRows(Nr, Nr) = Sn_Cons_time.middleRows(Nr, Nr) * (rhoU_max - rhoU_min) +
////                                                  Eigen::MatrixXd::Ones(Nr, 3) * rhoU_min;
////                Sn_Cons_time.middleRows(2 * Nr, Nr) = Sn_Cons_time.middleRows(2 * Nr, Nr) * (rhoV_max - rhoV_min) +
////                                                      Eigen::MatrixXd::Ones(Nr, 3) * rhoV_min;
////                Sn_Cons_time.middleRows(3 * Nr, Nr) = Sn_Cons_time.middleRows(3 * Nr, Nr) * (rhoE_max - rhoE_min) +
////                                                      Eigen::MatrixXd::Ones(Nr, 3) * rhoE_min;
//
//                if (settings.flag_mean == "IC") {
//                    for (int it = 0; it < 3; it++)
//                        Sn_Cons_time.col(it) += Ic;
//                } else if (settings.flag_mean == "YES") {
//                    for (int it = 0; it < 3; it++)
//                        Sn_Cons_time.col(it) += mean;
//                }
//
//                std::string mv_string = "mv history_rbm_00002.csv history_dmd_adapt_" +
//                            std::to_string(Nm) + "_" + std::to_string(settings.Dt_res[idtr]) + "_" + std::to_string(itr) +
//                            ".csv";
//
//                len_s = mv_string.length();
//                char mv_sys_call[len_s + 1];
//                strcpy(mv_sys_call, mv_string.c_str());
//                std::cout << "Writing time reconstruction " << std::endl;
//                // Write_Restart_Cons_Time( Sn_Cons_time, Coords, settings.out_file, t_evaluate.size(), nC, alpha );
//                Write_Restart_Cons_Time(Sn_Cons_time, Coords, settings.out_file, 3, nC, settings.alpha, settings.beta, binary);
//                //Executing SU2, removing all useless files, renaming files with residuals
//                std::cout << "Calling SU2 for residual evaluation and writing file to history " << std::endl;
//                auto otp = std::system(su2_sys_call);
//                otp = std::system(rmf_sys_call);
//                otp = std::system(mv_sys_call);
//            }
//        }
//    }


    std::cout << "MODES Adaptive Residual Evaluation ends" << std::endl << std::endl;
    return 0;

}
