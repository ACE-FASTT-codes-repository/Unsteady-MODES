/*
CODE FOR ERROR EVALUATION OF CONSERVATIVE QUANTITIES USING DIFFERENT RBM TECHNIQUES
!!!!!!(restart files need to be numbered without gaps)!!!!
INPUT ARGUMENTS
Config File RBM
ATTENTION ------> Set DS to the minimum possible in config file
 */

#include "Extract_Basis.hpp"
#include "read_Inputs.hpp"
#include "Pre-Process.hpp"
#include "System_Calls.hpp"
#include "Post-Process.hpp"
#include "Generate_snset.hpp"
#include "Reconstruction.hpp"
#include "write_Outputs.hpp"



int main( int argc, char *argv[] ) {

    std::cout << "-----------Pareto Front with direct error starts-----------" << std::endl << std::endl;

    std::cout << "Initializing common variables ... " << std::endl << std::endl;
    prob_settings settings;
    std::string filecfg = argv[1];

    std::vector<Eigen::VectorXd> Err_RBM_Nm_time;
    std::vector<Eigen::VectorXd> ErrP_RBM_Nm_time;
    std::vector<Eigen::VectorXd> EN;

    Read_cfg(filecfg, settings);


    int nC = settings.Cols.size();
    int Nr;
    Eigen::MatrixXd Coords;
    std::vector<double> t_vec(settings.Ns);
    common_vars(Nr, Coords, t_vec, settings);

    std::cout << "Storing snapshot Matrix ... \n ";
    Eigen::MatrixXd sn_set = generate_snap_matrix(Nr, settings);

    int Nsnap;
    if (settings.Ns % 2 == 0) {
        Nsnap = settings.Ns/2;
        t_vec.pop_back();
    }
    else Nsnap = std::floor(settings.Ns / 2) + 1;

    std::vector<double> t_train(Nsnap);

    t_train[0] = settings.nstart * settings.Dt_cfd;
    for (int i = 1; i < t_train.size(); i++)
        t_train[i] = t_train[i - 1] + settings.Dt_cfd * (double) settings.Ds * 2.0;// + settings.Dt_cfd*(double)settings.Ds/2.0;


    std::cout << "Computing mean/Initial Condition of CFD solution ... " << std::endl;
    //Defining Initial condition
    Eigen::VectorXd Ic = Eigen::VectorXd::Zero(nC * Nr);

    Eigen::VectorXd svd_cum_sum(settings.Ns);

    if (settings.flag_mean == "IC") {
        std::cout << "Subtracting Initial condition" << std::endl << std::endl;
        Ic = IC(sn_set, settings, nC, Nr);
    } else {
        std::cout << "Using data without subtracting any reference state" << std::endl << std::endl;
    }

    Eigen::MatrixXd sn_set_p(nC * Nr, Nsnap);
    int count = 0;
    for (int is = 0; is < settings.Ns; is += 2) {
        sn_set_p.col(count) = sn_set.col(is);
        count++;
    }

    Eigen::MatrixXd norm_sn_set = Eigen::MatrixXd::Zero(nC, t_vec.size());
    for (int j = 0; j < nC; j++) {
        for (int i = 0; i < t_vec.size(); i++)
            norm_sn_set(j, i) = sn_set.middleRows(j * Nr, Nr).col(i).norm();

    }

//Defining common scope for POD
    auto methods = settings.flag_method;
    std::vector<std::string>::iterator itPOD;
    itPOD = std::find(methods.begin(), methods.end(), "POD");
    if (itPOD != methods.end()) {

        Eigen::VectorXd lambda(Nsnap);
        Eigen::VectorXd K_pc(Nsnap);
        Eigen::MatrixXd eig_vec(Nsnap, Nsnap);
        std::vector<Eigen::MatrixXd> Phi(nC);
        std::vector<Eigen::MatrixXd> Coeffs_int(nC);
        std::vector<Eigen::VectorXd> eigs(nC);
        int Nm;
        std::vector<int> N_notZero(nC);
        //Check only for POD for now

        for (int ncons = 0; ncons < nC; ncons++) {
            std::cout << "Processing conservative variable " << ncons << std::endl;
            Phi[ncons] = SPOD_basis(sn_set_p.middleRows(ncons * Nr, Nr),
                                             lambda, K_pc, eig_vec,
                                             0,
                                             settings.flag_bc,
                                             settings.flag_filter,
                                             settings.sigma);
            N_notZero[ncons] = Phi[ncons].cols();

            std::vector<rbf> surr_coefs = getSurrCoefs(t_train, eig_vec, settings.flag_interp);
            Eigen::MatrixXd coef_t(t_vec.size(), N_notZero[ncons]);

            std::vector<double> tr(1);
            for (int j = 0; j < t_vec.size(); j++) {
                tr[0] = t_vec[j];
                for (int i = 0; i < N_notZero[ncons]; i++)
                    surr_coefs[i].evaluate(tr, coef_t(j, i));
            }

            Coeffs_int[ncons] = coef_t;
            eigs[ncons] = lambda;

        }

        for (int imode = 2; imode < sn_set_p.cols() + 1; imode++) {

            std::cout << "For iMode  " << imode << std::endl;
            std::cout << "Computing error of interpolation..." << "\t";

            for ( int icons = 0; icons < nC; icons++ ) {

                Nm = std::min(imode,N_notZero[icons]);
                Eigen::MatrixXd Sig = Eigen::MatrixXd::Zero(Nm, Nm);
                for (int i = 0; i < Nm; i++)
                    Sig(i, i) = std::sqrt(eigs[icons](i));

                Eigen::MatrixXd Err_SPOD_map = Eigen::MatrixXd::Zero(Nr, t_vec.size());
                Eigen::VectorXd Err_SPOD_Nm_time = Eigen::VectorXd::Zero(t_vec.size());

                Err_SPOD_map = sn_set.middleRows(icons * Nr, Nr).leftCols(t_vec.size()) -
                               Phi[icons].leftCols(Nm) * Sig * Coeffs_int[icons].leftCols(Nm).transpose();

                for (int i = 0; i < t_vec.size(); i++) {
                    for (int j = 0; j < Nr; j++)
                        Err_SPOD_Nm_time(i) += Err_SPOD_map(j, i) * Err_SPOD_map(j, i);

                    Err_SPOD_Nm_time(i) = std::sqrt(Err_SPOD_Nm_time(i)) / norm_sn_set(icons, i);

                }

                Err_RBM_Nm_time.push_back(Err_SPOD_Nm_time);
                EN.push_back(K_pc);

            }

            std::ofstream errfile;
            std::string file_err_name = "Err_POD_" + std::to_string(imode) + ".dat";
            errfile.open(file_err_name);

            for (int nm = 0; nm < t_vec.size(); nm++) {
                for (int j = 0; j < Err_RBM_Nm_time.size(); j++)
                    errfile << std::setprecision(8) << Err_RBM_Nm_time[j](nm) << "\t";

                errfile << std::endl;

            }

            errfile.close();

            Err_RBM_Nm_time.clear();
            ErrP_RBM_Nm_time.clear();
            EN.clear();
        }

    }

    //Defining common scope for RDMD
    std::vector<std::string>::iterator itRDMD;
    itRDMD = std::find(methods.begin(), methods.end(), "RDMD");
    if (itRDMD != methods.end()) {

        Eigen::VectorXd lambda =  Eigen::VectorXd::Zero(Nsnap);
        Eigen::VectorXd K_pc(Nsnap);
        Eigen::MatrixXd Coefs = Eigen::MatrixXd::Zero(settings.Ns, settings.Ns);
        std::vector<Eigen::MatrixXd> Phi(nC);
        std::vector<Eigen::MatrixXd> Coeffs_int(nC);

        int Nm;
        //Check only for POD for now

        for (int ncons = 0; ncons < nC; ncons++) {
            std::cout << "Processing conservative variable " << ncons << std::endl;
            Phi[ncons] = RDMD_modes_coefs ( sn_set_p.middleRows(ncons*Nr,Nr),
                    Coefs,
                    lambda,
                    K_pc,
                    -1,
                    Nsnap,
                    settings.En );

            std::vector<rbf> surr_coefs = getSurrCoefs(t_train, Coefs.transpose(), settings.flag_interp);
            Eigen::MatrixXd coef_t(t_vec.size(), Nsnap);

            std::vector<double> tr(1);
            for (int j = 0; j < t_vec.size(); j++) {
                tr[0] = t_vec[j];
                for (int i = 0; i < Nsnap; i++)
                    surr_coefs[i].evaluate(tr, coef_t(j, i));
            }

            Coeffs_int[ncons] = coef_t;

        }

        for (int imode = 2; imode < sn_set_p.cols() + 1; imode++) {

            std::cout << "For iMode  " << imode << std::endl;
            std::cout << "Computing error of interpolation..." << std::endl;

            for ( int icons = 0; icons < nC; icons++ ) {

                Nm = imode;

                Eigen::MatrixXd Err_RDMD_map = Eigen::MatrixXd::Zero(Nr, t_vec.size());
                Eigen::VectorXd Err_RDMD_Nm_time = Eigen::VectorXd::Zero(t_vec.size());

                Err_RDMD_map = sn_set.middleRows(icons * Nr, Nr).leftCols(t_vec.size()) -
                               Phi[icons].leftCols(Nm) * Coeffs_int[icons].leftCols(Nm).transpose();

                for (int i = 0; i < t_vec.size(); i++) {
                    for (int j = 0; j < Nr; j++)
                        Err_RDMD_Nm_time(i) += Err_RDMD_map(j, i) * Err_RDMD_map(j, i);

                    Err_RDMD_Nm_time(i) = std::sqrt(Err_RDMD_Nm_time(i)) / norm_sn_set(icons, i);

                }

                Err_RBM_Nm_time.push_back(Err_RDMD_Nm_time);
                EN.push_back(K_pc);

            }

            std::ofstream errfile;
            std::string file_err_name = "Err_RDMD_" + std::to_string(imode) + ".dat";
            errfile.open(file_err_name);

            for (int nm = 0; nm < t_vec.size(); nm++) {
                for (int j = 0; j < Err_RBM_Nm_time.size(); j++)
                    errfile << std::setprecision(8) << Err_RBM_Nm_time[j](nm) << "\t";

                errfile << std::endl;

            }

            errfile.close();

            Err_RBM_Nm_time.clear();
            ErrP_RBM_Nm_time.clear();
            EN.clear();
        }

    }

    //Defining common scope for DMD
    std::vector<std::string>::iterator itDMD;
    itDMD = std::find(methods.begin(), methods.end(), "DMD");
    if (itDMD != methods.end()) {

        Eigen::VectorXd lambda_POD;
        Eigen::MatrixXd eig_vec_POD;
        Eigen::VectorXcd lambda_DMD;
        Eigen::MatrixXcd eig_vec_DMD;

        for (int imode = 2; imode < sn_set_p.cols() ; imode++) {
            std::cout << "For iMode " << imode << std::endl;

            for (int ncons = 0; ncons < nC; ncons++) {

                Eigen::MatrixXcd Phi;
                Eigen::VectorXcd alfa;

                Phi = DMD_basis( sn_set_p.middleRows(ncons*Nr,Nr),
                                 lambda_DMD,
                                 eig_vec_DMD,
                                 lambda_POD,
                                 eig_vec_POD,
                                 imode );

                //         int Nm = Phi.cols();
                //         std::cout << "Number of modes extracted : " << Nm << std::endl;

                Eigen::VectorXcd omega(Phi.cols());
                for ( int i = 0; i < Phi.cols(); i++ )
                    omega(i) = std::log(lambda_DMD(i))/(settings.Dt_cfd*settings.Ds*2.0);

                // std::cout << "Calculating coefficients DMD ... " << "\t";
                alfa = Calculate_Coefs_DMD_exact (sn_set_p.middleRows(ncons*Nr,Nr).leftCols(Nsnap-1),
                                                   lambda_DMD,
                                                   Phi );


                Eigen::MatrixXcd V_and(lambda_DMD.size(), t_vec.size());
                for ( int i = 0; i < lambda_DMD.size(); i++ ) {
                    for ( int j = 0; j < t_vec.size(); j++ )
                        V_and(i,j) = std::pow(lambda_DMD(i), (double)j/2.0);
                }

                Eigen::MatrixXcd Psi = Eigen::MatrixXcd::Zero(alfa.size(), t_vec.size());
                for ( int i = 0; i < t_vec.size(); i++ )
                    Psi.col(i) = alfa.cwiseProduct(V_and.col(i));

                Eigen::MatrixXcd D_dmd = Phi*Psi;
                Eigen::MatrixXd Err_DMD_map = Eigen::MatrixXd::Zero(Nr, t_vec.size());
                Eigen::VectorXd Err_DMD_Nm_time = Eigen::VectorXd::Zero(t_vec.size());

                Err_DMD_map = sn_set.middleRows(ncons*Nr,Nr).leftCols(t_vec.size()) - D_dmd.real();

                for ( int i = 0; i < t_vec.size(); i++ ) {

                    for ( int j = 0; j < Nr; j++ ) {
                        Err_DMD_Nm_time(i) += Err_DMD_map(j,i)*Err_DMD_map(j,i);
                    }
                    Err_DMD_Nm_time(i) = std::sqrt(Err_DMD_Nm_time(i))/norm_sn_set(ncons,i);
                }
                Err_RBM_Nm_time.push_back(Err_DMD_Nm_time);
            }

            std::ofstream errfile;
            std::string file_err_name = "Err_DMD_" + std::to_string(imode) + ".dat";
            errfile.open(file_err_name);

            for (int nm = 0; nm < t_vec.size(); nm++) {
                for (int j = 0; j < Err_RBM_Nm_time.size(); j++)
                    errfile << std::setprecision(8) << Err_RBM_Nm_time[j](nm) << "\t";

                errfile << std::endl;

            }

            errfile.close();
            Err_RBM_Nm_time.clear();

        }
    }


    std::cout << "-----------Pareto Front with direct error ends-----------" << std::endl << std::endl;
    return 0;
}
