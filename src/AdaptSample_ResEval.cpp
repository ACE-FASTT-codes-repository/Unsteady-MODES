/*
CODE FOR RESIDUAL EVALUATION OF DIFFERENT RBM TECHNIQUES USING SU2 CODE
INPUT ARGUMENTS
Config File RBM + Config File SU2
*/

#include "Extract_Basis.hpp"
#include "read_Inputs.hpp"
#include "Pre-Process.hpp"
#include "Generate_snset.hpp"
#include "Reconstruction.hpp"
#include "write_Outputs.hpp"


int main( int argc, char *argv[] )
{

    std::cout << "-----------MODES Adaptive Residual Evaluation starts-------------" << std::endl << std::endl;
//
    std::cout << "Initializing common variables ... " << std::endl << std::endl;
    prob_settings settings;
    std::string filecfg = argv[1];
    std::string su2_conf = argv[2];
    std::string root_conf;
    root_conf.assign ( su2_conf, 0, su2_conf.size() - 4);
    std::string su2_conf_new = root_conf + "-reseval.cfg";
    std::string decision = argv[3];
    std::string su2dtr_string = "mpirun -np 6 ./SU2_DTR " + su2_conf_new + " > SU2.log"; // + " > resEval_su2.log";
    int len_s = su2dtr_string.length();
    char su2_sys_call[len_s + 1];
    strcpy(su2_sys_call, su2dtr_string.c_str());

    Read_cfg( filecfg, settings );
//
    std::string root_outputfile;
    root_outputfile.assign ( settings.out_file, 0, settings.out_file.size() - 4);
    std::string root_inputfile;
    root_inputfile.assign ( settings.in_file, 0, settings.in_file.size() - 4);
    std::string input_format;
    input_format.assign ( settings.in_file, settings.in_file.size() - 3, 3);
    std::string rmf_string = "rm -f " + root_outputfile + "_*";
    len_s = rmf_string.length();
    char rmf_sys_call[len_s + 1];
    strcpy(rmf_sys_call, rmf_string.c_str());

    int s_Nf = 1;   //Number of values for the SPOD filter (POD included)
    std::vector<int> Nf(s_Nf);
    Nf[0] = 0;
//    Nf[1] = std::ceil(settings.Ns/10.0);
//    Nf[2] = std::ceil(settings.Ns/2.0);
//    Nf[3] = std::ceil(2.0*settings.Ns/3.0);
//    Nf[4] = settings.Ns;

    int nC = settings.Cols.size();

    std::stringstream buffer;
    buffer << std::setfill('0') << std::setw(5) << std::to_string(settings.nstart);
    std::string file_1 = root_inputfile + "_" + buffer.str() + "." + input_format;

    std::cout << "Initializing Vector of times ... " << std::endl;
    std::vector<double> t_vec( settings.Ns );
    t_vec[0] = settings.nstart*settings.Dt_cfd;
    for ( int i = 1; i < settings.Ns; i++ )
        t_vec[i] = t_vec[i-1] + settings.Dt_cfd*settings.Ds;

    std::cout << t_vec[0] << "-------time of first snapshot" << std::endl;
    std::cout << t_vec[t_vec.size()-1] <<"-------time of last snapshot" << std::endl;
    std::cout << settings.t_res[0] <<"----time of first residual evaluation"<< std::endl;
    std::cout << settings.t_res[settings.t_res.size()-1]<< "----time of last residual evaluation"<< std::endl;

    for ( int idtr = 0; idtr < settings.Dt_res.size(); idtr++ ){
        std::cout<< "for this DT_res=" << settings.Dt_res[idtr]<<"......."<< std::endl;
        if ( ((settings.t_res[0] - (2.0 * settings.Dt_res[idtr])) < t_vec[0]) ||  (settings.t_res[settings.t_res.size()-1] > t_vec[t_vec.size() - 1]))
        {
            std::cout
                    << "Define proper Delta_t_res and T_RES vector " << std::endl;
            exit(EXIT_FAILURE);
        }else
        { std::cout << "Perfect usage of time... chapeau "<< std::endl;}
    }

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

    // std::vector<double> t_evaluate(2*settings.Ns-1);
    // t_evaluate[0] = settings.nstart*settings.Dt_cfd;
    // for ( int i = 1; i < t_evaluate.size(); i++)
    //     t_evaluate[i] = t_evaluate[i-1] + settings.Dt_cfd*(double)settings.Ds/2.0;

    std::cout << "Computing mean/Initial Condition of CFD solution ... " << std::endl;

    Eigen::VectorXd mean = sn_set.rowwise().mean();
    if ( settings.flag_mean == "YES" ) {
        for ( int it = 0; it < settings.Ns; it++ )
            sn_set.col(it) -= mean;
    }

    //Defining Initial condition
    Eigen::VectorXd Ic = Eigen::VectorXd::Zero(nC*Nr);

    std::string binary = "YES";
    Eigen::VectorXd svd_cum_sum(settings.Ns);

    //sn_set is modified according to flag
    std::string flag = "YES";
    if ( settings.flag_mean == "IC" ) Ic = IC(sn_set, settings, nC, Nr);
    std::cout << "Ic size :" << Ic.size() << std::endl;

    double rho_max, rho_min, rhoU_max, rhoU_min, rhoV_max, rhoV_min, rhoW_max, rhoW_min, rhoE_max, rhoE_min,
            tke_min, tke_max, omega_min, omega_max, nuTilde_min, nuTilde_max; //add turbulence

    get_MinMax_ConsVar(sn_set,settings,nC,rho_max, rho_min, rhoU_max, rhoU_min, rhoV_max, rhoV_min,
                       rhoW_max, rhoW_min, rhoE_max, rhoE_min,tke_min, tke_max, omega_min,
                       omega_max, nuTilde_min, nuTilde_max);

    if ( flag == "YES") Direct_Normalization(sn_set,settings,nC,rho_max, rho_min, rhoU_max, rhoU_min, rhoV_max, rhoV_min,
                                             rhoW_max, rhoW_min, rhoE_max, rhoE_min,tke_min, tke_max, omega_min,
                                             omega_max, nuTilde_min, nuTilde_max);


    //We start trying like this
    int Nsamples = settings.Ns;
    Eigen::VectorXi Ipos = Inverse_POS(sn_set, Nsamples);

    int Nm = Ipos.size(); //settings.t_pos.size();
    int nVar = Ipos.size();//settings.t_pos.size(); //that has to contain also first and last snapshot


    std::cout << "Positions sampled in time \n" << Ipos.transpose() << std::endl;
    std::cout << "From N = " << Nsamples << " uniform samples we get Nm = " << Nm << "non-uniform samples" << std::endl;

    //If gust is active we need free stream velocity
    double M = settings.Mach;
    double T = settings.T;
    double alpha = M_PI*settings.alpha/double(180);
    double beta = M_PI*settings.beta/double(180);
    if (settings.ndim == 2) beta = 0.0;
    double R = 287.058;
    double gamma = 1.4;

    double V_inf = M*std::sqrt(gamma*R*T)*std::cos(alpha);

    //Defining common scope for uniform sampling
    {

        //Vector of MatrixXd where to store the evolution in time of conservative variables
        Eigen::MatrixXd Sn_Cons_time = Eigen::MatrixXd::Zero(nC * Nr, 3);
        Eigen::MatrixXd Phi = Eigen::MatrixXd::Zero(Nr, settings.Ns);
        Eigen::VectorXd lambda_POD = Eigen::VectorXd::Zero(settings.Ns);
        Eigen::MatrixXcd Phi_DMD = Eigen::MatrixXd::Zero(Nr, settings.Ns);
        Eigen::VectorXcd lambda_DMD = Eigen::VectorXd::Zero(settings.Ns);
        Eigen::VectorXd K_pc(settings.Ns);
        Eigen::MatrixXd eig_vec(settings.Ns, settings.Ns);

        //To define for interp reconstruction formula
        std::vector<smartuq::surrogate::rbf> surr_coefs_POD;
        std::vector<smartuq::surrogate::rbf> surr_coefs_RDMD;
        std::vector<smartuq::surrogate::rbf> surr_coefs_DMD_r;
        std::vector<smartuq::surrogate::rbf> surr_coefs_DMD_i;
        //To define for exponential rec
        Eigen::VectorXcd alfa(Nm);

        int N_notZero;
        //Check only for POD for now

        if ( settings.flag_method[0] == "SPOD" ){

            std::cout << "Computing uniform SPOD modes with Nf : " << settings.Nf[0]; << "\n";
            Phi = SPOD_basis(sn_set,
                                 lambda_POD, K_pc, eig_vec,
                                 settings.Nf[0],
                                 settings.flag_bc,
                                 settings.flag_filter,
                                 settings.sigma);
            N_notZero = Phi.cols();
//            if (settings.r == 0) Nm = Nmod(settings.En, K_pc);
//            else Nm = std::min(settings.r, N_notZero);

            if ( Nm != nVar) {
                std::cout << "-----------------------------WARNING:------------------------"
                             "Truncation and adaptation don't have the same number of modes" << std::endl << std::endl;
                exit(EXIT_FAILURE);
            }
            std::cout << "Number of modes used in reconstruction " << Nm << std::endl;
            surr_coefs_POD = getSurrCoefs(t_vec, eig_vec, settings.flag_interp);

        } else if ( settings.flag_method[0] == "DMD" ) {

            Eigen::MatrixXcd eig_vec_DMD;
            std::cout << "Computing uniform DMD modes" << std::endl;
            Phi_DMD = DMD_basis(sn_set,
                                   lambda_DMD,
                                   eig_vec_DMD,
                                   lambda_POD,
                                   eig_vec,
                                   nVar);

            Eigen::MatrixXcd PhiTPhi = Phi_DMD.leftCols(Nm).transpose()*Phi_DMD.leftCols(Nm);
            Eigen::MatrixXcd Coeffs = PhiTPhi.inverse()*(Phi_DMD.leftCols(Nm).transpose()*sn_set);
            surr_coefs_DMD_r = getSurrCoefs(t_vec, Coeffs.real().transpose(), settings.flag_interp);
            surr_coefs_DMD_i = getSurrCoefs(t_vec, Coeffs.imag().transpose(), settings.flag_interp);
            alfa = Calculate_Coefs_DMD_exact(sn_set.leftCols(settings.Ns - 1),
                                                                  lambda_DMD,
                                                                  Phi_DMD);

        } else if ( settings.flag_method[0] == "RDMD" ){

            std::cout << "Computing uniform RDMD modes \n";
            Eigen::MatrixXd Coeffs = Eigen::MatrixXd::Zero(settings.Ns,settings.Ns);
            Phi = RDMD_modes_coefs ( sn_set,Coeffs,lambda_POD,K_pc,-1,nVar,settings.En );
            N_notZero = Phi.cols();
//            if (settings.r == 0) Nm = Nmod(settings.En, K_pc);
//            else Nm = std::min(settings.r, N_notZero);

            if ( Nm != nVar) {
                std::cout << "-----------------------------WARNING:------------------------"
                             "Truncation and adaptation don't have the same number of modes" << std::endl << std::endl;
                exit(EXIT_FAILURE);
            }
            std::cout << "Number of modes used in reconstruction " << Nm << std::endl;
            surr_coefs_RDMD = getSurrCoefs(t_vec, Coeffs.transpose(), settings.flag_interp);

        }

        if ( decision == "-e") {

            for ( int itr = 0; itr < settings.t_res.size(); itr++ ) {
                std::cout << " Computing residuals at time : " << settings.t_res[itr] << std::endl;
                Modify_su2_cfg ( su2_conf, su2_conf_new, settings.Dt_res[0], settings.t_res[itr], V_inf );
                if ( settings.flag_method[0] == "SPOD" ) {

                    Eigen::MatrixXd coef_t(3, Nm);

                    std::vector<double> tr(1);
                    std::vector<double> t_evaluate = {settings.t_res[itr] - 2.0 * settings.Dt_res[0],
                                                      settings.t_res[itr] - settings.Dt_res[0],
                                                      settings.t_res[itr]};

                    for (int j = 0; j < 3; j++) {
                        tr[0] = t_evaluate[j];
                        for (int i = 0; i < Nm; i++)
                            surr_coefs_POD[i].evaluate(tr, coef_t(j, i));
                    }


                    Eigen::MatrixXd Sig = Eigen::MatrixXd::Zero(Nm, Nm);
                    for (int i = 0; i < Nm; i++)
                        Sig(i, i) = std::sqrt(lambda_POD(i));
                    Sn_Cons_time = Phi.leftCols(Nm) * Sig * coef_t.transpose();

                } else if ( settings.flag_method[0] == "DMD" ) {

                    Eigen::MatrixXcd coef_t(3, Nm);

                    std::vector<double> tr(1);
                    std::vector<double> t_evaluate = {settings.t_res[itr] - 2.0 * settings.Dt_res[0],
                                                      settings.t_res[itr] - settings.Dt_res[0],
                                                      settings.t_res[itr]};

                    if ( settings.flag_rec == "UNIFORM" ) {
                        double tmp_r, tmp_i;
                        //DMD using same rec formula as POD
//                        --------------------------------------------------------------------------------------------
                        for (int j = 0; j < 3; j++) {
                            tr[0] = t_evaluate[j];
                            for (int i = 0; i < Nm; i++) {
                                surr_coefs_DMD_r[i].evaluate(tr, tmp_r);
                                surr_coefs_DMD_i[i].evaluate(tr, tmp_i);
                                std::complex<double> c(tmp_r,tmp_i);
                                coef_t(j,i) = c;
                            }
                        }
                        Eigen::MatrixXcd Appo = Phi_DMD.leftCols(Nm) * coef_t.transpose();
                        Sn_Cons_time = Appo.real();
//                      ----------------------------------------------------------------------------------------------
                    }

                    //DMD using classical reconstruction formula
                    if ( settings.flag_rec == "EXP" ) {
                        for (int j = 0; j < 3; j++) {

                            Eigen::MatrixXcd Rec = Reconstruction_DMD(t_evaluate[j],
                                                                      settings.Dt_cfd * settings.Ds,
                                                                      alfa,
                                                                      Phi_DMD,
                                                                      lambda_DMD,
                                                                      "SCALAR");

                            Sn_Cons_time.col(j) = Rec.real();

                        }
                    }
                } else if ( settings.flag_method[0] == "RDMD" ) {

                    Eigen::MatrixXd coef_t(3, Nm);

                    std::vector<double> tr(1);
                    std::vector<double> t_evaluate = {settings.t_res[itr] - 2.0 * settings.Dt_res[0],
                                                      settings.t_res[itr] - settings.Dt_res[0],
                                                      settings.t_res[itr]};

                    for (int j = 0; j < 3; j++) {
                        tr[0] = t_evaluate[j];
                        for (int i = 0; i < Nm; i++)
                            surr_coefs_RDMD[i].evaluate(tr, coef_t(j, i));
                    }

                    Sn_Cons_time = Phi.leftCols(Nm) * coef_t.transpose();

                }

                //Introduce an if on the number of conservative variables
                if ( flag == "YES" ) {
                    Inverse_Normalization(Sn_Cons_time, settings, nC, rho_max, rho_min, rhoU_max, rhoU_min, rhoV_max,
                                          rhoV_min,rhoW_max, rhoW_min, rhoE_max, rhoE_min, tke_min, tke_max, omega_min,
                                          omega_max, nuTilde_min, nuTilde_max);
                }
                if (settings.flag_mean == "IC") {
                    for (int it = 0; it < 3; it++)
                        Sn_Cons_time.col(it) += Ic;
                } else if (settings.flag_mean == "YES") {
                    for (int it = 0; it < 3; it++)
                        Sn_Cons_time.col(it) += mean;
                }

                std::string mv_string = "mv history_rbm_00002.csv history_uniform_" + std::to_string(itr) + ".csv";
                len_s = mv_string.length();
                char mv_sys_call[len_s + 1];
                strcpy(mv_sys_call, mv_string.c_str());
                std::cout << "Writing time reconstruction " << std::endl;
                // Write_Restart_Cons_Time( Sn_Cons_time, Coords, settings.out_file, t_evaluate.size(), nC, alpha );
                Write_Restart_Cons_Time(Sn_Cons_time, Coords, settings.out_file, 3, nC, settings.alpha, settings.beta, binary);
                //Executing SU2, removing all useless files, renaming files with residuals
                std::cout << "Calling SU2 for residual evaluation and writing file to history " << std::endl;
                auto otp = std::system(su2_sys_call);
                otp = std::system(rmf_sys_call);
                otp = std::system(mv_sys_call);
            }
        } else if ( decision == "-r") {
            for ( int nt = 0; nt < settings.t_rec.size(); nt++ ) {

                std::cout << "Reconstructing Momentum at time : " << settings.t_rec[nt] << "\t";
                Eigen::MatrixXd Rec = Eigen::MatrixXd::Zero(2*Nr, 2);

                if ( settings.flag_method[0] == "SPOD" ) {
                     Rec = Reconstruction_S_POD(t_vec,
                                           K_pc, lambda_POD, eig_vec.transpose(),
                                           Phi.middleRows(Nr, 2 * Nr), settings.t_rec[nt],
                                           Nm,
                                          "VECTOR-2D",
                                           settings.flag_interp);

                } else if ( settings.flag_method[0] == "DMD" ) {
                    Eigen::MatrixXcd PhiTPhi = Phi_DMD.leftCols(Nm).transpose()*Phi_DMD.leftCols(Nm);
                    Eigen::MatrixXcd Coeffs = PhiTPhi.inverse()*(Phi_DMD.leftCols(Nm).transpose()*sn_set);
                    Rec =  Reconstruction_DMD_Interp ( settings.t_rec[nt],
                                                        t_vec,
                                                       Coeffs,
                                                         Phi_DMD.middleRows(Nr, 2 * Nr),
                                                         "VECTOR-2D",
                                                         settings.flag_interp );
                }
                std::cout << "Done" << std::endl;

                Rec.col(0) = Rec.col(0) * (rhoU_max - rhoU_min) + Eigen::MatrixXd::Ones(Nr, 1)*rhoU_min;
                Rec.col(1) = Rec.col(1) * (rhoV_max - rhoV_min) + Eigen::MatrixXd::Ones(Nr, 1)*rhoV_min;

                if (settings.flag_mean == "IC") {
                    Rec.col(0) = Rec.col(0) + Ic.middleRows(Nr, Nr);
                    Rec.col(1) = Rec.col(1) + Ic.middleRows(2*Nr, Nr);
                }
                if (settings.flag_mean == "YES") {
                    Rec.col(0) = Rec.col(0) + mean.middleRows(Nr, Nr);
                    Rec.col(1) = Rec.col(1) + mean.middleRows(2*Nr, Nr);
                }

                std::cout << "Writing reconstructed field ..." << "\t";
                std::string filename = root_outputfile + "_uni.dat";

                write_Reconstructed_fields(Rec, Coords,
                                           filename,
                                           "VECTOR-2D", nt);

                std::cout << "Done" << std::endl << std::endl;
            }

        } else {
            std::cout << "Nothing else to do!" << std::endl;
        }
 //                 Sn_Cons_time.middleRows(ncons*Nr,Nr) = Phi.leftCols(Nm)*Sig*eig_vec.leftCols(Nm).transpose();
    }


//    Defining common scope for adaptive sampling
    {
        if ( nVar != Nm ) std::cout << "-----------------WARNING-----------------\n "
                                       "comparison with different number of modes" << std::endl << std::endl;

//        Eigen::VectorXi Ipos = Eigen::VectorXi::Zero(nVar);
//        for ( int ipos = 0; ipos < nVar; ipos++ ) Ipos(ipos) = settings.t_pos[ipos];

//        t_pos << 0, 1, 6, 10, 19, 28, 39, 44, 52, 67, 83,  99;

        //Vector of MatrixXd where to store the evolution in time of conservative variables
        Eigen::MatrixXd Sn_Cons_time = Eigen::MatrixXd::Zero(nC*Nr, 3);
        Eigen::MatrixXd Phi = Eigen::MatrixXd::Zero(Nr,nVar);
        Eigen::VectorXd lambda_POD = Eigen::VectorXd::Zero(nVar);
        Eigen::MatrixXcd Phi_DMD = Eigen::MatrixXd::Zero(Nr,settings.Ns);
        Eigen::VectorXcd lambda_DMD = Eigen::VectorXd::Zero(settings.Ns);
        Eigen::VectorXd K_pc(nVar);
        Eigen::MatrixXd eig_vec(nVar, nVar);

        std::vector<smartuq::surrogate::rbf> surr_coefs_POD;
        std::vector<smartuq::surrogate::rbf> surr_coefs_RDMD;
        std::vector<smartuq::surrogate::rbf> surr_coefs_DMD_r;
        std::vector<smartuq::surrogate::rbf> surr_coefs_DMD_i;

        int N_notZero;

        Eigen::MatrixXd sub_sn_set = indexing(sn_set, Eigen::ArrayXi::LinSpaced(nC*Nr,0,nC*Nr-1),Ipos);

        if ( settings.flag_method[0] == "SPOD") {
            std::cout << "Computing adaptive SPOD modes with Nf : " << settings.Nf[0] << "\n";
            Phi = SPOD_basis(sub_sn_set,
                                 lambda_POD, K_pc, eig_vec,
                                 settings.Nf[0],
                                 settings.flag_bc,
                                 settings.flag_filter,
                                 settings.sigma);
            Nm = Phi.cols();
            Eigen::MatrixXd Coeffs = Phi.transpose() * sn_set;
            surr_coefs_POD = getSurrCoefs(t_vec, Coeffs.transpose(), settings.flag_interp);
        } else if ( settings.flag_method[0] == "DMD") {
            Eigen::MatrixXcd eig_vec_DMD;
            std::cout << "Computing adaptive DMD modes" << std::endl;
//            if (settings.r == 0) Nm = Nmod(settings.En, K_pc);
//            else Nm = std::min(settings.r, settings.Ns);

            //Using Constant delta t shifts
//            Phi_DMD = DMD_Adaptive_basis(sn_set,
//                                lambda_DMD,
//                                eig_vec_DMD,
//                                lambda_POD,
//                                eig_vec,
//                                Ipos);

            //Using non-uniform delta t shifts
            Phi_DMD = DMD_basis(sub_sn_set,
                                lambda_DMD,
                                eig_vec_DMD,
                                lambda_POD,
                                eig_vec,
                                -1);

            Nm = Phi_DMD.cols();
            Eigen::MatrixXcd PhiTPhi = Phi_DMD.transpose()*Phi_DMD;
            Eigen::MatrixXcd Coeffs = PhiTPhi.inverse()*(Phi_DMD.transpose()*sn_set);
            surr_coefs_DMD_r = getSurrCoefs(t_vec, Coeffs.real().transpose(), settings.flag_interp);
            surr_coefs_DMD_i = getSurrCoefs(t_vec, Coeffs.imag().transpose(), settings.flag_interp);
        } else if ( settings.flag_method[0] == "RDMD") {

            std::cout << "Computing adaptive RDMD modes \n";
            Eigen::MatrixXd Coeffs = Eigen::MatrixXd::Zero(settings.Ns,settings.Ns);
            Phi = RDMD_Adaptive_basis ( sn_set,Coeffs,K_pc,Ipos);
            surr_coefs_RDMD = getSurrCoefs(t_vec, Coeffs.transpose(), settings.flag_interp);
        }

        if ( decision == "-e" ) {
            for ( int itr = 0; itr < settings.t_res.size(); itr++ ) {
                std::cout << " Computing residuals at time : " << settings.t_res[itr] << std::endl;
                Modify_su2_cfg ( su2_conf, su2_conf_new, settings.Dt_res[0], settings.t_res[itr], V_inf );
                if ( settings.flag_method[0] == "SPOD" ) {
                    Eigen::MatrixXd coef_t(3, Nm);

                    std::vector<double> tr(1);
                    std::vector<double> t_evaluate = {settings.t_res[itr] - 2.0 * settings.Dt_res[0],
                                                      settings.t_res[itr] - settings.Dt_res[0],
                                                      settings.t_res[itr]};

                    for (int j = 0; j < 3; j++) {
                        tr[0] = t_evaluate[j];
                        for (int i = 0; i < Nm; i++)
                            surr_coefs_POD[i].evaluate(tr, coef_t(j, i));

                    }

                    Sn_Cons_time = Phi * coef_t.transpose();
                } else if ( settings.flag_method[0] == "DMD" ) {

                    Eigen::MatrixXcd coef_t(3, Nm);

                    std::vector<double> tr(1);
                    std::vector<double> t_evaluate = {settings.t_res[itr] - 2.0 * settings.Dt_res[0],
                                                      settings.t_res[itr] - settings.Dt_res[0],
                                                      settings.t_res[itr]};

                    double tmp_r, tmp_i;
                    for (int j = 0; j < 3; j++) {
                        tr[0] = t_evaluate[j];
                        for (int i = 0; i < Nm; i++) {
                            surr_coefs_DMD_r[i].evaluate(tr, tmp_r);
                            surr_coefs_DMD_i[i].evaluate(tr, tmp_i);
                            std::complex<double> c(tmp_r,tmp_i);
                            coef_t(j,i) = c;
                        }
                    }

                    Eigen::MatrixXcd Appo = Phi_DMD * coef_t.transpose();
                    Sn_Cons_time = Appo.real();
                } else if ( settings.flag_method[0] == "RDMD" ) {
                    Eigen::MatrixXd coef_t(3, Nm);

                    std::vector<double> tr(1);
                    std::vector<double> t_evaluate = {settings.t_res[itr] - 2.0 * settings.Dt_res[0],
                                                      settings.t_res[itr] - settings.Dt_res[0],
                                                      settings.t_res[itr]};

                    for (int j = 0; j < 3; j++) {
                        tr[0] = t_evaluate[j];
                        for (int i = 0; i < Nm; i++)
                            surr_coefs_RDMD[i].evaluate(tr, coef_t(j, i));

                    }

                    Sn_Cons_time = Phi * coef_t.transpose();
                }

                if ( flag == "YES" ) {
                    Inverse_Normalization(Sn_Cons_time, settings, nC, rho_max, rho_min, rhoU_max, rhoU_min, rhoV_max,
                                          rhoV_min,rhoW_max, rhoW_min, rhoE_max, rhoE_min, tke_min, tke_max, omega_min,
                                          omega_max, nuTilde_min, nuTilde_max);
                }

                if (settings.flag_mean == "IC") {
                    for (int it = 0; it < 3; it++)
                        Sn_Cons_time.col(it) += Ic;
                } else if (settings.flag_mean == "YES") {
                    for (int it = 0; it < 3; it++)
                        Sn_Cons_time.col(it) += mean;
                }

                std::string mv_string = "mv history_rbm_00002.csv history_adaptive_" + std::to_string(itr) + ".csv";
                len_s = mv_string.length();
                char mv_sys_call[len_s + 1];
                strcpy(mv_sys_call, mv_string.c_str());
                std::cout << "Writing time reconstruction " << std::endl;
                // Write_Restart_Cons_Time( Sn_Cons_time, Coords, settings.out_file, t_evaluate.size(), nC, alpha );
                Write_Restart_Cons_Time(Sn_Cons_time, Coords, settings.out_file, 3, nC, settings.alpha, settings.beta, binary);
                //Executing SU2, removing all useless files, renaming files with residuals
                std::cout << "Calling SU2 for residual evaluation and writing file to history " << std::endl;
                auto otp = std::system(su2_sys_call);
                otp = std::system(rmf_sys_call);
                otp = std::system(mv_sys_call);
            }
        }
        else if ( decision == "-r" ) {
            for ( int nt = 0; nt < settings.t_rec.size(); nt++ ) {

                std::cout << "Reconstructing Momentum at time : " << settings.t_rec[nt] << "\t";
                Eigen::MatrixXd Rec = Eigen::MatrixXd::Zero(2*Nr, 2);

                if ( settings.flag_method[0] == "SPOD" ) {
                    Eigen::MatrixXd Coeffs = Phi.transpose() * sn_set;
                    for (int i = 0; i < Coeffs.rows(); i++) Coeffs.row(i) = Coeffs.row(i) / std::sqrt(lambda_POD(i));

                    Eigen::MatrixXd Rec = Reconstruction_S_POD(t_vec,
                                                               K_pc, lambda_POD, Coeffs,
                                                               Phi.middleRows(Nr, 2 * Nr), settings.t_rec[nt],
                                                               Nm,
                                                               "VECTOR-2D",
                                                               settings.flag_interp);
                } else if ( settings.flag_method[0] == "DMD" ) {
                    Eigen::MatrixXcd PhiTPhi = Phi_DMD.leftCols(Nm).transpose()*Phi_DMD.leftCols(Nm);
                    Eigen::MatrixXcd Coeffs = PhiTPhi.inverse()*(Phi_DMD.leftCols(Nm).transpose()*sn_set);
                    Rec =  Reconstruction_DMD_Interp ( settings.t_rec[nt],
                                                       t_vec,
                                                       Coeffs,
                                                       Phi_DMD.middleRows(Nr, 2 * Nr),
                                                       "VECTOR-2D",
                                                       settings.flag_interp );
                }

                std::cout << "Done" << std::endl;

                Rec.col(0) = Rec.col(0) * (rhoU_max - rhoU_min) + Eigen::MatrixXd::Ones(Nr,1) * rhoU_min;
                Rec.col(1) = Rec.col(1) * (rhoV_max - rhoV_min) + Eigen::MatrixXd::Ones(Nr,1) * rhoV_min;

                if (settings.flag_mean == "IC") {
                    Rec.col(0) = Rec.col(0) + Ic.middleRows(Nr, Nr);
                    Rec.col(1) = Rec.col(1) + Ic.middleRows(2*Nr, Nr);
                }
                if (settings.flag_mean == "YES") {
                    Rec.col(0) = Rec.col(0) + mean.middleRows(Nr, Nr);
                    Rec.col(1) = Rec.col(1) + mean.middleRows(2*Nr, Nr);
                }

                std::cout << "Writing reconstructed field ..." << "\t";
                std::string filename = root_outputfile + "_adapt.dat";

                write_Reconstructed_fields(Rec, Coords,
                                           filename,
                                           "VECTOR-2D", nt);

                std::cout << "Done" << std::endl << std::endl;
            }
        }
        else {
            std::cout << "Nothing to do!" << std::endl;
        }
    }

    std::cout << "MODES Adaptive Residual Evaluation ends" << std::endl << std::endl;
    return 0;

}



