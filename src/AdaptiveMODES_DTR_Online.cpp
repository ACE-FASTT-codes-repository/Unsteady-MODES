/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------ Copyright (C) 2021 University of Strathclyde and Authors ------
-------------------- e-mail: gaetano.pascarella.ac.uk ----------------
----------------------- Author: Gaetano Pascarella -----------------------

Code for adaptive reconstruction based on residual evaluation
Input config file + error file (+ Modes,Coefs and Encontent RDMD if already available)

Output reconstructed field at the desired time instants with the adaptive technique
based on residual evaluation
*/

#include "Generate_snset.hpp"
#include "Reconstruction.hpp"
#include "write_Outputs.hpp"
#include "Pre-Process.hpp"

int main( int argc, char *argv[] )
{
    std::cout << "Adaptive Reconstruction with ResEval starts " << std::endl;

    prob_settings settings;
    std::string filecfg = argv[1];

    //Reading configuration file
    Read_cfg( filecfg, settings );
    if ( settings.flag_prob != "CONSERVATIVE"){
        std::cout << "Reconstruction with residual evaluation only implemented for COnservative Variables flag \n "
                     "FLAG_PROB must be CONSERVATIVE\n Exiting ..." << std::endl;
        exit(EXIT_FAILURE);
    }
    double t_0 = settings.nstart*settings.Dt_cfd;
    double alpha = settings.alpha;
    double beta = settings.beta;
    if (settings.ndim == 2) beta = 0.0;

    int s_Nf = 1;
    int Nmethods = s_Nf + 2;
    std::vector<int> Nf(s_Nf);
    Nf[0] = 0;
//    Nf[1] = std::ceil(settings.Ns/10.0);
//    Nf[2] = std::ceil(settings.Ns/2.0);
//    Nf[3] = std::ceil(2.0*settings.Ns/3.0);
//    Nf[4] = settings.Ns;
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

    // Eigen::MatrixXd sn_set = Eigen::MatrixXd::Zero(settings.ndim*Nr, settings.Ns);
    std::cout << "Computing mean/Initial Condition of CFD solution ... " << std::endl;
    //Defining Mean/Initial condition
    Eigen::VectorXd mean = sn_set.rowwise().mean();
    int nC = settings.Cols.size();
    Eigen::VectorXd Ic = Eigen::VectorXd::Zero(nC*Nr);

    if ( settings.flag_mean == "YES" ) {
        for ( int it = 0; it < settings.Ns; it++ )
            sn_set.col(it) -= mean;
    }

    if ( settings.flag_mean == "IC" ) {
        Ic = IC(sn_set,settings,nC,Nr);
    }


    std::cout << "Storing database for all methods" << std::endl;
    //--------------------------------------------------------//
    //POD
    //--------------------------------------------------------//
    std::cout << "POD..." << std::endl;
    std::vector<Eigen::VectorXd> lambdaPOD(settings.ndim + 2);
    std::vector<Eigen::VectorXd> K_pcPOD(settings.ndim + 2);
    std::vector<Eigen::MatrixXd> eig_vecPOD(settings.ndim + 2);
    std::vector<Eigen::MatrixXd> PhiPOD(settings.ndim + 2);
    Eigen::VectorXi NmPOD(settings.ndim + 2);

    for (int iDim = 0; iDim < settings.ndim + 2; iDim++ ){
        lambdaPOD[iDim] = Eigen::VectorXd::Zero(settings.Ns);
        K_pcPOD[iDim] = Eigen::VectorXd::Zero(settings.Ns);
        eig_vecPOD[iDim] = Eigen::MatrixXd::Zero(settings.Ns,settings.Ns);
        PhiPOD[iDim] = SPOD_basis( sn_set.middleRows(iDim*Nr,Nr),
                                   lambdaPOD[iDim], K_pcPOD[iDim], eig_vecPOD[iDim],
                                   Nf_SPOD,
                                   settings.flag_bc,
                                   settings.flag_filter,
                                   settings.sigma);
        if ( settings.r == 0 ) {
            NmPOD[iDim] = Nmod(settings.En, K_pcPOD[iDim]);
            std::cout << "Number of modes for desired energetic content: " << NmPOD[iDim] << std::endl;
        } else {
            NmPOD[iDim] = std::min(settings.r,settings.Ns);
            std::cout << "Number of modes (fixed): " << NmPOD[iDim] << std::endl;
        }
    }

    //Vector of POD times for interpolation
    std::vector<double> t_v( settings.Ns );
    t_v[0] = (double)settings.nstart*settings.Dt_cfd;
    for ( int kt = 1; kt < settings.Ns; kt++ )
        t_v[kt] = t_v[kt-1] + settings.Dt_cfd*(double)settings.Ds;

    //--------------------------------------------------------//
    //DMD
    //--------------------------------------------------------//
    std::cout << "DMD..." << std::endl;
    Eigen::VectorXd lambda_POD;
    Eigen::MatrixXd eig_vec_POD;
    Eigen::MatrixXcd eig_vec_DMD;
    std::vector<Eigen::VectorXcd> lambda_DMD(settings.ndim + 2);
    std::vector<Eigen::MatrixXcd> PhiDMD(settings.ndim + 2);
    std::vector<Eigen::VectorXcd> alfaDMD(settings.ndim + 2);

    for (int iDim = 0; iDim < settings.ndim + 2; iDim++){

        if ( settings.r == 0 ) PhiDMD[iDim] = DMD_basis( sn_set.middleRows(iDim*Nr,Nr), lambda_DMD[iDim], eig_vec_DMD, lambda_POD, eig_vec_POD, -1 );
        else PhiDMD[iDim] = DMD_basis( sn_set.middleRows(iDim*Nr,Nr), lambda_DMD[iDim], eig_vec_DMD, lambda_POD, eig_vec_POD, settings.r );

        alfaDMD[iDim] = Calculate_Coefs_DMD_exact ( sn_set.middleRows(iDim*Nr,Nr).leftCols(settings.Ns-1), lambda_DMD[iDim], PhiDMD[iDim] );

        // std::cout << "Reordering modes DMD ... " << "\t";
        Eigen::VectorXd En = Eigen::VectorXd::Zero(PhiDMD[iDim].cols());
        double T = t_v[t_v.size()-1];

        Eigen::VectorXcd omega(PhiDMD[iDim].cols());
        for ( int idmd = 0; idmd < PhiDMD[iDim].cols(); idmd++ )
            omega(idmd) = std::log(lambda_DMD[iDim](idmd))/(settings.Dt_cfd*(double)settings.Ds);

        for ( int idmd = 0 ; idmd < PhiDMD[iDim].cols(); idmd ++ ) {
            double alfa_i = alfaDMD[iDim](idmd).imag();
            double alfa_r = alfaDMD[iDim](idmd).real();
            double sigma = omega(idmd).real();
            En(idmd) = (alfa_r*alfa_r + alfa_i*alfa_i)*(std::exp(2.0*sigma*T) - 1.0)/(2.0*sigma);
        }
        dmd_sort( En, PhiDMD[iDim], lambda_DMD[iDim], alfaDMD[iDim]);
    }

    //--------------------------------------------------------//
    //RDMD
    //--------------------------------------------------------//
    std::cout << "RDMD..." << std::endl;
    Eigen::VectorXd lambda = Eigen::VectorXd::Zero(settings.Ns);
    Eigen::VectorXd K_pc = Eigen::VectorXd::Zero(settings.Ns);
    std::vector<Eigen::MatrixXd> PhiRDMD(settings.ndim + 2);
    std::vector<Eigen::MatrixXd> CoefsRDMD(settings.ndim + 2);
    Eigen::VectorXi NmRDMD(settings.ndim + 2);

    if ( settings.flag_wdb_be == "READ") {

        for (int iDim = 0; iDim < settings.ndim + 2; iDim++) {
            PhiRDMD[iDim] = Eigen::MatrixXd::Zero(Nr,settings.r_RDMD);
            CoefsRDMD[iDim] = Eigen::MatrixXd::Zero(settings.Ns, settings.r_RDMD);
            NmRDMD[iDim] = settings.r_RDMD;
            std::cout << "Number of Modes to use in Reconstruction: " << NmRDMD[iDim] << std::endl;
        }
        ModeDB_Read("ModesRDMD_", "CoefsRDMD_", PhiRDMD, CoefsRDMD, settings);
    } else {

        for (int iDim = 0; iDim < settings.ndim + 2; iDim++ ){
            PhiRDMD[iDim] = Eigen::MatrixXd::Zero(Nr,settings.r_RDMD);
            CoefsRDMD[iDim] = Eigen::MatrixXd::Zero(settings.Ns, settings.r_RDMD);
            NmRDMD[iDim] = settings.r_RDMD;
            std::cout << "Number of Modes to use in Reconstruction: " << NmRDMD[iDim] << std::endl;
            PhiRDMD[iDim] = RDMD_modes_coefs(sn_set.middleRows(iDim * Nr, Nr), CoefsRDMD[iDim], lambda, K_pc, -1, NmRDMD[iDim], settings.En);
        }
    }

    //---------------------------------------------------------//
    //---------------------------------------------------------//

    std::cout << "Reading Residuals ... " << std::endl;

    std::vector<std::string> resfilename = {"history_pod.csv", "history_dmd.csv", "history_rdmd.csv"};

    Eigen::MatrixXd Err_RBM = Eigen::MatrixXd::Zero(settings.t_res.size(), Nmethods);
    Eigen::MatrixXd Err_RBM_rho = Eigen::MatrixXd::Zero(settings.t_res.size(), Nmethods);
    Eigen::MatrixXd Err_RBM_rhoV = Eigen::MatrixXd::Zero(settings.t_res.size(), Nmethods);
    Eigen::MatrixXd Err_RBM_rhoU = Eigen::MatrixXd::Zero(settings.t_res.size(), Nmethods);
    Eigen::MatrixXd Err_RBM_rhoW = Eigen::MatrixXd::Zero(settings.t_res.size(), Nmethods);
    Eigen::MatrixXd Err_RBM_rhoE = Eigen::MatrixXd::Zero(settings.t_res.size(), Nmethods);

    for ( int i = 0; i < resfilename.size(); i++ ) {
        std::ifstream file_data;
        file_data.open( resfilename[i] );
        if ( !file_data.is_open() ) {
            std::cout << "File : " << resfilename[i] << " not found" << std::endl;
            exit (EXIT_FAILURE);
        }

        std::string line_flow_data ;
        int n_row = 0, count = 0;
        // Reading row of headers
        getline( file_data, line_flow_data );

        while ( getline( file_data, line_flow_data ) && n_row <  settings.t_res.size()) {
            std::istringstream iss(line_flow_data);
            std::string token;
            double err;
            int idx;
            count = 0;
            while( getline( iss, token, ',') ) {
                err = std::stod(token);
                if (settings.ndim == 3) {
                    if (count == 17) Err_RBM_rho(n_row, i) = std::pow(10.0, err);
                    if (count == 18) Err_RBM_rhoU(n_row, i) = std::pow(10.0, err);
                    if (count == 19) Err_RBM_rhoV(n_row, i) = std::pow(10.0, err);
                    if (count == 20) Err_RBM_rhoW(n_row, i) = std::pow(10.0, err);
                    if (count == 21) Err_RBM_rhoE(n_row, i) = std::pow(10.0, err);
                }

                if (settings.ndim == 2) {
                    if (count == 17) Err_RBM_rho(n_row, i) = std::pow(10.0, err);
                    if (count == 18) Err_RBM_rhoU(n_row, i) = std::pow(10.0, err);
                    if (count == 19) Err_RBM_rhoV(n_row, i) = std::pow(10.0, err);
                    if (count == 20) Err_RBM_rhoE(n_row, i) = std::pow(10.0, err);
                }
                count ++;
            }
            n_row++;
        }
        file_data.close();
    }

//Adaptive reconstruction on each selected time step
    int best_method_idx;

    std::cout << "Initializing Vector of time ... " << std::endl;
    Eigen::VectorXd t_vec( settings.t_res.size());
    for ( int it = 0; it < settings.t_res.size(); it++ ) t_vec(it) = settings.t_res[it];

    double tol = settings.tol;
    int index1, index2;
    Eigen::VectorXd Err_interp(Nmethods);

    std::vector<double> t_pod = {};
    std::vector<double> t_dmd = {};
    std::vector<double> t_rdmd = {};

    for ( int i = 0; i < settings.t_rec.size(); i++ ) {

        Eigen::VectorXd Rec_rho(Nr);
        Eigen::VectorXd Rec_rhoU(Nr);
        Eigen::VectorXd Rec_rhoV(Nr);
        Eigen::VectorXd Rec_rhoW(Nr);
        Eigen::VectorXd Rec_rhoE(Nr);

        std::vector<int> pos = {};
        std::cout << " Adaptive reconstruction at time : " << settings.t_rec[i] << std::endl;

        index1 = 0;
        index2 = 0;
        for ( int nt = 0; nt < t_vec.size()-1; nt ++ ) {
            if ( (settings.t_rec[i] >= t_vec(nt)) && (settings.t_rec[i] <= t_vec(nt+1)) ) {
                index1 = nt;
                index2 = nt+1;
                break;
            }
        }

        if ( index1 == index2 ) {
            std::cout << "Time for reconstruction out of interval!" << std::endl;
            continue;
        }

        for ( int iDim = 0; iDim < settings.ndim + 2; iDim ++ ) {

            if (iDim == 0) Err_RBM = Err_RBM_rho;
            if (iDim == 1) Err_RBM = Err_RBM_rhoU;
            if (iDim == 2) Err_RBM = Err_RBM_rhoV;
            if (iDim == 3 && settings.ndim == 2) Err_RBM = Err_RBM_rhoE;
            if (iDim == 3 && settings.ndim == 3) Err_RBM = Err_RBM_rhoW;
            if (iDim == 4 ) Err_RBM = Err_RBM_rhoE;

            int count = 0;
            double Dt = t_vec[index2] - t_vec[index1];
            for (int k = 0; k < Nmethods; k++) {
                Err_interp(k) = Err_RBM(index1, k) + (Err_RBM(index2, k) - Err_RBM(index1, k)) /
                                                     Dt * (settings.t_rec[i] - t_vec[index1]);
            }
            std::cout << std::endl;
            double eps = Err_interp.minCoeff(&best_method_idx);

            //FIX THIS FUNCTION
            std::string method = method_selected(best_method_idx, Nf_SPOD, Nf);
            std::cout << "Best method is " << method << std::endl;
            std::cout << " Error : " << Err_interp(best_method_idx) << std::endl;

            std::cout << "Computing Reconstruction using selected methods " << std::endl;

            if ( method == "SPOD" ) {
                Eigen::MatrixXd Rec = Reconstruction_S_POD ( t_v,
                                                             K_pcPOD[iDim], lambdaPOD[iDim], eig_vecPOD[iDim].transpose(),
                                                             PhiPOD[iDim], settings.t_rec[i],
                                                             NmPOD[iDim],
                                                             "SCALAR",
                                                             settings.flag_interp ) ;

                if ( iDim == 0 && settings.flag_mean == "NO" ) Rec_rho = Rec.col(0);
                if ( iDim == 1 && settings.flag_mean == "NO" ) Rec_rhoU = Rec.col(0);
                if ( iDim == 2 && settings.flag_mean == "NO" ) Rec_rhoV = Rec.col(0);
                if ( iDim == 3 && settings.ndim == 2 && settings.flag_mean == "NO" ) Rec_rhoE = Rec.col(0);
                if ( iDim == 3 && settings.ndim == 3 && settings.flag_mean == "NO" ) Rec_rhoW = Rec.col(0);
                if ( iDim == 4 && settings.flag_mean == "NO" ) Rec_rhoE = Rec.col(0);

                if ( iDim == 0 && settings.flag_mean == "IC" ) Rec_rho = Rec.col(0) + Ic.middleRows(0,Nr);
                if ( iDim == 1 && settings.flag_mean == "IC" ) Rec_rhoU = Rec.col(0) + Ic.middleRows(Nr,Nr);
                if ( iDim == 2 && settings.flag_mean == "IC" ) Rec_rhoV = Rec.col(0) + Ic.middleRows(2*Nr,Nr);
                if ( iDim == 3 && settings.ndim == 2 && settings.flag_mean == "IC" ) Rec_rhoE = Rec.col(0) + Ic.middleRows(3*Nr,Nr);
                if ( iDim == 3 && settings.ndim == 3 && settings.flag_mean == "IC" ) Rec_rhoW = Rec.col(0) + Ic.middleRows(3*Nr,Nr);
                if ( iDim == 4 && settings.flag_mean == "IC" ) Rec_rhoE = Rec.col(0) + Ic.middleRows(4*Nr,Nr);

                if ( iDim == 0 && settings.flag_mean == "YES" ) Rec_rho = Rec.col(0) + mean.middleRows(0,Nr);
                if ( iDim == 1 && settings.flag_mean == "YES" ) Rec_rhoU = Rec.col(0) + mean.middleRows(Nr,Nr);
                if ( iDim == 2 && settings.flag_mean == "YES" ) Rec_rhoV = Rec.col(0) + mean.middleRows(2*Nr,Nr);
                if ( iDim == 3 && settings.ndim == 2 && settings.flag_mean == "YES" ) Rec_rhoE = Rec.col(0) + mean.middleRows(3*Nr,Nr);
                if ( iDim == 3 && settings.ndim == 3 && settings.flag_mean == "YES" ) Rec_rhoW = Rec.col(0) + mean.middleRows(3*Nr,Nr);
                if ( iDim == 4 && settings.flag_mean == "YES" ) Rec_rhoE = Rec.col(0) + mean.middleRows(4*Nr,Nr);
            }

            if ( method == "DMD" ) {

                double t_dmd = settings.t_rec[i] - (double)settings.nstart*settings.Dt_cfd;//t_vec(0);
                Eigen::MatrixXcd Rec = Reconstruction_DMD ( t_dmd,
                                                            settings.Dt_cfd*settings.Ds,
                                                            alfaDMD[iDim],
                                                            PhiDMD[iDim],
                                                            lambda_DMD[iDim],
                                                            "SCALAR" );

                if ( iDim == 0 && settings.flag_mean == "NO" ) Rec_rho = Rec.real().col(0);
                if ( iDim == 1 && settings.flag_mean == "NO" ) Rec_rhoU = Rec.real().col(0);
                if ( iDim == 2 && settings.flag_mean == "NO" ) Rec_rhoV = Rec.real().col(0);
                if ( iDim == 3 && settings.ndim == 2 && settings.flag_mean == "NO" ) Rec_rhoE = Rec.real().col(0);
                if ( iDim == 3 && settings.ndim == 3 && settings.flag_mean == "NO" ) Rec_rhoW = Rec.real().col(0);
                if ( iDim == 4 && settings.flag_mean == "NO" ) Rec_rhoE = Rec.real().col(0);

                if ( iDim == 0 && settings.flag_mean == "IC" ) Rec_rho = Rec.real().col(0) + Ic.middleRows(0,Nr);
                if ( iDim == 1 && settings.flag_mean == "IC" ) Rec_rhoU = Rec.real().col(0) + Ic.middleRows(Nr,Nr);
                if ( iDim == 2 && settings.flag_mean == "IC" ) Rec_rhoV = Rec.real().col(0) + Ic.middleRows(2*Nr,Nr);
                if ( iDim == 3 && settings.ndim == 2 && settings.flag_mean == "IC" ) Rec_rhoE = Rec.real().col(0) + Ic.middleRows(3*Nr,Nr);
                if ( iDim == 3 && settings.ndim == 3 && settings.flag_mean == "IC" ) Rec_rhoW = Rec.real().col(0) + Ic.middleRows(3*Nr,Nr);
                if ( iDim == 4 && settings.flag_mean == "IC" ) Rec_rhoE = Rec.real().col(0) + Ic.middleRows(4*Nr,Nr);

                if ( iDim == 0 && settings.flag_mean == "YES" ) Rec_rho = Rec.real().col(0) + mean.middleRows(0,Nr);
                if ( iDim == 1 && settings.flag_mean == "YES" ) Rec_rhoU = Rec.real().col(0) + mean.middleRows(Nr,Nr);
                if ( iDim == 2 && settings.flag_mean == "YES" ) Rec_rhoV = Rec.real().col(0) + mean.middleRows(2*Nr,Nr);
                if ( iDim == 3 && settings.ndim == 2 && settings.flag_mean == "YES" ) Rec_rhoE = Rec.real().col(0) + mean.middleRows(3*Nr,Nr);
                if ( iDim == 3 && settings.ndim == 3 && settings.flag_mean == "YES" ) Rec_rhoW = Rec.real().col(0) + mean.middleRows(3*Nr,Nr);
                if ( iDim == 4 && settings.flag_mean == "YES" ) Rec_rhoE = Rec.real().col(0) + mean.middleRows(4*Nr,Nr);

            }

            if ( method == "RDMD" ) {

                Eigen::MatrixXd Coefs;
                if (settings.flag_wdb_be == "READ") Coefs = CoefsRDMD[iDim].transpose();
                else Coefs = CoefsRDMD[iDim];

                Eigen::MatrixXd Rec = Reconstruction_RDMD ( settings.t_rec[i],
                                                            t_v,
                                                            Coefs.topRows(NmRDMD[iDim]),
                                                            PhiRDMD[iDim].leftCols(NmRDMD[iDim]),
                                                            "SCALAR",
                                                            settings.flag_interp );

                if ( iDim == 0 && settings.flag_mean == "NO" ) Rec_rho = Rec.col(0);
                if ( iDim == 1 && settings.flag_mean == "NO" ) Rec_rhoU = Rec.col(0);
                if ( iDim == 2 && settings.flag_mean == "NO" ) Rec_rhoV = Rec.col(0);
                if ( iDim == 3 && settings.ndim == 2 && settings.flag_mean == "NO" ) Rec_rhoE = Rec.col(0);
                if ( iDim == 3 && settings.ndim == 3 && settings.flag_mean == "NO" ) Rec_rhoW = Rec.col(0);
                if ( iDim == 4 && settings.flag_mean == "NO" ) Rec_rhoE = Rec.col(0);

                if ( iDim == 0 && settings.flag_mean == "IC" ) Rec_rho = Rec.col(0) + Ic.middleRows(0,Nr);
                if ( iDim == 1 && settings.flag_mean == "IC" ) Rec_rhoU = Rec.col(0) + Ic.middleRows(Nr,Nr);
                if ( iDim == 2 && settings.flag_mean == "IC" ) Rec_rhoV = Rec.col(0) + Ic.middleRows(2*Nr,Nr);
                if ( iDim == 3 && settings.ndim == 2 && settings.flag_mean == "IC" ) Rec_rhoE = Rec.col(0) + Ic.middleRows(3*Nr,Nr);
                if ( iDim == 3 && settings.ndim == 3 && settings.flag_mean == "IC" ) Rec_rhoW = Rec.col(0) + Ic.middleRows(3*Nr,Nr);
                if ( iDim == 4 && settings.flag_mean == "IC" ) Rec_rhoE = Rec.col(0) + Ic.middleRows(4*Nr,Nr);

                if ( iDim == 0 && settings.flag_mean == "YES" ) Rec_rho = Rec.col(0) + mean.middleRows(0,Nr);
                if ( iDim == 1 && settings.flag_mean == "YES" ) Rec_rhoU = Rec.col(0) + mean.middleRows(Nr,Nr);
                if ( iDim == 2 && settings.flag_mean == "YES" ) Rec_rhoV = Rec.col(0) + mean.middleRows(2*Nr,Nr);
                if ( iDim == 3 && settings.ndim == 2 && settings.flag_mean == "YES" ) Rec_rhoE = Rec.col(0) + mean.middleRows(3*Nr,Nr);
                if ( iDim == 3 && settings.ndim == 3 && settings.flag_mean == "YES" ) Rec_rhoW = Rec.col(0) + mean.middleRows(3*Nr,Nr);
                if ( iDim == 4 && settings.flag_mean == "YES" ) Rec_rhoE = Rec.col(0) + mean.middleRows(4*Nr,Nr);

            }

        }

        Eigen::MatrixXd Rec_M(Nr, settings.ndim + 2);

        if ( settings.ndim == 2) {
            Rec_M.col(0) = Rec_rho;
            Rec_M.col(1) = Rec_rhoU;
            Rec_M.col(2) = Rec_rhoV;
            Rec_M.col(3) = Rec_rhoE;
        } else {
            Rec_M.col(0) = Rec_rho;
            Rec_M.col(1) = Rec_rhoU;
            Rec_M.col(2) = Rec_rhoV;
            Rec_M.col(3) = Rec_rhoW;
            Rec_M.col(4) = Rec_rhoE;
        }
        std::cout << "Writing reconstructed field ..." << "\t";
        write_Reconstructed_fields ( Rec_M, Coords, settings.out_file, "CONSERVATIVE", i );
        std::cout << "Done" << std::endl << std::endl << std::endl;
    }

    std::cout << "Adaptive Reconstruction MODES ends " << std::endl;

    return 0;
}
