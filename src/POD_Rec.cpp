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
    std::cout << "POD Reconstruction starts" << std::endl;

    prob_settings settings;
    std::string filecfg = argv[1];

    //Reading configuration file
    Read_cfg( filecfg, settings );
    if ( settings.flag_prob != "CONSERVATIVE"){
        std::cout << "Reconstruction with residual evaluation only implemented for Conservative Variables flag \n "
                     "FLAG_PROB must be CONSERVATIVE\n Exiting ..." << std::endl;
        exit(EXIT_FAILURE);
    }
    double t_0 = settings.nstart*settings.Dt_cfd;
    double alpha = settings.alpha;
    double beta  = settings.beta;

    int Nf_SPOD = 0;

    std::string root_inputfile;
    root_inputfile.assign ( settings.in_file, 0, settings.in_file.size() - 4);
    std::string input_format;
    input_format.assign ( settings.in_file, settings.in_file.size() - 3, 3);

    std::stringstream buffer;
    buffer << std::setfill('0') << std::setw(5) << std::to_string(settings.nstart);
    std::string file_1 = root_inputfile + "_" + buffer.str() + "." + input_format;

    int nC = settings.Cols.size();
    int Nr = N_gridpoints ( file_1 );
    std::cout << "Number of grid points : " << Nr << std::endl;

    std::cout << "Reading Coordinates ... \t ";
    Eigen::MatrixXd Coords = read_col( file_1, Nr, settings.Cols_coords );
    std::cout << "Done " << std::endl;

    std::cout << "Storing snapshot Matrix ... \n ";
    Eigen::MatrixXd sn_set = generate_snap_matrix( Nr, settings);
    // Eigen::MatrixXd sn_set = Eigen::MatrixXd::Zero(settings.ndim*Nr, settings.Ns);
    std::cout << "Computing mean/Initial Condition of CFD solution ... " << std::endl;
    //Defining Initial condition

    Eigen::VectorXd Ic = Eigen::VectorXd::Zero(nC*Nr);
    Eigen::VectorXd mean = sn_set.rowwise().mean();

    if ( settings.flag_mean == "YES" ) {
        for ( int it = 0; it < settings.Ns; it++ )
            sn_set.col(it) -= mean;
    }

    if ( settings.flag_mean == "IC" ) {
        Ic = IC(sn_set, settings, nC, Nr);
    }
//POD basis computation

    std::vector<Eigen::VectorXd> lambda(settings.ndim+2);
    std::vector<Eigen::VectorXd> K_pc(settings.ndim+2);
    std::vector<Eigen::MatrixXd> eig_vec(settings.ndim+2);
    std::vector<Eigen::MatrixXd> Phi(settings.ndim+2);
    std::vector<int> Nm(settings.ndim+2);

    for ( int iDim = 0; iDim < settings.ndim+2; iDim++){

        lambda[iDim] = Eigen::VectorXd::Zero(settings.Ns);
        K_pc[iDim] = Eigen::VectorXd::Zero(settings.Ns);
        eig_vec[iDim] = Eigen::MatrixXd::Zero(settings.Ns,settings.Ns);
        Phi[iDim] = SPOD_basis( sn_set.middleRows(iDim*Nr,Nr),
                                          lambda[iDim], K_pc[iDim], eig_vec[iDim],
                                          Nf_SPOD,
                                          settings.flag_bc,
                                          settings.flag_filter,
                                          settings.sigma);

        if ( settings.r == 0 ) {
            Nm[iDim] = Nmod(settings.En, K_pc[iDim]);
            std::cout << "Number of modes for desired energetic content: " << Nm[iDim] << std::endl;
        } else {
            int Nmax = Phi[iDim].cols();
            Nm[iDim] = std::min(settings.r,Nmax);
            std::cout << "Number of modes (fixed): " << Nm[iDim] << std::endl;
        }
    }

    std::vector<double> t_v( settings.Ns );
    t_v[0] = (double)settings.nstart*settings.Dt_cfd;

    for ( int kt = 1; kt < settings.Ns; kt++ )
        t_v[kt] = t_v[kt-1] + settings.Dt_cfd*(double)settings.Ds;


//POD reconstruction on each selected time step

    for ( int i = 0; i < settings.t_rec.size(); i++ )
    {

        Eigen::VectorXd Rec_rho(Nr);
        Eigen::VectorXd Rec_rhoU(Nr);
        Eigen::VectorXd Rec_rhoV(Nr);
        Eigen::VectorXd Rec_rhoW(Nr);
        Eigen::VectorXd Rec_rhoE(Nr);

        std::vector<int> pos = {};
        std::cout << " POD reconstruction at time : " << settings.t_rec[i] << std::endl;

        for ( int iDim = 0; iDim < settings.ndim + 2; iDim ++ ){

            Eigen::MatrixXd Rec = Reconstruction_S_POD ( t_v,
                                                         K_pc[iDim], lambda[iDim], eig_vec[iDim].transpose(),
                                                         Phi[iDim], settings.t_rec[i],
                                                         Nm[iDim],
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

        Eigen::MatrixXd Rec_M(Nr, settings.ndim + 2);

        if ( settings.ndim == 2)
        {
            Rec_M.col(0) = Rec_rho;
            Rec_M.col(1) = Rec_rhoU;
            Rec_M.col(2) = Rec_rhoV;
            Rec_M.col(3) = Rec_rhoE;
        } else
        {
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

    std::cout << "POD Reconstruction ends " << std::endl;

    return 0;
}