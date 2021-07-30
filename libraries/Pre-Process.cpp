#include "Pre-Process.hpp"

void common_vars( int &Nr, Eigen::MatrixXd &Coords, std::vector<double> &t_vec, prob_settings settings){

    std::cout << "Initializing Vector of times ... " << std::endl;

    t_vec[0] = (double)settings.nstart*settings.Dt_cfd;
    for ( int i = 1; i < settings.Ns; i++ )
        t_vec[i] = t_vec[i-1] + settings.Dt_cfd*settings.Ds;

    std::cout << t_vec[0] << "-------time of first snapshot" << std::endl;
    std::cout << t_vec[t_vec.size()-1] <<"-------time of last snapshot" << std::endl;
    std::cout << settings.t_res[0] <<"----time of first residual evaluation"<< std::endl;
    std::cout << settings.t_res[settings.t_res.size()-1]<< "----time of last residual evaluation"<< std::endl;

    for ( int idtr = 0; idtr < settings.Dt_res.size(); idtr++ ){
        if ( ((settings.t_res[0] - (2.0 * settings.Dt_res[idtr])) < t_vec[0]) ||  (settings.t_res[settings.t_res.size()-1] > t_vec[t_vec.size() - 1]))
        {
            std::cout
                    << "Define proper Delta_t_res and T_RES vector " << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    std::string root_inputfile;
    root_inputfile.assign ( settings.in_file, 0, settings.in_file.size() - 4);
    std::string input_format;
    input_format.assign ( settings.in_file, settings.in_file.size() - 3, 3);

    std::stringstream buffer;
    buffer << std::setfill('0') << std::setw(5) << std::to_string(settings.nstart);
    std::string file_1 = root_inputfile + "_" + buffer.str() + "." + input_format;

    // Calculate number of grid points
    Nr = N_gridpoints ( file_1 );
    std::cout << "Number of grid points : " << Nr << std::endl;

    std::cout << "Reading Coordinates ... \t ";
    Coords = read_col( file_1, Nr, settings.Cols_coords );
    std::cout << "Done " << std::endl;

}

void get_MinMax_ConsVar (const Eigen::MatrixXd sn_set, const prob_settings &settings, const int nC, double &rho_max,
                         double &rho_min, double &rhoU_max, double &rhoU_min, double &rhoV_max, double &rhoV_min,
                         double &rhoW_max, double &rhoW_min, double &rhoE_max, double &rhoE_min, double &tke_min,
                         double &tke_max, double &omega_min, double &omega_max, double &nuTilde_min, double &nuTilde_max) {

    int Nr = sn_set.rows()/nC;
    if (settings.ndim == 2 && nC == 4) {

        rho_max = sn_set.middleRows(0, Nr).maxCoeff();
        rho_min = sn_set.middleRows(0, Nr).minCoeff();
        rhoU_max = sn_set.middleRows(Nr, Nr).maxCoeff();
        rhoU_min = sn_set.middleRows(Nr, Nr).minCoeff();
        rhoV_max = sn_set.middleRows(2*Nr, Nr).maxCoeff();
        rhoV_min = sn_set.middleRows(2*Nr, Nr).minCoeff();
        rhoE_max = sn_set.middleRows(3*Nr, Nr).maxCoeff();
        rhoE_min = sn_set.middleRows(3*Nr, Nr).minCoeff();

    }
    else if (settings.ndim == 2 && nC == 5) {

        rho_max = sn_set.middleRows(0, Nr).maxCoeff();
        rho_min = sn_set.middleRows(0, Nr).minCoeff();
        rhoU_max = sn_set.middleRows(Nr, Nr).maxCoeff();
        rhoU_min = sn_set.middleRows(Nr, Nr).minCoeff();
        rhoV_max = sn_set.middleRows(2*Nr, Nr).maxCoeff();
        rhoV_min = sn_set.middleRows(2*Nr, Nr).minCoeff();
        rhoE_max = sn_set.middleRows(3*Nr, Nr).maxCoeff();
        rhoE_min = sn_set.middleRows(3*Nr, Nr).minCoeff();
        nuTilde_max = sn_set.middleRows(4*Nr, Nr).maxCoeff();
        nuTilde_min = sn_set.middleRows(4*Nr, Nr).minCoeff();

    }
    else if (settings.ndim == 2 && nC == 6) {

        rho_max = sn_set.middleRows(0, Nr).maxCoeff();
        rho_min = sn_set.middleRows(0, Nr).minCoeff();
        rhoU_max = sn_set.middleRows(Nr, Nr).maxCoeff();
        rhoU_min = sn_set.middleRows(Nr, Nr).minCoeff();
        rhoV_max = sn_set.middleRows(2*Nr, Nr).maxCoeff();
        rhoV_min = sn_set.middleRows(2*Nr, Nr).minCoeff();
        rhoE_max = sn_set.middleRows(3*Nr, Nr).maxCoeff();
        rhoE_min = sn_set.middleRows(3*Nr, Nr).minCoeff();
        tke_max = sn_set.middleRows(4*Nr, Nr).maxCoeff();
        tke_min = sn_set.middleRows(4*Nr, Nr).minCoeff();
        omega_max = sn_set.middleRows(5*Nr, Nr).maxCoeff();
        omega_min = sn_set.middleRows(5*Nr, Nr).minCoeff();

    }
    else if (settings.ndim == 3 && nC == 7) {

        rho_max = sn_set.middleRows(0, Nr).maxCoeff();
        rho_min = sn_set.middleRows(0, Nr).minCoeff();
        rhoU_max = sn_set.middleRows(Nr, Nr).maxCoeff();
        rhoU_min = sn_set.middleRows(Nr, Nr).minCoeff();
        rhoV_max = sn_set.middleRows(2*Nr, Nr).maxCoeff();
        rhoV_min = sn_set.middleRows(2*Nr, Nr).minCoeff();
        rhoW_max = sn_set.middleRows(3*Nr, Nr).maxCoeff();
        rhoW_min = sn_set.middleRows(3*Nr, Nr).minCoeff();
        rhoE_max = sn_set.middleRows(4*Nr, Nr).maxCoeff();
        rhoE_min = sn_set.middleRows(4*Nr, Nr).minCoeff();
        tke_max = sn_set.middleRows(5*Nr, Nr).maxCoeff();
        tke_min = sn_set.middleRows(5*Nr, Nr).minCoeff();
        omega_max = sn_set.middleRows(6*Nr, Nr).maxCoeff();
        omega_min = sn_set.middleRows(6*Nr, Nr).minCoeff();

    }
    else {
        std::cout << "Combination of N_DIM and Conservative Variable not available " << std::endl << std::endl;
        exit(EXIT_FAILURE);
    }
}


Eigen::VectorXd IC ( Eigen::MatrixXd &sn_set, prob_settings settings, int nC, int Nr ) {

    double M = settings.Mach;
    double Re = settings.Re;
    double T = settings.T;
    double alpha = M_PI*settings.alpha/double(180);
    double beta =  M_PI*settings.beta/double(180);
    if ( settings.ndim == 2 ) beta = 0.0;

    double length = 1.0;
    double R = 287.058;
    double gamma = 1.4;
    double mu_ref = 1.716E-5;
    double T_ref = 273.15;
    double S = 110.4;

    double mu = mu_ref*std::pow(T/T_ref,1.5)*(T_ref + S)/(T + S);
    double V_magn = M*std::sqrt(gamma*R*T);
    double rho = Re*mu/(V_magn*length);
    double rhoU = rho*V_magn*std::cos(alpha)*std::cos(beta);
    double rhoV = rho*V_magn*std::sin(alpha);
    double rhoW = rho*V_magn*std::cos(alpha)*std::sin(beta);
    double rhoE = rho*(R/(gamma-1)*T + 0.5*V_magn*V_magn);

    //That only matters for turbulent calculation
    //since these values are used as default values in SU2, they are not present in config file but hard coded
    double n_turb = 0.05;
    double mu_turb2lam_ratio = 10.0;
    double nuTilde=0.1*mu/rho;

    double tke = 1.5*n_turb*n_turb*V_magn*V_magn;
    double omega = rho*tke/(std::max(mu*mu_turb2lam_ratio,1.e-25));
    // double rhotke = rho*tke;
    // double rhoomega = rho*omega;
    double rhotke = tke;
    double rhoomega = omega;

    Eigen::VectorXd Ic = Eigen::VectorXd::Zero(nC*Nr);

    if ( nC == 2 && settings.ndim == 2 )
    {
        Ic.head(Nr) = rhoU*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(Nr, Nr) = rhoV*Eigen::MatrixXd::Ones(Nr,1);
        for ( int it = 0; it < settings.Ns; it++ )
            sn_set.col(it) -= Ic;

    } else if ( nC == 3 && settings.ndim == 3 )
    {
        Ic.head(Nr) = rhoU*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(Nr, Nr) = rhoV*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(2*Nr, Nr) = rhoW*Eigen::MatrixXd::Ones(Nr,1);
        for ( int it = 0; it < settings.Ns; it++ )
            sn_set.col(it) -= Ic;

    } else if ( nC == 4 && settings.ndim == 2) //Laminar 2D Navier-Stokes
    {
        Ic.head(Nr) = rho*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(Nr, Nr) = rhoU*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(2*Nr, Nr) = rhoV*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(3*Nr, Nr) = rhoE*Eigen::MatrixXd::Ones(Nr,1);
        for ( int it = 0; it < settings.Ns; it++ )
            sn_set.col(it) -= Ic;

    } else if ( nC == 5 && settings.ndim == 3 ) //Laminar 3D Navier-Stokes
    {
        Ic.head(Nr) = rho*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(Nr, Nr) = rhoU*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(2*Nr, Nr) = rhoV*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(3*Nr, Nr) = rhoW*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(4*Nr, Nr) = rhoE*Eigen::MatrixXd::Ones(Nr,1);
        for ( int it = 0; it < settings.Ns; it++ )
            sn_set.col(it) -= Ic;

    } else if ( nC== 5 && settings.ndim == 2 ) // Turbolent 2D Spalart Allmaras
    {
        Ic.head(Nr) = rho*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(Nr, Nr) = rhoU*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(2*Nr, Nr) = rhoV*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(3*Nr, Nr) = rhoE*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(4*Nr, Nr) = nuTilde*Eigen::MatrixXd::Ones(Nr,1);
        for ( int it = 0; it < settings.Ns; it++ )
            sn_set.col(it) -= Ic;

    } else if ( nC == 6 && settings.ndim == 2 ) //Turbulent 2D SST
    {
        Ic.head(Nr) = rho*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(Nr, Nr) = rhoU*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(2*Nr, Nr) = rhoV*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(3*Nr, Nr) = rhoE*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(4*Nr, Nr) = rhotke*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(5*Nr, Nr) = rhoomega*Eigen::MatrixXd::Ones(Nr,1);
        for ( int it = 0; it < settings.Ns; it++ )
            sn_set.col(it) -= Ic;

    } else if ( nC == 7 ) //Turbulent 3D Navier-Stokes
    {
        Ic.head(Nr) = rho*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(Nr, Nr) = rhoU*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(2*Nr, Nr) = rhoV*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(3*Nr, Nr) = rhoW*Eigen::MatrixXd::Ones(Nr,1); //no sideslip angle
        Ic.segment(4*Nr, Nr) = rhoE*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(5*Nr, Nr) = rhotke*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(6*Nr, Nr) = rhoomega*Eigen::MatrixXd::Ones(Nr,1);
        for ( int it = 0; it < settings.Ns; it++ )
            sn_set.col(it) -= Ic;

    } else {
        std::cout << "Set well number of Variables for subtracting initial condition" << std::endl;
        exit(EXIT_FAILURE);
    }

    return Ic;
}


void Direct_Normalization(Eigen::MatrixXd &sn_set, const prob_settings &settings, const int nC, double &rho_max,
                          double &rho_min, double &rhoU_max, double &rhoU_min, double &rhoV_max, double &rhoV_min,
                          double &rhoW_max, double &rhoW_min, double &rhoE_max, double &rhoE_min, double &tke_min,
                          double &tke_max, double &omega_min, double &omega_max, double &nuTilde_min, double &nuTilde_max) {


    int Nr = sn_set.rows()/nC;
    if ( settings.ndim == 2 && nC == 4 ) {

        sn_set.middleRows(0, Nr) = (sn_set.middleRows(0, Nr) - Eigen::MatrixXd::Ones(Nr, settings.Ns)*rho_min )/(rho_max - rho_min);
        sn_set.middleRows(Nr, Nr) = (sn_set.middleRows(Nr, Nr) - Eigen::MatrixXd::Ones(Nr, settings.Ns)*rhoU_min)/(rhoU_max - rhoU_min);
        sn_set.middleRows(2*Nr, Nr) = (sn_set.middleRows(2*Nr, Nr) - Eigen::MatrixXd::Ones(Nr, settings.Ns)*rhoV_min)/(rhoV_max - rhoV_min);
        sn_set.middleRows(3*Nr, Nr) = (sn_set.middleRows(3*Nr, Nr) - Eigen::MatrixXd::Ones(Nr, settings.Ns)*rhoE_min)/(rhoE_max - rhoE_min);

    }
    else if ( settings.ndim == 2 && nC == 5 ) {

        sn_set.middleRows(0, Nr) = (sn_set.middleRows(0, Nr) - Eigen::MatrixXd::Ones(Nr, settings.Ns)*rho_min )/(rho_max - rho_min);
        sn_set.middleRows(Nr, Nr) = (sn_set.middleRows(Nr, Nr) - Eigen::MatrixXd::Ones(Nr, settings.Ns)*rhoU_min)/(rhoU_max - rhoU_min);
        sn_set.middleRows(2*Nr, Nr) = (sn_set.middleRows(2*Nr, Nr) - Eigen::MatrixXd::Ones(Nr, settings.Ns)*rhoV_min)/(rhoV_max - rhoV_min);
        sn_set.middleRows(3*Nr, Nr) = (sn_set.middleRows(3*Nr, Nr) - Eigen::MatrixXd::Ones(Nr, settings.Ns)*rhoE_min)/(rhoE_max - rhoE_min);
        sn_set.middleRows(4*Nr, Nr) = (sn_set.middleRows(4*Nr, Nr) - Eigen::MatrixXd::Ones(Nr, settings.Ns)*nuTilde_min)/(nuTilde_max - nuTilde_min);

    }
    else if ( settings.ndim == 2 && nC == 6 ) {

        sn_set.middleRows(0, Nr) = (sn_set.middleRows(0, Nr) - Eigen::MatrixXd::Ones(Nr, settings.Ns)*rho_min )/(rho_max - rho_min);
        sn_set.middleRows(Nr, Nr) = (sn_set.middleRows(Nr, Nr) - Eigen::MatrixXd::Ones(Nr, settings.Ns)*rhoU_min)/(rhoU_max - rhoU_min);
        sn_set.middleRows(2*Nr, Nr) = (sn_set.middleRows(2*Nr, Nr) - Eigen::MatrixXd::Ones(Nr, settings.Ns)*rhoV_min)/(rhoV_max - rhoV_min);
        sn_set.middleRows(3*Nr, Nr) = (sn_set.middleRows(3*Nr, Nr) - Eigen::MatrixXd::Ones(Nr, settings.Ns)*rhoE_min)/(rhoE_max - rhoE_min);
        sn_set.middleRows(4*Nr, Nr) = (sn_set.middleRows(4*Nr, Nr) - Eigen::MatrixXd::Ones(Nr, settings.Ns)*tke_min)/(tke_max - tke_min);
        sn_set.middleRows(5*Nr, Nr) = (sn_set.middleRows(5*Nr, Nr) - Eigen::MatrixXd::Ones(Nr, settings.Ns)*omega_min)/(omega_max - omega_min);

    }
    else if ( settings.ndim == 3 && nC == 7 ) {

        sn_set.middleRows(0, Nr) = (sn_set.middleRows(0, Nr) - Eigen::MatrixXd::Ones(Nr, settings.Ns)*rho_min )/(rho_max - rho_min);
        sn_set.middleRows(Nr, Nr) = (sn_set.middleRows(Nr, Nr) - Eigen::MatrixXd::Ones(Nr, settings.Ns)*rhoU_min)/(rhoU_max - rhoU_min);
        sn_set.middleRows(2*Nr, Nr) = (sn_set.middleRows(2*Nr, Nr) - Eigen::MatrixXd::Ones(Nr, settings.Ns)*rhoV_min)/(rhoV_max - rhoV_min);
        sn_set.middleRows(3*Nr, Nr) = (sn_set.middleRows(3*Nr, Nr) - Eigen::MatrixXd::Ones(Nr, settings.Ns)*rhoW_min)/(rhoW_max - rhoW_min);
        sn_set.middleRows(4*Nr, Nr) = (sn_set.middleRows(4*Nr, Nr) - Eigen::MatrixXd::Ones(Nr, settings.Ns)*rhoE_min)/(rhoE_max - rhoE_min);
        sn_set.middleRows(5*Nr, Nr) = (sn_set.middleRows(5*Nr, Nr) - Eigen::MatrixXd::Ones(Nr, settings.Ns)*tke_min)/(tke_max - tke_min);
        sn_set.middleRows(6*Nr, Nr) = (sn_set.middleRows(6*Nr, Nr) - Eigen::MatrixXd::Ones(Nr, settings.Ns)*omega_min)/(omega_max - omega_min);

    }
    else {
        std::cout << "Combination of N_DIM and Conservative Variable not available " << std::endl << std::endl;
        exit(EXIT_FAILURE);
    }
}


void Inverse_Normalization(Eigen::MatrixXd &sn_set, const prob_settings &settings, const int nC, double &rho_max,
                           double &rho_min, double &rhoU_max, double &rhoU_min, double &rhoV_max, double &rhoV_min,
                           double &rhoW_max, double &rhoW_min, double &rhoE_max, double &rhoE_min, double &tke_min,
                           double &tke_max, double &omega_min, double &omega_max, double &nuTilde_min, double &nuTilde_max) {


    int Nr = sn_set.rows()/nC;
    int ntimes = sn_set.cols();

    if (settings.ndim == 2 && nC == 4) {
        sn_set.middleRows(0, Nr) = sn_set.middleRows(0, Nr) * (rho_max - rho_min) + Eigen::MatrixXd::Ones(Nr, ntimes)*rho_min;
        sn_set.middleRows(Nr, Nr) = sn_set.middleRows(Nr, Nr) * (rhoU_max - rhoU_min) + Eigen::MatrixXd::Ones(Nr, ntimes)*rhoU_min;
        sn_set.middleRows(2 * Nr, Nr) = sn_set.middleRows(2 * Nr, Nr) * (rhoV_max - rhoV_min) + Eigen::MatrixXd::Ones(Nr, ntimes)*rhoV_min;
        sn_set.middleRows(3 * Nr, Nr) = sn_set.middleRows(3 * Nr, Nr) * (rhoE_max - rhoE_min) + Eigen::MatrixXd::Ones(Nr, ntimes)*rhoE_min;
    }
    else if (settings.ndim == 2 && nC == 5) {
        sn_set.middleRows(0, Nr) = sn_set.middleRows(0, Nr) * (rho_max - rho_min) + Eigen::MatrixXd::Ones(Nr, ntimes)*rho_min;
        sn_set.middleRows(Nr, Nr) = sn_set.middleRows(Nr, Nr) * (rhoU_max - rhoU_min) + Eigen::MatrixXd::Ones(Nr, ntimes)*rhoU_min;
        sn_set.middleRows(2 * Nr, Nr) = sn_set.middleRows(2 * Nr, Nr) * (rhoV_max - rhoV_min) + Eigen::MatrixXd::Ones(Nr, ntimes)*rhoV_min;
        sn_set.middleRows(3 * Nr, Nr) = sn_set.middleRows(3 * Nr, Nr) * (rhoE_max - rhoE_min) + Eigen::MatrixXd::Ones(Nr, ntimes)*rhoE_min;
        sn_set.middleRows(4 * Nr, Nr) = sn_set.middleRows(4 * Nr, Nr) * (nuTilde_max - nuTilde_min) + Eigen::MatrixXd::Ones(Nr, ntimes)*nuTilde_min;
    }
    else if (settings.ndim == 2 && nC == 6) {
        sn_set.middleRows(0, Nr) = sn_set.middleRows(0, Nr) * (rho_max - rho_min) + Eigen::MatrixXd::Ones(Nr, ntimes)*rho_min;
        sn_set.middleRows(Nr, Nr) = sn_set.middleRows(Nr, Nr) * (rhoU_max - rhoU_min) + Eigen::MatrixXd::Ones(Nr, 3)*rhoU_min;
        sn_set.middleRows(2 * Nr, Nr) = sn_set.middleRows(2 * Nr, Nr) * (rhoV_max - rhoV_min) + Eigen::MatrixXd::Ones(Nr, ntimes)*rhoV_min;
        sn_set.middleRows(3 * Nr, Nr) = sn_set.middleRows(3 * Nr, Nr) * (rhoE_max - rhoE_min) + Eigen::MatrixXd::Ones(Nr, ntimes)*rhoE_min;
        sn_set.middleRows(4 * Nr, Nr) = sn_set.middleRows(4 * Nr, Nr) * (tke_max - tke_min) + Eigen::MatrixXd::Ones(Nr, ntimes)*tke_min;
        sn_set.middleRows(5 * Nr, Nr) = sn_set.middleRows(5 * Nr, Nr) * (omega_max - omega_min) + Eigen::MatrixXd::Ones(Nr, ntimes)*omega_min;
    }
    else if (settings.ndim == 3 && nC == 7) {
        sn_set.middleRows(0, Nr) = sn_set.middleRows(0, Nr) * (rho_max - rho_min) + Eigen::MatrixXd::Ones(Nr, ntimes)*rho_min;
        sn_set.middleRows(Nr, Nr) = sn_set.middleRows(Nr, Nr) * (rhoU_max - rhoU_min) + Eigen::MatrixXd::Ones(Nr, 3)*rhoU_min;
        sn_set.middleRows(2 * Nr, Nr) = sn_set.middleRows(2 * Nr, Nr) * (rhoV_max - rhoV_min) + Eigen::MatrixXd::Ones(Nr, ntimes)*rhoV_min;
        sn_set.middleRows(3 * Nr, Nr) = sn_set.middleRows(3 * Nr, Nr) * (rhoW_max - rhoW_min) + Eigen::MatrixXd::Ones(Nr, ntimes)*rhoW_min;
        sn_set.middleRows(4 * Nr, Nr) = sn_set.middleRows(4 * Nr, Nr) * (rhoE_max - rhoE_min) + Eigen::MatrixXd::Ones(Nr, ntimes)*rhoE_min;
        sn_set.middleRows(5 * Nr, Nr) = sn_set.middleRows(5 * Nr, Nr) * (tke_max - tke_min) + Eigen::MatrixXd::Ones(Nr, ntimes)*tke_min;
        sn_set.middleRows(6 * Nr, Nr) = sn_set.middleRows(6 * Nr, Nr) * (omega_max - omega_min) + Eigen::MatrixXd::Ones(Nr, ntimes)*omega_min;

    }
    else {
        std::cout << "Combination of N_DIM and Conservative Variable not available " << std::endl << std::endl;
        exit(EXIT_FAILURE);
    }
}


Eigen::VectorXi Inverse_POS (const Eigen::MatrixXd &sn_set, int Nsamples) {

    int Ns = sn_set.cols();
    int Ndof = sn_set.rows();
    Eigen::VectorXd norm_sn_set = Eigen::VectorXd::Zero(Ns);
    for ( int i = 0; i < Ns; i++ )
        norm_sn_set(i) = sn_set.col(i).norm()/(double)Ndof;

    Eigen::VectorXd ysamples = Eigen::VectorXd::LinSpaced(Nsamples, norm_sn_set.minCoeff(), norm_sn_set.maxCoeff());
    std::vector<int> I_POS = {};

    for (int i = 0; i < Nsamples; i++) {
        Eigen::VectorXd ftime = norm_sn_set - Eigen::VectorXd::Ones(Ns)*ysamples(i);
        for ( int j = 0; j < Ns-1; j++) {
            if ( (ftime(j)*ftime(j+1)) < 0.0 )
                I_POS.push_back(j);
        }
    }

    std::sort(I_POS.begin(),I_POS.end());
    I_POS.erase(std::unique(I_POS.begin(),I_POS.end()),I_POS.end());

    Eigen::Map<Eigen::VectorXi> Ipos(I_POS.data(), I_POS.size());

    return Ipos;

}