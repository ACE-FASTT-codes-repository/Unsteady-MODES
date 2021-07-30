#include "Opt_struct.hpp"

#ifdef __USE_PAGMO

//Leave this blank:
std::pair<std::vector<double>, std::vector<double> > SPOD_Adapt_Samp::get_bounds() const {
    return {problemBounds_[0], problemBounds_[1]};
}

//Fitness Function:
std::vector<double> SPOD_Adapt_Samp::fitness(const std::vector<double> &variables) const {
    
    // Initializing meaningful quantities
    double Dt_cfd = m_settings.Dt_cfd;
    double D_Samp = (double)m_settings.Ds;
    int Np = m_sn_set.rows();
    int Ns = m_settings.Ns;

    Eigen::VectorXd norm_sn_set = Eigen::VectorXd::Zero(Ns);
    for ( int it = 0; it < Ns; it ++) 
        norm_sn_set(it) = m_sn_set.col(it).norm(); 

    //Generating the vector of column indices to read in the full snapshot matrix
    int N_Var = static_cast<int>(variables.size());

    std::vector<int> ci_vec(N_Var+2); //+2 takes into account first and last snapshots (always included in the sampling)
    ci_vec[0] = 0;
    ci_vec[N_Var+1] = Ns-1;

    for ( int iVar = 1; iVar<(N_Var+1); iVar++ ) 
        ci_vec[iVar] = static_cast<int>(std::round(variables[iVar-1]/(Dt_cfd*D_Samp)));

    std::sort(ci_vec.begin(),ci_vec.end());
    ci_vec.erase(std::unique(ci_vec.begin(),ci_vec.end()),ci_vec.end());

    N_Var = static_cast<int>(ci_vec.size())-2;
    Eigen::Map<Eigen::VectorXi> ci(ci_vec.data(), N_Var+2);

// for (int iappo = 0 ; iappo < N_Var+2; iappo++ )
//     std::cout << "civec[" << iappo << "] = " << ci_vec[iappo] << ", ci_eigen[" << iappo << "] = " << ci(iappo) << std::endl;

// std::cout << std::endl;
    Eigen::MatrixXd sn_set_test = indexing(m_sn_set, Eigen::ArrayXi::LinSpaced(Np,0,Np-1),ci);

    //Performing basis extraction
    Eigen::VectorXd lambda = Eigen::VectorXd::Zero(N_Var+2);
    Eigen::VectorXd K_pc = Eigen::VectorXd::Zero(N_Var+2);
    Eigen::MatrixXd eig_vec = Eigen::MatrixXd::Zero(N_Var+2,N_Var+2);

    Eigen::MatrixXd Phi = SPOD_basis( sn_set_test,
                        lambda, K_pc, eig_vec,
                        m_Nf,
                        m_settings.flag_bc, 
                        m_settings.flag_filter,  
                        m_settings.sigma);            
    int Nm = Phi.cols();

    //Computing projection error
    Eigen::MatrixXd Sig = Eigen::MatrixXd::Zero(Nm, Nm);

    for ( int iMode = 0; iMode < Nm; iMode++ )
        Sig(iMode,iMode) = std::sqrt(lambda(iMode));

    Eigen::MatrixXd ErrP_SPOD_map = Eigen::MatrixXd::Zero(Np, Ns);
    Eigen::VectorXd ErrP_SPOD_time = Eigen::VectorXd::Zero(Ns);

    Eigen::MatrixXd PhiTPhi = Phi.transpose()*Phi;
    Eigen::MatrixXd dumCoefs = Phi.transpose()*m_sn_set;

    ErrP_SPOD_map = m_sn_set - Phi*(PhiTPhi.inverse()*dumCoefs);

    for ( int it = 0; it < Ns; it++ ) {
        int count = 0;
        for ( int iP = 0; iP < Np; iP++ )
            ErrP_SPOD_time(it) += ErrP_SPOD_map(iP,it)*ErrP_SPOD_map(iP,it);
        
        ErrP_SPOD_time(it) = std::sqrt(ErrP_SPOD_time(it))/norm_sn_set(it);
    }     

    //Computing final value of the objective function
    std::vector<double> fitness_vector;
    fitness_vector.push_back(ErrP_SPOD_time.maxCoeff());

    return fitness_vector;

    }




std::pair<std::vector<double>, std::vector<double> > SPOD_Adapt_Samp_::get_bounds() const {
    return {problemBounds_[0], problemBounds_[1]};
}

//Fitness Function:
std::vector<double> SPOD_Adapt_Samp_::fitness(const std::vector<double> &variables) const {

    // Initializing meaningful quantities
    double Dt_cfd = m_settings.Dt_cfd;
    double D_Samp = (double)m_settings.Ds;
    int Np = m_sn_set.rows();
    int Ns = m_settings.Ns;

    Eigen::VectorXd norm_sn_set = Eigen::VectorXd::Zero(Ns);
    for ( int it = 0; it < Ns; it ++)
        norm_sn_set(it) = m_sn_set.col(it).norm();

    //Generating the vector of column indices to read in the full snapshot matrix
    int N_Var = static_cast<int>(variables.size());

    std::vector<int> ci_vec(N_Var+2); //+2 takes into account first and last snapshots (always included in the sampling)
    ci_vec[0] = 0;
    ci_vec[N_Var+1] = Ns-1;

    for ( int iVar = 1; iVar<(N_Var+1); iVar++ )
        ci_vec[iVar] = ci_vec[iVar-1] + static_cast<int>(round(variables[iVar-1]/(Dt_cfd*D_Samp)));

//    for (int iappo = 0 ; iappo < N_Var+2; iappo++ )
//        std::cout << "civec[" << iappo << "] = " << ci_vec[iappo] << ", ci_eigen[" << iappo << "] = " << ci(iappo) << std::endl;

    for ( int it = 0; it < ci_vec.size(); it++ ) {
        if (ci_vec[it] > (Ns-1)) ci_vec[it] = 0;
    }

    std::sort(ci_vec.begin(),ci_vec.end());
    ci_vec.erase(std::unique(ci_vec.begin(),ci_vec.end()),ci_vec.end());
    N_Var = static_cast<int>(ci_vec.size())-2;
    Eigen::Map<Eigen::VectorXi> ci(ci_vec.data(), N_Var+2);

    Eigen::MatrixXd sn_set_test = indexing(m_sn_set, Eigen::ArrayXi::LinSpaced(Np,0,Np-1),ci);

    //Performing basis extraction
    Eigen::VectorXd lambda = Eigen::VectorXd::Zero(N_Var+2);
    Eigen::VectorXd K_pc = Eigen::VectorXd::Zero(N_Var+2);
    Eigen::MatrixXd eig_vec = Eigen::MatrixXd::Zero(N_Var+2,N_Var+2);

    Eigen::MatrixXd Phi = SPOD_basis( sn_set_test,
                                      lambda, K_pc, eig_vec,
                                      m_Nf,
                                      m_settings.flag_bc,
                                      m_settings.flag_filter,
                                      m_settings.sigma);
    int Nm = Phi.cols();
    //Computing projection error
    Eigen::MatrixXd Sig = Eigen::MatrixXd::Zero(Nm, Nm);

    for ( int iMode = 0; iMode < Nm; iMode++ )
        Sig(iMode,iMode) = std::sqrt(lambda(iMode));

    Eigen::MatrixXd ErrP_SPOD_map = Eigen::MatrixXd::Zero(Np, Ns);
    Eigen::VectorXd ErrP_SPOD_time = Eigen::VectorXd::Zero(Ns);

    Eigen::MatrixXd PhiTPhi = Phi.transpose()*Phi;
    Eigen::MatrixXd dumCoefs = Phi.transpose()*m_sn_set;
//    ErrP_SPOD_map = m_sn_set - Phi*(PhiTPhi.inverse()*dumCoefs); //for SPOD with any filter value
    ErrP_SPOD_map = m_sn_set - Phi * dumCoefs; //Only for POD

    for ( int it = 0; it < Ns; it++ ) {
        int count = 0;
        for ( int iP = 0; iP < Np; iP++ )
            ErrP_SPOD_time(it) += ErrP_SPOD_map(iP,it)*ErrP_SPOD_map(iP,it);

        ErrP_SPOD_time(it) = std::sqrt(ErrP_SPOD_time(it))/norm_sn_set(it);
    }

    //Computing final value of the objective function
    std::vector<double> fitness_vector;
    fitness_vector.push_back(ErrP_SPOD_time.maxCoeff());
    double sum_variables = 0.0;
    for ( int it = 0; it < variables.size(); it++ ) sum_variables += variables[it];

    sum_variables -= Dt_cfd*(Ns-2);
//    std::cout << "Sum variables = " << sum_variables << std::endl;
//    std::cout << std::endl;
//    std::cout << std::endl;
    fitness_vector.push_back(sum_variables);
//    fitness_vector.push_back(static_cast<double>(ci_vec[N_Var] - (Ns-1)));
    return fitness_vector;

}




std::pair<std::vector<double>, std::vector<double> > DMD_Adapt_Samp::get_bounds() const {
    return {problemBounds_[0], problemBounds_[1]};
}

//Fitness Function:
std::vector<double> DMD_Adapt_Samp::fitness(const std::vector<double> &variables) const {

    // Initializing meaningful quantities
    double Dt_cfd = m_settings.Dt_cfd;
    double D_Samp = (double)m_settings.Ds;
    int Np = m_sn_set.rows();
    int Ns = m_settings.Ns;

    Eigen::VectorXd norm_sn_set = Eigen::VectorXd::Zero(Ns);
    for ( int it = 0; it < Ns; it ++)
        norm_sn_set(it) = m_sn_set.col(it).norm();

    //Generating the vector of column indices to read in the full snapshot matrix
    int N_Var = static_cast<int>(variables.size());

    std::vector<int> ci_vec(N_Var+2); //+2 takes into account first and last snapshots (always included in the sampling)
    ci_vec[0] = 0;
    ci_vec[N_Var+1] = Ns-2;

    for ( int iVar = 1; iVar<(N_Var+1); iVar++ )
        ci_vec[iVar] = static_cast<int>(std::round(variables[iVar-1]/(Dt_cfd*D_Samp)));

    std::sort(ci_vec.begin(),ci_vec.end());
    ci_vec.erase(std::unique(ci_vec.begin(),ci_vec.end()),ci_vec.end());

    N_Var = static_cast<int>(ci_vec.size())-2;
    Eigen::VectorXi t_pos(ci_vec.size());
    for ( int i = 0; i < ci_vec.size(); i++ )
        t_pos(i) = ci_vec[i];

    //Performing basis extraction
    Eigen::VectorXcd lambda_DMD = Eigen::VectorXd::Zero(N_Var+2);
    Eigen::MatrixXcd eig_vec_DMD = Eigen::MatrixXd::Zero(N_Var+2,N_Var+2);
    Eigen::VectorXd lambda_POD = Eigen::VectorXd::Zero(N_Var+2);
    Eigen::MatrixXd eig_vec_POD = Eigen::MatrixXd::Zero(N_Var+2,N_Var+2);

    Eigen::MatrixXcd Phi = DMD_Adaptive_basis( m_sn_set,
                            lambda_DMD,
                            eig_vec_DMD,
                            lambda_POD,
                            eig_vec_POD,
                            t_pos );

    int Nm = Phi.cols();

    //Computing projection error

    Eigen::MatrixXd ErrP_DMD_map = Eigen::MatrixXd::Zero(Np, Ns);
    Eigen::VectorXd ErrP_DMD_time = Eigen::VectorXd::Zero(Ns);

    Eigen::MatrixXcd PhiTPhi = Phi.transpose()*Phi;
    Eigen::MatrixXcd dumCoefs = Phi.transpose()*m_sn_set;
//    Eigen::MatrixXcd P_u = Phi*(PhiTPhi.inverse()*dumCoefs); //Using pseudo(Moore-Penrose)-inverse

    //Using linear system solving
    Eigen::MatrixXcd Coeffs = PhiTPhi.colPivHouseholderQr().solve(dumCoefs);
    Eigen::MatrixXcd P_u = Phi * Coeffs;

    ErrP_DMD_map = m_sn_set - P_u.real();

    for ( int it = 0; it < Ns; it++ ) {
        int count = 0;
        for ( int iP = 0; iP < Np; iP++ )
            ErrP_DMD_time(it) += ErrP_DMD_map(iP,it)*ErrP_DMD_map(iP,it);

        ErrP_DMD_time(it) = std::sqrt(ErrP_DMD_time(it))/norm_sn_set(it);
    }

    //Computing final value of the objective function
//    double sum = 0;
//    ErrP_DMD_time(0) = ErrP_DMD_time(0)/2.0;
//    ErrP_DMD_time(ErrP_DMD_time.size()-1) = ErrP_DMD_time(ErrP_DMD_time.size()-1)/2.0;
//    for ( int it = 0; it < ErrP_DMD_time.size() ; it++ )
//        sum += ErrP_DMD_time(it);
//    sum *= Dt_cfd*D_Samp;

    std::vector<double> fitness_vector;
    //Integral
    fitness_vector.push_back(ErrP_DMD_time.norm());
    //Maximum
//    fitness_vector.push_back(ErrP_DMD_time.maxCoeff());


    return fitness_vector;

}






















std::pair<std::vector<double>, std::vector<double> > SPOD_Adapt_Samp_Int::get_bounds() const {
    return {problemBounds_[0], problemBounds_[1]};
}

//Fitness Function:
std::vector<double> SPOD_Adapt_Samp_Int::fitness(const std::vector<double> &variables) const {

    // Initializing meaningful quantities
    int Np = m_sn_set.rows();
    int Ns = m_settings.Ns;

    Eigen::VectorXd norm_sn_set = Eigen::VectorXd::Zero(Ns);
    for ( int it = 0; it < Ns; it ++)
        norm_sn_set(it) = m_sn_set.col(it).norm();

    //Generating the vector of column indices to read in the full snapshot matrix
    int N_Var = static_cast<int>(variables.size());

    std::vector<int> ci_vec(N_Var+2); //+2 takes into account first and last snapshots (always included in the sampling)
    ci_vec[0] = 0;
    ci_vec[N_Var+1] = Ns-1;

    for ( int iVar = 1; iVar<(N_Var+1); iVar++ )
        ci_vec[iVar] = variables[iVar-1];

//    for (int iappo = 0 ; iappo < N_Var+2; iappo++ )
//        std::cout << "civec[" << iappo << "] = " << ci_vec[iappo] << ", ci_eigen[" << iappo << "] = " << ci(iappo) << std::endl;
//    for (int iappo = 0 ; iappo < N_Var+2; iappo++ )
//        std::cout << "civec[" << iappo << "] = " << ci_vec[iappo]<< std::endl;

//    std::sort(ci_vec.begin(),ci_vec.end());
    ci_vec.erase(std::unique(ci_vec.begin(),ci_vec.end()),ci_vec.end());
    N_Var = static_cast<int>(ci_vec.size())-2;
    Eigen::Map<Eigen::VectorXi> ci(ci_vec.data(), N_Var+2);

    Eigen::MatrixXd sn_set_test = indexing(m_sn_set, Eigen::ArrayXi::LinSpaced(Np,0,Np-1),ci);

    //Performing basis extraction
    Eigen::VectorXd lambda = Eigen::VectorXd::Zero(N_Var+2);
    Eigen::VectorXd K_pc = Eigen::VectorXd::Zero(N_Var+2);
    Eigen::MatrixXd eig_vec = Eigen::MatrixXd::Zero(N_Var+2,N_Var+2);

    Eigen::MatrixXd Phi = SPOD_basis( sn_set_test,
                                      lambda, K_pc, eig_vec,
                                      m_Nf,
                                      m_settings.flag_bc,
                                      m_settings.flag_filter,
                                      m_settings.sigma);
    int Nm = Phi.cols();

    Eigen::MatrixXd ErrP_SPOD_map = Eigen::MatrixXd::Zero(Np, Ns);
    Eigen::VectorXd ErrP_SPOD_time = Eigen::VectorXd::Zero(Ns);

    Eigen::MatrixXd PhiTPhi = Phi.transpose()*Phi;
    Eigen::MatrixXd dumCoefs = Phi.transpose()*m_sn_set;
    ErrP_SPOD_map = m_sn_set - Phi*(PhiTPhi.inverse()*dumCoefs);

    for ( int it = 0; it < Ns; it++ ) {

        for ( int iP = 0; iP < Np; iP++ )
            ErrP_SPOD_time(it) += ErrP_SPOD_map(iP,it)*ErrP_SPOD_map(iP,it);

        ErrP_SPOD_time(it) = std::sqrt(ErrP_SPOD_time(it))/norm_sn_set(it);
    }

    //Computing final value of the objective function
    std::vector<double> fitness_vector;
    fitness_vector.push_back(ErrP_SPOD_time.maxCoeff());

    return fitness_vector;

}






std::pair<std::vector<double>, std::vector<double> > SPOD_Adapt_Samp_Int_::get_bounds() const {
    return {problemBounds_[0], problemBounds_[1]};
}

//Fitness Function:
std::vector<double> SPOD_Adapt_Samp_Int_::fitness(const std::vector<double> &variables) const {

    // Initializing meaningful quantities
    int Np = m_sn_set.rows();
    int Ns = m_settings.Ns;

    Eigen::VectorXd norm_sn_set = Eigen::VectorXd::Zero(Ns);
    for ( int it = 0; it < Ns; it ++)
        norm_sn_set(it) = m_sn_set.col(it).norm();

    //Generating the vector of column indices to read in the full snapshot matrix
    int N_Var = static_cast<int>(variables.size());

    std::vector<int> ci_vec(N_Var+2); //+2 takes into account first and last snapshots (always included in the sampling)
    ci_vec[0] = 0;
    ci_vec[N_Var+1] = Ns-1;

    for ( int iVar = 1; iVar<(N_Var+1); iVar++ )
        ci_vec[iVar] = ci_vec[iVar-1] + variables[iVar-1];

//    for (int iappo = 0 ; iappo < N_Var+2; iappo++ )
//        std::cout << "civec[" << iappo << "] = " << ci_vec[iappo] << ", ci_eigen[" << iappo << "] = " << ci(iappo) << std::endl;
//    for (int iappo = 0 ; iappo < N_Var+2; iappo++ )
//        std::cout << "civec[" << iappo << "] = " << ci_vec[iappo]<< std::endl;

    for ( int it = 0; it < ci_vec.size(); it++ ) {
        if (ci_vec[it] > (Ns-1)) ci_vec[it] = 0;
    }
//    std::sort(ci_vec.begin(),ci_vec.end());
    ci_vec.erase(std::unique(ci_vec.begin(),ci_vec.end()),ci_vec.end());
    N_Var = static_cast<int>(ci_vec.size())-2;
    Eigen::Map<Eigen::VectorXi> ci(ci_vec.data(), N_Var+2);

    Eigen::MatrixXd sn_set_test = indexing(m_sn_set, Eigen::ArrayXi::LinSpaced(Np,0,Np-1),ci);

    //Performing basis extraction
    Eigen::VectorXd lambda = Eigen::VectorXd::Zero(N_Var+2);
    Eigen::VectorXd K_pc = Eigen::VectorXd::Zero(N_Var+2);
    Eigen::MatrixXd eig_vec = Eigen::MatrixXd::Zero(N_Var+2,N_Var+2);

    Eigen::MatrixXd Phi = SPOD_basis( sn_set_test,
                                      lambda, K_pc, eig_vec,
                                      m_Nf,
                                      m_settings.flag_bc,
                                      m_settings.flag_filter,
                                      m_settings.sigma);
    int Nm = Phi.cols();

    Eigen::MatrixXd ErrP_SPOD_map = Eigen::MatrixXd::Zero(Np, Ns);
    Eigen::VectorXd ErrP_SPOD_time = Eigen::VectorXd::Zero(Ns);

    Eigen::MatrixXd PhiTPhi = Phi.transpose()*Phi;
    Eigen::MatrixXd dumCoefs = Phi.transpose()*m_sn_set;
    ErrP_SPOD_map = m_sn_set - Phi*(PhiTPhi.inverse()*dumCoefs);

    for ( int it = 0; it < Ns; it++ ) {

        for ( int iP = 0; iP < Np; iP++ )
            ErrP_SPOD_time(it) += ErrP_SPOD_map(iP,it)*ErrP_SPOD_map(iP,it);

        ErrP_SPOD_time(it) = std::sqrt(ErrP_SPOD_time(it))/norm_sn_set(it);
    }

    //Computing final value of the objective function
    std::vector<double> fitness_vector;
    fitness_vector.push_back(ErrP_SPOD_time.maxCoeff());
    double sum_variables = 0.0;
    for ( int it = 0; it < variables.size(); it++ ) sum_variables += variables[it];

    sum_variables -= Ns-2;
//    std::cout << "Sum variables = " << sum_variables << std::endl;
//    std::cout << std::endl;
//    std::cout << std::endl;
    fitness_vector.push_back(sum_variables);
//    fitness_vector.push_back(static_cast<double>(ci_vec[N_Var] - (Ns-1)));

    return fitness_vector;

}



std::pair<std::vector<double>, std::vector<double> > DMD_Best_Modes::get_bounds() const {
    return {problemBounds_[0], problemBounds_[1]};
}

//Fitness Function:
std::vector<double> DMD_Best_Modes::fitness(const std::vector<double> &variables) const {

    int Np = m_sn_set.rows();
    int Ns = m_sn_set.cols();
    int N_Var = static_cast<int>(variables.size());
    std::vector<int> bool_mode(N_Var);
    Eigen::VectorXd norm_sn_set = Eigen::VectorXd::Zero(Ns);
    for ( int it = 0; it < Ns; it ++)
        norm_sn_set(it) = m_sn_set.col(it).norm();

    for ( int iMode = 0; iMode < N_Var; iMode++ ) bool_mode[iMode] = static_cast<int> (variables[iMode]);

    int Nm = 0;
    for ( int ibool = 0; ibool < N_Var; ibool++) Nm += bool_mode[ibool];
//    std::cout << "Nm = " << Nm << std::endl;
    Eigen::MatrixXcd Phi_Best = Eigen::MatrixXcd::Zero(Np,Nm);
    Eigen::MatrixXcd Coefs_Best = Eigen::MatrixXcd::Zero(Nm,Ns);
    std::vector<int> idx_Best = {};
    int iMBest = 0;
    for ( int iMode = 0; iMode < N_Var; iMode++ ){
        if ( bool_mode[iMode] == 1) {
            Phi_Best.col(iMBest) = m_Phi.col(iMode);
            Coefs_Best.row(iMBest) = m_Coefs.row(iMode);
            idx_Best.push_back(iMode);
            iMBest++;
        }
    }
//    std::cout << "In fitness " << std::endl;
//    if ( iMBest > Nm-1 ) std::cout << "Warning! Constraint not verified. Nm = " << iMBest << std::endl;

    //Computing projection error

    Eigen::MatrixXd ErrP_DMD_map = Eigen::MatrixXd::Zero(Np, Ns);
    Eigen::VectorXd ErrP_DMD_time = Eigen::VectorXd::Zero(Ns);

    Eigen::MatrixXcd P_u = Phi_Best * Coefs_Best;
    ErrP_DMD_map = m_sn_set - P_u.real();

    for ( int it = 0; it < Ns; it++ ) {
        int count = 0;
        for ( int iP = 0; iP < Np; iP++ )
            ErrP_DMD_time(it) += ErrP_DMD_map(iP,it)*ErrP_DMD_map(iP,it);

        ErrP_DMD_time(it) = std::sqrt(ErrP_DMD_time(it))/norm_sn_set(it);
    }

    std::vector<double> fitness_vector;
    //Integral
//    fitness_vector.push_back(ErrP_DMD_time.norm());
    //Maximum
    fitness_vector.push_back(ErrP_DMD_time.maxCoeff());

//    fitness_vector.push_back(Nm - m_settings.r);

    return fitness_vector;

}


std::pair<std::vector<double>, std::vector<double> > Best_Modes::get_bounds() const {
    return {problemBounds_[0], problemBounds_[1]};
}

//Fitness Function:
std::vector<double> Best_Modes::fitness(const std::vector<double> &variables) const {

    int Np = m_sn_set.rows();
    int Ns = m_sn_set.cols();
    int N_Var = static_cast<int>(variables.size());
    std::vector<int> bool_mode(N_Var);
    Eigen::VectorXd norm_sn_set = Eigen::VectorXd::Zero(Ns);
    for ( int it = 0; it < Ns; it ++)
        norm_sn_set(it) = m_sn_set.col(it).norm();

    for ( int iMode = 0; iMode < N_Var; iMode++ ) bool_mode[iMode] = static_cast<int> (variables[iMode]);

    int Nm = 0;
    for ( int ibool = 0; ibool < N_Var; ibool++) Nm += bool_mode[ibool];
//    std::cout << "Nm = " << Nm << std::endl;
    Eigen::MatrixXd Phi_Best = Eigen::MatrixXd::Zero(Np,Nm);
    Eigen::MatrixXd Coefs_Best = Eigen::MatrixXd::Zero(Nm,Ns);
    std::vector<int> idx_Best = {};
    int iMBest = 0;
    for ( int iMode = 0; iMode < N_Var; iMode++ ){
        if ( bool_mode[iMode] == 1) {
            Phi_Best.col(iMBest) = m_Phi.col(iMode);
            Coefs_Best.row(iMBest) = m_Coefs.row(iMode);
            idx_Best.push_back(iMode);
            iMBest++;
        }
    }
//    std::cout << "In fitness " << std::endl;
//    if ( iMBest > Nm-1 ) std::cout << "Warning! Constraint not verified. Nm = " << iMBest << std::endl;

    //Computing projection error

    Eigen::MatrixXd ErrP_DMD_map = Eigen::MatrixXd::Zero(Np, Ns);
    Eigen::VectorXd ErrP_DMD_time = Eigen::VectorXd::Zero(Ns);

    Eigen::MatrixXd P_u = Phi_Best * Coefs_Best;
    ErrP_DMD_map = m_sn_set - P_u;

    for ( int it = 0; it < Ns; it++ ) {
        int count = 0;
        for ( int iP = 0; iP < Np; iP++ )
            ErrP_DMD_time(it) += ErrP_DMD_map(iP,it)*ErrP_DMD_map(iP,it);

        ErrP_DMD_time(it) = std::sqrt(ErrP_DMD_time(it))/norm_sn_set(it);
    }

    std::vector<double> fitness_vector;
    //Integral
//    fitness_vector.push_back(ErrP_DMD_time.norm());
    //Maximum
    fitness_vector.push_back(ErrP_DMD_time.maxCoeff());

//    fitness_vector.push_back(Nm - m_settings.r);

    return fitness_vector;

}



#endif