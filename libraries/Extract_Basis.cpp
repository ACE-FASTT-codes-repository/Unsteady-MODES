#include "Extract_Basis.hpp"

int Nmod ( double En, Eigen::VectorXd K_pc )
{

    double En_content = 0.0;
    int count = 0;

    while ( En_content < En && count < K_pc.size() )
    {
        En_content = K_pc(count);
        count++;
    }

    return count;

}



void eig_sort( Eigen::VectorXd &lam, Eigen::MatrixXd &eig_vec ) {

    unsigned int swap_count = 1;
    double temp;
    Eigen::VectorXd temp_vec(eig_vec.rows());

    while (swap_count > 0)
    {

        swap_count = 0;

        for(unsigned int index = 1; index < lam.size(); index++)
        {

            if ( lam(index) > lam(index-1) )
            {

                temp = lam(index-1);
                lam(index-1) = lam(index);
                lam(index) = temp;

                temp_vec = eig_vec.col(index-1);
                eig_vec.col(index-1) = eig_vec.col(index);
                eig_vec.col(index) = temp_vec;

                swap_count++;

            }
        }
    }
}


void dmd_sort( Eigen::VectorXd &En, Eigen::MatrixXcd &Phi, Eigen::VectorXcd &lam, Eigen::VectorXcd &alfa ) {

    unsigned int swap_count = 1;
    double temp;
    std::complex<double> temp1, temp2;

    Eigen::VectorXcd temp_vec(Phi.rows()); 

    while (swap_count > 0)
    {

        swap_count = 0;

        for(unsigned int index = 1; index < En.size(); index++)
        {

            if ( En(index) > En(index-1) )
            {

                temp = En(index-1);
                En(index-1) = En(index);
                En(index) = temp;

                temp1 = alfa(index-1);
                alfa(index-1) = alfa(index);
                alfa(index) = temp1;

                temp2 = lam(index-1);
                lam(index-1) = lam(index);
                lam(index) = temp2;

                temp_vec = Phi.col(index-1);
                Phi.col(index-1) = Phi.col(index);
                Phi.col(index) = temp_vec;

                swap_count++;

            }
        }
    }
}



void dmd_tenvelope_sort( Eigen::MatrixXcd &Phi,
                         Eigen::VectorXcd &omega,
                         Eigen::VectorXcd &alfa,
                         std::vector<double> t_vec )
{
    double dt = t_vec[1] - t_vec[0];
    int Nm = Phi.cols();
    int Np = Phi.rows();
    int Nt = t_vec.size();
    Eigen::MatrixXd Tdyn = Eigen::MatrixXd::Zero(Nm,Nt);

    //Building matrix for time dynamics (absolute values)
    for ( int i = 0; i < Nm; i++ ){
        for ( int j = 0; j < Nt; j++ ){
            Tdyn(i,j) = std::abs(alfa(i)*std::exp(omega(i)*t_vec[j]));
        }
    }

//    std::cout << "Tdyn : \n " << Tdyn << std::endl;
    std::vector<int> idx_sort = {};

    //Ordering indices on the basis of time envelope

    int idxi;
    for ( int i_nm = 0; i_nm < Nm; i_nm++ ) {
        for (int i_time = 0; i_time < Nt; i_time++) {
            Eigen::VectorXd temp = Tdyn.col(i_time);
            temp.maxCoeff(&idxi);
            Tdyn(idxi, i_time) = 0.0;
            if (std::find(idx_sort.begin(), idx_sort.end(), idxi) != idx_sort.end()) {
                continue;
            } else {
                idx_sort.push_back(idxi);
            }
        }
    }
//    std::cout << "Tdyn should be all zero after the procedure, verify : \n " << Tdyn << std::endl << std::endl;
//    std::cout << "Modes DMD reordered with t-envelope:\n ";
//    for ( int it = 0; it < idx_sort.size(); it++ ) std::cout << idx_sort[it] << " ";
//
//    std::cout << std::endl;
    //Redistributing modes and coeffs on the basis of time envelope
    Eigen::VectorXcd temp_alfa = alfa;
    Eigen::VectorXcd temp_omega = omega;
    Eigen::MatrixXcd temp_Phi = Phi;
    for ( int i = 0; i < Nm; i++ ){
        alfa(i) = temp_alfa(idx_sort[i]);
        omega(i) = temp_omega(idx_sort[i]);
        Phi.col(i) = temp_Phi.col(idx_sort[i]);
    }

}


// Should work only with singular values ( Not eigenvalues problems )
int SVHT ( Eigen::VectorXd lam, int m, int n )
{

    double median;
    int n_sv = lam.size();

    if ( n_sv%2 != 0 )
        median = lam((n_sv-1)/2);
    else
        median = 0.5*(lam(n_sv/2) + lam(n_sv/2-1));

    double beta = (double)std::min(m,n)/(double)std::max(n,m);
    double omega = 0.56*std::pow(beta,3.0) - 0.95*std::pow(beta,2.0) 
                    + 1.82*beta + 1.43;
    double tau = omega*median;

    double eps = 1000;
    int Nm = 0;

    do
    {
        eps = lam(Nm);
        Nm++;
    } while ( eps > tau );

    return Nm-1;

}


int Not_zero( Eigen::VectorXd lam )
{
    int count = 0;
    double eps = 1.0;

    while ( eps > 1e-15 && count < lam.size() )
    {
        eps = lam(count);
        count ++;
    }

    if ( eps > 1e-15 )
        return count;
    
    return (count - 1);
}

Eigen::MatrixXd nullspace(Eigen::VectorXd s, Eigen::MatrixXd vh, const double atol, const double rtol )
{

    // Eigen::VectorXd s = svd.singularValues();
    // Eigen::MatrixXd vh = svd.matrixV();
    // eig_sort( s, vh);
    double tol = std::max(atol, rtol * s(0));
    
    int sum = 0;
    for ( int i = 0; i < s.size(); i++ )
    {
        if ( s(i) > tol )
            sum++;
    }

    Eigen::MatrixXd ns = vh.rightCols(vh.cols()-sum);
    
    return ns;
}



bool check_linear_consistency( Eigen::MatrixXd X, Eigen::MatrixXd Y, Eigen::MatrixXd Nullspace, const double atol, const double rtol )
{

    bool test;
    // Eigen::MatrixXd Nullspace = nullspace(X);
    int total = Nullspace.cols();
    int fails = 0;

    for ( int i = 0; i < total; i++ )
    { 
        Eigen::MatrixXd vec_check = Y*Nullspace.col(i);
        if ( vec_check.norm() > atol )
            fails += 1;
    }

    if ( fails > 0 )
    {   
        test = 0;
        std::cout << "Linear consistency check failed " << fails << " out of " << total << std::endl;
        std::cout << "Results might be misleading" << std::endl;
    }
    else
    {
        test = 1;
        std::cout << "Linear consistency ok!" << std::endl;
    }

    return test; 

}


std::string method_selected ( int n, int &Nf_SPOD, std::vector<int> Nf )
{
    if ( Nf.size() == 1 ) { //Only POD
        if ( n == 0 )
            return "SPOD";
        else if ( n == 1)
            return "DMD";
        else if ( n == 2 )
            return "RDMD";
        else {
            std::cout << "Bad value for index " << std::endl;
            return "NONE";
        }
    } else {
        if ((n > -1) && (n < (Nf.size()))) {
            for (int i = 0; i < Nf.size(); i++) {
                if (n == i)
                    Nf_SPOD = Nf[i];
            }
            return "SPOD";
        } else if (n == Nf.size())
            return "DMD";
        else if (n == (Nf.size() + 1))
            return "RDMD";
        else {
            std::cout << "Bad value for index " << std::endl;
            return "NONE";
        }
    }

    return "NONE";
}



Eigen::MatrixXd SPOD_basis( const Eigen::MatrixXd &snap_set,
                                Eigen::VectorXd &lam,
                                Eigen::VectorXd &K_pc,
                                Eigen::MatrixXd &eig_vec,
                                const int Nf,    
                                std::string bc_flag , 
                                std::string filter_flag ,  
                                double sigma )
{

    int count;
    int Nr = snap_set.rows();
    int Ns = snap_set.cols();

    Eigen::MatrixXd R_f(Ns, Ns);
    Eigen::MatrixXd phi_c(Nr, Ns);
    Eigen::VectorXd mean(Nr);
    Eigen::VectorXd g(2*Nf+1);

    R_f.setZero(Ns, Ns);

    if ( Nf == 0) 
    {                           //Calculating R POD
        R_f = snap_set.transpose()*snap_set;

    } else
    {                                      //Calculating R SPOD
    
        if ( filter_flag == "BOX" )
        {             //with BOX filter  

            for (int k = 0; k <= 2*Nf; k++ )
                g(k) = 1.0/(2.0*(double)Nf+1.0); 

        } else if ( filter_flag == "GAUSSIAN" )
        {             //with GAUSSIAN filter
            
            if ( sigma == 0.0)
            {
                std::cout << "sigma = 0 then only POD could be performed" << std::endl;
                exit (EXIT_FAILURE);
            }

            for ( int k = 0; k <= 2*Nf; k++ )
                g(k) = exp(-k*k/(2.0*sigma*sigma));
            

        } else
        {
            std::cout << "Filter selected not available" << std::endl;
            exit (EXIT_FAILURE);
        }

        Eigen::MatrixXd R = snap_set.transpose()*snap_set;

        //Zero-padded boundary conditions
        if (bc_flag == "ZERO"){

            for ( int i = 0; i < Ns; i++ )
            {
                for ( int j = 0; j < Ns; j++ )
                {

                    count = 0;

                    for ( int k = -Nf; k <= Nf; k++ ){
                        if ( i + k < 0 || j + k < 0 || i + k >= Ns || j + k >= Ns )
                            R_f(i,j) += 0.0;
                        else
                            R_f(i,j) += g(count)*R(i+k, j+k);
                        
                        count++; 
        
                    }
                }
            }
        } else 
        {
            std::cout << "Boundary condition not implemented " << std::endl;
            exit (EXIT_FAILURE);
        }
    }


    Eigen::EigenSolver<Eigen::MatrixXd> es(R_f); 
    lam = es.eigenvalues().real();
    eig_vec = es.eigenvectors().real();
    eig_sort( lam, eig_vec);

    double sum = 0;

    for (int i = 0; i < Ns; i++){
        sum += lam(i)/lam.sum();
        K_pc(i) = sum;
    }

    double tol = lam(0)*1e-12;
    // double tol = 1e-16;
    phi_c = snap_set*eig_vec;

    count = 0;
    while ( count < lam.size() && lam(count) > tol)
            count++;

    Eigen::MatrixXd phi(Nr,count);
    for ( int i = 0 ; i < count ; i++ )
        phi.col(i) = phi_c.col(i)/sqrt(lam(i));


    return phi;   


}



Eigen::MatrixXcd DMD_basis ( const Eigen::MatrixXd &snap_set,
                            Eigen::VectorXcd &lam,
                            Eigen::MatrixXcd &eig_vec,
                            Eigen::VectorXd &lam_POD,
                            Eigen::MatrixXd &eig_vec_POD,
                            const int r )
{   

    int Ns = snap_set.cols() - 1;
    int Nm;
    //                                Eigen::MatrixXd U = SPOD_basis(snap_set.leftCols(Ns), lam_POD, K_pc, eig_vec_POD );
    Eigen::BDCSVD<Eigen::MatrixXd> svd( snap_set.leftCols(Ns), 
                                        Eigen::ComputeThinU | Eigen::ComputeThinV );

    lam_POD = svd.singularValues();
    eig_vec_POD = svd.matrixV();
    eig_sort(lam_POD, eig_vec_POD);
    Eigen::MatrixXd Nullspace = nullspace( lam_POD, eig_vec_POD );
    bool test = check_linear_consistency( snap_set.leftCols(Ns), snap_set.rightCols(Ns), Nullspace );


    Eigen::MatrixXd U = svd.matrixU();                         
    Eigen::MatrixXd Sig_inv = Eigen::MatrixXd::Zero(U.cols(), U.cols());
    //std::cout << "Number of non-zero modes : " << U.cols() << std::endl;

    // for ( int i = 0; i < U.cols(); i++ )
        // Sig_inv(i, i) = 1.0/std::sqrt(lam_POD(i)); 

    for ( int i = 0; i < U.cols(); i++ )
        Sig_inv(i, i) = 1.0/lam_POD(i); 

    // int Nm = Nmod( En, K_pc );

    // if ( En == 1.0 )
    //     Nm = U.cols();

    if ( r == 0)
    {
        Nm = SVHT ( lam_POD, Ns, snap_set.rows() );
        std::cout << "DMD-rank from SVHT : " << Nm << std::endl;
    }
    else if ( r > 0 )
    {                    
        Nm = std::min(r, Not_zero ( lam_POD ));
        std::cout << "DMD user-defined rank : " << Nm << std::endl;
    }
    else
    {
        Nm = Not_zero ( lam_POD );
        std::cout << "DMD rank based on non zero singular values : " << Nm << std::endl;
    }
    
    Eigen::MatrixXcd phi = Eigen::MatrixXcd::Zero(snap_set.rows(), Nm);
    if ( Nm == 0)
    {
        std::cout << " Rank is zero, bad data to define low order evolution " << std::endl;
        lam = Eigen::VectorXd::Zero(0);
        eig_vec = Eigen::MatrixXd::Zero(0,0);

        return phi;

    }

    // std::cout << "Singular Values : \n" << lam_POD.head(Nm) << std::endl;
    Eigen::MatrixXd Atilde = U.leftCols(Nm).transpose()*snap_set.rightCols(Ns)*
                                eig_vec_POD.leftCols(Nm)*Sig_inv.block(0,0,Nm,Nm);
    
    // std::cout << "Atilde :\n " << Atilde << std::endl << std::endl; 

    //Eigen::VectorXcd Full_lam;
    if ( Atilde.size() == 1 && std::abs(Atilde(0,0)) < 1e-15 )
    {
        lam = Eigen::VectorXcd::Zero(1);
        eig_vec = Eigen::MatrixXcd::Ones(1,1);
        return phi;
    }
    else 
    {
        Eigen::EigenSolver<Eigen::MatrixXd> es(Atilde); 
        // Full_lam = es.eigenvalues();
        lam = es.eigenvalues();
        eig_vec = es.eigenvectors();
    }

    Eigen::MatrixXcd appo = snap_set.rightCols(Ns)*
                            eig_vec_POD.leftCols(Nm)*Sig_inv.block(0,0,Nm,Nm);                              

    
    //This is not Working well!!!
    //-----------Return only modes with frequencies : f < f_sampling/10 with f_sampling reading snapshots frequency------
    //Attenzione!!! Anche in tal caso la matrice di autovettori  rimane sempre la stessa (N_rank_DMDxN_rank_DMD)
    //quindi il calcolo dei coefficienti va effettuato tramite la funzione Calculate_Coefficients_DMD_Exact
    // Eigen::MatrixXcd Phi = appo*eig_vec;
    // std::complex<double> omegaj;

    // std::vector<int> idx_ny_modes = {};
    // std::cout << "Full Lam : \n" << Full_lam << std::endl;
    // for ( int i = 0; i < Full_lam.size(); i++ )
    // {
    //     omegaj = std::log(Full_lam(i));

    //     if ( (std::abs(omegaj.imag())) <= 1.0/8.0 )
    //         idx_ny_modes.push_back(i);
    // }

    // Eigen::MatrixXcd Phi_nyq(snap_set.rows(),idx_ny_modes.size());
    // lam = Eigen::VectorXcd::Zero(idx_ny_modes.size());

    // for ( int i = 0; i < idx_ny_modes.size(); i++ )
    // {
    //     Phi_nyq.col(i) = Phi.col(idx_ny_modes[i]);
    //     lam(i) = Full_lam(idx_ny_modes[i]);
    // } 

    // std::cout << "Final number of non-spurious modes (f < f_sampling/10) : " << Phi_nyq.cols() << std::endl;

    // return Phi_nyq;


    //This is the only choice for now
    //------------Return Modes considering also spurious ones (f > Nyquist frequency)----------------
    // Non divido per lambda
    //Se non divido per lambda mi trovo un timestep avanti (perchè? non so ancora)
    // return appo*eig_vec;

    // Divido per lambda

                                    
    for (int i = 0; i < Nm; i++)
    {
        phi.col(i) = 1.0/lam(i)*appo*eig_vec.col(i);
        phi.col(i) = phi.col(i)/phi.col(i).norm();
    }
    // 

    return phi;

    //Standard DMD
    // return U*eig_vec;
}


//The vector of lam and eig_vec DMD has to contain only the selected modes for reconstruction
Eigen::VectorXcd Calculate_Coefs_DMD ( const Eigen::MatrixXcd &eig_vec,
                                    const Eigen::MatrixXcd &eig_vec_POD,
                                    const Eigen::VectorXcd &lam,
                                    const Eigen::VectorXcd &lam_POD,
                                    const int Ns)   //Number of original snapshots minus 1
{
    
    int Nm = lam.size();
    Eigen::MatrixXcd Sigma = Eigen::MatrixXcd::Zero(Nm, Nm);
    Eigen::MatrixXcd V_and(Nm, Ns); 

    for ( int i = 0; i < Nm; i++)
        Sigma(i, i) = lam_POD(i);

    for ( int i = 0; i < Nm; i++ )
    {
        for ( int j = 0; j < Ns; j++ )
            V_and(i,j) = std::pow(lam(i), j);                                                                                         
    }

    Eigen::MatrixXcd Y_sq = eig_vec.conjugate().transpose()*eig_vec;
    Eigen::MatrixXcd V_and_sq = V_and*V_and.conjugate().transpose();
    V_and_sq = V_and_sq.conjugate();

    Eigen::MatrixXcd M1(Nm, Nm);

    for ( int i = 0; i < Nm; i++ )
    {

        for ( int j = 0; j < Nm; j++)
            M1(i,j) = Y_sq(i,j)*V_and_sq(i,j);

    } 

    Eigen::MatrixXcd M2 = V_and*eig_vec_POD.leftCols(Nm)*Sigma.conjugate()*eig_vec;
    Eigen::VectorXcd dum = M2.diagonal();
    dum = dum.conjugate();

    return M1.inverse()*dum;
    

} 


Eigen::VectorXcd Calculate_Coefs_DMD_exact ( const Eigen::MatrixXd &sn_set,  //matrix of first Ns-1 snaps 
                                            const Eigen::VectorXcd &lam,  //slow eigenvalues
                                            const Eigen::MatrixXcd &Phi ) //slow exact DMD modes
{

    int Nm = lam.size();
    int Ns = sn_set.cols();
    Eigen::MatrixXcd V_and(Nm, Ns);

    for ( int i = 0; i < Nm; i++ )
    {
        for ( int j = 0; j < Ns; j++ )
            V_and(i,j) = std::pow(lam(i), j);                                                                                         
    }

    Eigen::MatrixXcd Phi_sq = Phi.transpose().conjugate()*Phi;
    Eigen::MatrixXcd V_sq = V_and*V_and.conjugate().transpose();
    Eigen::MatrixXcd V_sqc = V_sq.conjugate();

    Eigen::MatrixXcd P = Phi_sq.cwiseProduct(V_sqc);
    Eigen::MatrixXcd appo = V_and*sn_set.transpose()*Phi;
    Eigen::VectorXcd dum = appo.diagonal();
    Eigen::VectorXcd dumc = dum.conjugate();

    return P.inverse()*dumc;

}


Eigen::MatrixXcd Calculate_Coefs_Matrix_DMD ( const Eigen::MatrixXd &sn_set,
                                            const Eigen::MatrixXcd &Phi,
                                            const Eigen::VectorXcd &omega,
                                            const double t_0,
                                            const double dt_dmd)
{
    int Ns = sn_set.cols();
    int Np = Phi.rows();
    int r = Phi.cols();
    Eigen::VectorXd t(Ns);
    t(0) = t_0;

    for ( int i = 1; i < Ns; i++ )
        t(i) = t(i-1) + dt_dmd;

    Eigen::MatrixXcd alfas(Ns,r);

    for ( int i = 0; i < Ns; i++ )
    {
        std::cout << "Solving LS at time step :" << i << std::endl;
        Eigen::MatrixXcd Phit(Np,r);
        for ( int j = 0; j < r; j++ )
            Phit.col(j) = Phi.col(j)*std::exp(omega(j)*t(i));

        Eigen::VectorXcd b = Eigen::VectorXcd::Zero(Np);
        for ( int k = 0; k < Np; k++ )
            b(k).real(sn_set(k,i));

        alfas.row(i) = Phit.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
        // std::cout << "Coefs : " << alfas.row(i).transpose() << std::endl << std::endl;
    }

    return alfas;

}


std::vector<node_mrDMD> mrDMD_basis( Eigen::MatrixXd &snap_set,
                                    std::vector<node_mrDMD> &nodes,                          
                                    const int r,
                                    double dts,                                                
                                    double t_0,
                                    int level,
                                    int bin_num,
                                    int offset,
                                    int max_levels,
                                    int max_cycles,
                                    std::string flag_coefs)
{
    std::cout << "--------LEVEL " << level << "---------------" << std::endl << std::endl;                                        
    const double PI = 3.1415926535;
    int nyq = 1*max_cycles;                         //number of snaps needed to capture cycles

    int N = snap_set.rows();                        //dimension of each snapshot
    int bin_size = snap_set.cols();                 //number of total snapshots available for the particular time bin

    if ( bin_size < nyq )
    {
        std::cout << "Max resolution possible reached ..." <<
                    "\n Returning to above level" << std::endl << std::endl;
        return nodes;
    }
    int step = std::floor( (double)bin_size/(double)(nyq - 1) );
    Eigen::MatrixXd _snap_set(N, nyq);

    for ( int i = 0, jj = 0; i < bin_size; i += step, jj++ )
        _snap_set.col(jj) = snap_set.col(i);

    Eigen::VectorXd lam_POD;
    Eigen::VectorXcd lam_DMD;
    Eigen::MatrixXd eig_vec_POD;
    Eigen::MatrixXcd eig_vec_DMD;

    //Perform pure DMD in the time bin
    Eigen::MatrixXcd Phi = DMD_basis( _snap_set,
                                        lam_DMD,
                                        eig_vec_DMD,
                                        lam_POD,
                                        eig_vec_POD,
                                        r );

    //select only the slow modes (only if max_level > 0 (not pure DMD))
    double rho = (double)max_cycles/(double)bin_size;
    if ( max_levels == 0 )
        rho = 1e6;

    std::vector<int> slow_idx = {};

    for ( int i =0; i < lam_DMD.size(); i++ )
    {

        if ( std::abs(std::log(lam_DMD(i))/(2.0*PI*(double)step)) <= rho )  //define PI
            slow_idx.push_back(i);
    }

    int n = slow_idx.size();

    
    Eigen::MatrixXcd Modes = Eigen::MatrixXcd::Zero(N,n);
    Eigen::VectorXcd lam_slw = Eigen::VectorXcd::Zero(n) ;
    Eigen::MatrixXcd eig_vec_slw = Eigen::MatrixXcd::Zero(eig_vec_DMD.rows(), n);
    Eigen::VectorXcd b_opt;
    Eigen::MatrixXcd alfas;    
    Eigen::MatrixXcd Psi = Eigen::MatrixXcd::Zero(n, bin_size);           //time evolution matrix
    Eigen::MatrixXcd D_dmd;
    Eigen::VectorXcd omega;

    //Calculate the correspondent coefficients and time evolution
    if ( n > 0 )
    {

        Eigen::MatrixXcd V_and(lam_slw.size(), bin_size);

        for ( int i = 0; i < n ; i++ )
        {

            Modes.col(i) = Phi.col(slow_idx[i]);                                                
            lam_slw(i) = lam_DMD(slow_idx[i]);
            eig_vec_slw.col(i) = eig_vec_DMD.col(slow_idx[i]);
            omega = Eigen::VectorXcd::Zero(lam_slw.size());

            for ( int rn = 0; rn < lam_slw.size(); rn++ )
            {
                omega(rn) = std::log(lam_slw(rn))/(dts*(double)step);
            }
        }

        b_opt = Calculate_Coefs_DMD_exact ( snap_set.leftCols(bin_size - 1),
                                            lam_slw,
                                            Modes );

        // alfas = Calculate_Coefs_Matrix_DMD ( snap_set.leftCols(bin_size - 1),
        //                                     Modes,
        //                                     omega,
        //                                     0.0,
        //                                     dts*step); 

        for ( int i = 0; i < lam_slw.size(); i++ )
        {
            for ( int j = 0; j < bin_size; j++ )
                V_and(i,j) = std::pow(std::pow(lam_slw(i), 1.0/(double)step), (double)j);                                                                                         
        }                                                                                    
//-----------Remove after check----------------------------------------------------------
// std::cout << " size V_and : [" << V_and.rows() << ", " << V_and.cols() << "]" << std::endl;
// std::cout << " size b_opt : [" << b_opt.rows() << ", " << b_opt.cols() << "]" << std::endl;

        for ( int i = 0; i < bin_size; i++ )
        {
            // if ( flag_coefs == "OPT" )
            // {
            Psi.col(i) = b_opt.cwiseProduct(V_and.col(i));
            // }
            // else if ( flag_coefs == "HYBRID" )
            // {
            //     Eigen::VectorXcd b = alfas.row(i);
            //     Psi.col(i) = b.cwiseProduct(V_and.col(i));
            // }
        }
        D_dmd = Modes*Psi;
        
//------------Remove after check-------------------------------------------------------------
// std::cout << " size Modes : [" << Modes.rows() << ", " << Modes.cols() << "]" << std::endl;
// std::cout << " size Psi : [" << Psi.rows() << ", " << Psi.cols() << "]" << std::endl;                                     
// std::cout << " size D_dmd : [" << D_dmd.rows() << ", " << D_dmd.cols() << "]" << std::endl;
    }
    else        //initialize b_opt  and Psi as empty vector and matrix
    {

        b_opt = Eigen::VectorXcd::Zero(lam_slw.size());
        // alfas = Eigen::MatrixXcd::Zero(0,0);
        Psi = Eigen::MatrixXcd::Zero(lam_slw.size(),bin_size);
        D_dmd = Eigen::MatrixXcd::Zero(N,bin_size);

    }

    //Subtracting influence of slow modes
    Eigen::MatrixXd Sub = D_dmd.real();
    snap_set = snap_set - Sub;

    //Storing all the necessary information in the node
    node_mrDMD node;
    node.l = level;
    node.bin_num = bin_num;
    node.bin_size = bin_size;      
    node.start = offset;           
    node.stop = offset + bin_size - 1; 
    node.step = step;          
    node.rho = rho;                
    node.r = lam_DMD.size();                    
    node.n = n;                    
    node.t_begin = t_0;
    node.t_end = t_0 + (bin_size - 1)*dts;
    node.dt = dts*step;
    node.lam = lam_slw;                  
    node.Modes = Modes;                
    node.Psi = Psi;         //time-dynamics transpose               
    node.Coefs = b_opt; 

    //Output nodes information
    std::cout << "----> Node  Info : " << std::endl << std::endl;
    std::cout <<  "Time interval : [" << node.t_begin << 
                    ", " << node.t_end << "]" << std::endl;
    std::cout << "Singular values : \n " << lam_POD << std::endl;    
    std::cout << " Snapshot interval (snaps index) : " << 
                    node.start << "\t" << node.stop << std::endl;
    std::cout << " bin_size : " << node.bin_size << std::endl;
    std::cout << " bin_num : " << node.bin_num << std::endl;
    std::cout << " DMD-rank : " << node.r << std::endl;
    std::cout << " Number of slow modes : " << node.n << std::endl;
    std::cout << std::endl << std::endl;

    nodes.push_back(node);

    //Apply DMD recursively
    int split = ceil(bin_size/2);

    if ( level < max_levels && split > 1 )
    {    

        Eigen::MatrixXd snap_set1;
        Eigen::MatrixXd snap_set2;
        
        if ( bin_size%2 != 0 )
        {
            snap_set1 = snap_set.leftCols(split + 1);
            nodes = mrDMD_basis( snap_set1,
                                    nodes,                         
                                    r,
                                    dts,
                                    t_0,
                                    level+1,
                                    2*bin_num,
                                    offset,
                                    max_levels,
                                    max_cycles,
                                    flag_coefs);
            snap_set2 = snap_set.rightCols(split + 1);
            nodes = mrDMD_basis( snap_set2, 
                                nodes,                         
                                r,
                                dts,
                                t_0 + (double)split*dts,
                                level+1,
                                2*bin_num+1,
                                offset+split,
                                max_levels,
                                max_cycles,
                                flag_coefs);

        }
        else
        {
            
            snap_set1 = snap_set.leftCols(split + 1);
            nodes = mrDMD_basis( snap_set1,
                                    nodes,                         
                                    r,
                                    dts,
                                    t_0,
                                    level+1,
                                    2*bin_num,
                                    offset,
                                    max_levels,
                                    max_cycles,
                                    flag_coefs);

            snap_set2 = snap_set.rightCols(split);
            nodes = mrDMD_basis( snap_set2,
                                nodes,                         
                                r,
                                dts,
                                t_0 + (double)(split)*dts,
                                level+1,
                                2*bin_num+1,
                                offset+split,
                                max_levels,
                                max_cycles,
                                flag_coefs);

        }
    }

    return nodes;

}



Eigen::MatrixXd RDMD_modes_coefs ( const Eigen::MatrixXd &sn_set,
                                    Eigen::MatrixXd &Coefs,
                                    Eigen::VectorXd &lambda,
                                    Eigen::VectorXd &K_pc,
                                    const int r,
                                    int &rdmd,
                                    double En )
{
    int Np = sn_set.rows();
    int Ns = sn_set.cols();

    Eigen::MatrixXd Phi_RDMD = Eigen::MatrixXd::Zero(Np, Ns);
    Eigen::MatrixXd res_set = sn_set;
    Eigen::VectorXd lam_POD;
    Eigen::VectorXcd lam_DMD;
    Eigen::MatrixXd eig_vec_POD;
    Eigen::MatrixXcd eig_vec_DMD;

    Eigen::VectorXd residual_time(Np);
    Eigen::VectorXd residual_time_norm(Ns);

    double eps = 0.0;
    Eigen::VectorXd svd_new = Eigen::VectorXd::Zero(Ns);
    Eigen::VectorXd svd_old = Eigen::VectorXd::Zero(Ns);
    double count = 0;

    if ( rdmd > Ns ) {
        std::cout << "Rank RDMD too high for the number of snapshots available. Resetting it to maximum value admissable" << std::endl;
        rdmd = Ns; 
    }    

    if ( rdmd == 0 ) {
        while ( eps < En && count < Ns+1 ) {
            //Perform pure DMD
            Eigen::MatrixXcd Phi = DMD_basis( res_set,
                                                lam_DMD,
                                                eig_vec_DMD,
                                                lam_POD,
                                                eig_vec_POD,
                                                r );
            
            if ( Phi.cols() == 0 )
                break;

            if ( count == 0)
                svd_old = lam_POD.cwiseProduct(lam_POD);

            svd_new = lam_POD.cwiseProduct(lam_POD);
            eps = (svd_old.sum() - svd_new.sum())/svd_old.sum();

            Eigen::MatrixXd Phi_r = Phi.real();
            Eigen::MatrixXd coef_mod(Phi.cols(),Ns);

            for ( int j = 0; j < Phi.cols(); j++ ) {
                double sum = 0.0;
                for ( int k = 0; k < Np; k++ )
                    sum += Phi_r(k,j)*Phi_r(k,j);

                if ( std::abs(sum) < 1e-15 ) {
                    Phi_r.col(j) = Eigen::MatrixXd::Zero(Np,1);
                } else {
                    Phi_r.col(j) = Phi_r.col(j)/std::sqrt(sum);
                }
            }

            Eigen::VectorXd residual_average(Phi.cols());
            int min_idx;

            for ( int r_dmd = 0; r_dmd < Phi.cols(); r_dmd++ ) {
                for ( int nt = 0; nt < Ns; nt++ ) {
                    coef_mod(r_dmd, nt) = res_set.col(nt).transpose()*Phi_r.col(r_dmd);
                    residual_time = res_set.col(nt) - coef_mod(r_dmd, nt)*Phi_r.col(r_dmd);
                    residual_time_norm(nt) = residual_time.norm();
                }
                double mean = 0.0;
                for ( int m = 0; m < Ns; m++ )
                    mean += residual_time_norm(m); 

                residual_average(r_dmd) = mean/Ns;
            }

            double min_Val = residual_average.minCoeff( &min_idx );

            Phi_RDMD.col(count) = Phi_r.col(min_idx);
            Coefs.row(count) = coef_mod.row(min_idx);
            lambda(count) = lam_DMD(min_idx).real();
            res_set = res_set - Phi_RDMD.col(count)*Coefs.row(count);
            std::cout << "Energy level at iteration " << count << " : " << std::setprecision(12) << eps*100 << "%" << std::endl; 

            count ++;
        }
        rdmd = count;

    } else {
         for ( int i = 0; i <= rdmd; i++ ){ //Considering also i = rdmd only for computing energy at the last iteration
            //Perform pure DMD
            Eigen::MatrixXcd Phi = DMD_basis( res_set,
                                                lam_DMD,
                                                eig_vec_DMD,
                                                lam_POD,
                                                eig_vec_POD,
                                                r );

            if ( Phi.cols() == 0 )  break;
            if ( count == 0)    svd_old = lam_POD.cwiseProduct(lam_POD);

            svd_new = lam_POD.cwiseProduct(lam_POD);
            eps = (svd_old.sum() - svd_new.sum())/svd_old.sum();
            
            if ( i > 0 )    K_pc(i-1) = eps;
            
            count ++;
            std::cout << "Energy at Iteration " << i << " : " << std::setprecision(12) << eps*100 << "%" << std::endl;
            
            if ( i == rdmd )    break;

            Eigen::MatrixXd Phi_r = Phi.real();
            Eigen::MatrixXd coef_mod(Phi.cols(),Ns);
        // Real part Modes normalization
            for ( int j = 0; j < Phi.cols(); j++ ) {
                double sum = 0.0;
                for ( int k = 0; k < Np; k++ )
                    sum += Phi_r(k,j)*Phi_r(k,j);

                if ( std::abs(sum) < 1e-15 ) {
                    Phi_r.col(j) = Eigen::MatrixXd::Zero(Np,1);
                } else {
                    Phi_r.col(j) = Phi_r.col(j)/std::sqrt(sum);
                }
            }
            Eigen::VectorXd residual_average(Phi.cols());
            int min_idx;

            for ( int r_dmd = 0; r_dmd < Phi.cols(); r_dmd++ ) {
                for ( int nt = 0; nt < Ns; nt++ ) {
                    coef_mod(r_dmd, nt) = res_set.col(nt).transpose()*Phi_r.col(r_dmd);
                    residual_time = res_set.col(nt) - coef_mod(r_dmd, nt)*Phi_r.col(r_dmd);
                    residual_time_norm(nt) = residual_time.norm();
                }

                double mean = 0.0;
                for ( int m = 0; m < Ns; m++ )
                    mean += residual_time_norm(m);

                residual_average(r_dmd) = mean/Ns;

            }

            double min_Val = residual_average.minCoeff( &min_idx );
            Phi_RDMD.col(i) = Phi_r.col(min_idx);
            Coefs.row(i) = coef_mod.row(min_idx);
            lambda(i) = lam_DMD(min_idx).real();
            res_set = res_set - Phi_RDMD.col(i)*Coefs.row(i);
        }
    }

    return Phi_RDMD;

}


Eigen::MatrixXd RDMD_lsq_basis ( const Eigen::MatrixXd &sn_set,
                                   Eigen::MatrixXd &Coefs,
                                   Eigen::VectorXd &lambda,
                                   Eigen::VectorXd &K_pc,
                                   const int r,
                                   int &rdmd,
                                   double En )
{
    int Np = sn_set.rows();
    int Ns = sn_set.cols();

    Eigen::MatrixXd Phi_RDMD = Eigen::MatrixXd::Zero(Np, Ns);
    Eigen::MatrixXd res_set = sn_set;
    Eigen::VectorXd lam_POD;
    Eigen::VectorXcd lam_DMD;
    Eigen::MatrixXd eig_vec_POD;
    Eigen::MatrixXcd eig_vec_DMD;

    Eigen::VectorXd residual_time(Np);
    Eigen::VectorXd residual_time_norm(Ns);

    double eps = 0.0;
    Eigen::VectorXd svd_new = Eigen::VectorXd::Zero(Ns);
    Eigen::VectorXd svd_old = Eigen::VectorXd::Zero(Ns);
    double count = 0;

    if ( rdmd > Ns ) {
        std::cout << "Rank RDMD too high for the number of snapshots available. Resetting it to maximum value admissable" << std::endl;
        rdmd = Ns;
    }

    if ( rdmd == 0 ) {
        while ( eps < En && count < Ns+1 ) {
            //Perform pure DMD
            Eigen::MatrixXcd Phi = DMD_basis( res_set,
                                              lam_DMD,
                                              eig_vec_DMD,
                                              lam_POD,
                                              eig_vec_POD,
                                              r );

            if ( Phi.cols() == 0 )
                break;

            if ( count == 0)
                svd_old = lam_POD.cwiseProduct(lam_POD);

            svd_new = lam_POD.cwiseProduct(lam_POD);
            eps = (svd_old.sum() - svd_new.sum())/svd_old.sum();

            Eigen::MatrixXd Phi_r = Phi.real();
            Eigen::MatrixXd coef_mod(Phi.cols(),Ns);

            for ( int j = 0; j < Phi.cols(); j++ ) {
                double sum = 0.0;
                for ( int k = 0; k < Np; k++ )
                    sum += Phi_r(k,j)*Phi_r(k,j);

                Phi_r.col(j) = Phi_r.col(j)/std::sqrt(sum);
            }

            Eigen::MatrixXd PhiTU = Phi_r.transpose()*res_set;
            Eigen::MatrixXd PhiTPhi = Phi_r.transpose()*Phi_r;
            coef_mod = PhiTPhi.colPivHouseholderQr().solve(PhiTU);

            Eigen::VectorXd residual_average(Phi.cols());

            int min_idx;

            for ( int r_dmd = 0; r_dmd < Phi.cols(); r_dmd++ ) {
                for ( int nt = 0; nt < Ns; nt++ ) {
                    residual_time = res_set.col(nt) - coef_mod(r_dmd, nt)*Phi_r.col(r_dmd);
                    residual_time_norm(nt) = residual_time.norm();
                }
                double mean = 0.0;
                for ( int m = 0; m < Ns; m++ )
                    mean += residual_time_norm(m);

                residual_average(r_dmd) = mean/Ns;
            }

            double min_Val = residual_average.minCoeff( &min_idx );

            Phi_RDMD.col(count) = Phi_r.col(min_idx);
            Coefs.row(count) = coef_mod.row(min_idx);
            lambda(count) = lam_DMD(min_idx).real();
            res_set = res_set - Phi_RDMD.col(count)*Coefs.row(count);
            std::cout << "Energy level at iteration " << count << " : " << std::setprecision(12) << eps*100 << "%" << std::endl;

            count ++;
        }
        rdmd = count;

    } else {
        for ( int i = 0; i <= rdmd; i++ ){ //Considering also i = rdmd only for computing energy at the last iteration
            //Perform pure DMD
            Eigen::MatrixXcd Phi = DMD_basis( res_set,
                                              lam_DMD,
                                              eig_vec_DMD,
                                              lam_POD,
                                              eig_vec_POD,
                                              r );

            if ( Phi.cols() == 0 )  break;
            if ( count == 0)    svd_old = lam_POD.cwiseProduct(lam_POD);

            svd_new = lam_POD.cwiseProduct(lam_POD);
            eps = (svd_old.sum() - svd_new.sum())/svd_old.sum();

            if ( i > 0 )    K_pc(i-1) = eps;

            count ++;
            std::cout << "Energy at Iteration " << i << " : " << std::setprecision(12) << eps*100 << "%" << std::endl;

            if ( i == rdmd )    break;

            Eigen::MatrixXd Phi_r = Phi.real();
            Eigen::MatrixXd coef_mod(Phi.cols(),Ns);
            // Real part Modes normalization
            for ( int j = 0; j < Phi.cols(); j++ ) {
                double sum = 0.0;
                for ( int k = 0; k < Np; k++ )
                    sum += Phi_r(k,j)*Phi_r(k,j);

                Phi_r.col(j) = Phi_r.col(j)/std::sqrt(sum);
            }

            Eigen::MatrixXd PhiTU = Phi_r.transpose()*res_set;
            Eigen::MatrixXd PhiTPhi = Phi_r.transpose()*Phi_r;
            coef_mod = PhiTPhi.colPivHouseholderQr().solve(PhiTU);

            Eigen::VectorXd residual_average(Phi.cols());
            int min_idx;

            for ( int r_dmd = 0; r_dmd < Phi.cols(); r_dmd++ ) {
                for ( int nt = 0; nt < Ns; nt++ ) {
                    residual_time = res_set.col(nt) - coef_mod(r_dmd, nt)*Phi_r.col(r_dmd);
                    residual_time_norm(nt) = residual_time.norm();
                }

                double mean = 0.0;
                for ( int m = 0; m < Ns; m++ )
                    mean += residual_time_norm(m);

                residual_average(r_dmd) = mean/Ns;

            }

            double min_Val = residual_average.minCoeff( &min_idx );

            Phi_RDMD.col(i) = Phi_r.col(min_idx);
            Coefs.row(i) = coef_mod.row(min_idx);
            lambda(i) = lam_DMD(min_idx).real();
            res_set = res_set - Phi_RDMD.col(i)*Coefs.row(i);

        }
    }

    return Phi_RDMD;

}









//Eigen::MatrixXcd fbDMD_basis ( const Eigen::MatrixXd &snap_set,
//                            Eigen::VectorXcd &lam,
//                            Eigen::MatrixXcd &eig_vec,
//                            const int r )
//{
//
//    int Ns = snap_set.cols() - 1;
//    int Nm;
//
//    Eigen::BDCSVD<Eigen::MatrixXd> svdX( snap_set.leftCols(Ns),
//                                        Eigen::ComputeThinU | Eigen::ComputeThinV );
//    Eigen::VectorXd lam_PODX = svdX.singularValues();
//    Eigen::MatrixXd eig_vec_PODX = svdX.matrixV();
//    eig_sort(lam_PODX, eig_vec_PODX);
//    Eigen::MatrixXd Nullspace = nullspace( lam_PODX, eig_vec_PODX );
//    bool test = check_linear_consistency( snap_set.leftCols(Ns), snap_set.rightCols(Ns), Nullspace );
//
//
//    Eigen::MatrixXd UX = svdX.matrixU();
//    Eigen::MatrixXd Sig_invX = Eigen::MatrixXd::Zero(UX.cols(), UX.cols());
//
//    for ( int i = 0; i < UX.cols(); i++ )
//        Sig_invX(i, i) = 1.0/lam_PODX(i);
//
//
//    Eigen::BDCSVD<Eigen::MatrixXd> svdY( snap_set.rightCols(Ns),
//                                        Eigen::ComputeThinU | Eigen::ComputeThinV );
//    Eigen::VectorXd lam_PODY = svdY.singularValues();
//    Eigen::MatrixXd eig_vec_PODY = svdY.matrixV();
//    eig_sort(lam_PODY, eig_vec_PODY);
//
//
//    Eigen::MatrixXd UY = svdY.matrixU();
//    Eigen::MatrixXd Sig_invY = Eigen::MatrixXd::Zero(UY.cols(), UY.cols());
//
//    for ( int i = 0; i < UX.cols(); i++ )
//        Sig_invY(i, i) = 1.0/lam_PODY(i);
//
//    if ( r == 0)
//    {
//        Nm = SVHT ( lam_PODX, Ns, snap_set.rows() );
//        std::cout << "DMD-rank from SVHT : " << Nm << std::endl;
//    }
//    else
//    {
//        Nm = std::min(r, Ns);
//        std::cout << "DMD user-defined rank : " << Nm << std::endl;
//    }
//    std::cout << "Singular Values : \n" << lam_PODX.head(Nm) << std::endl;
//    Eigen::MatrixXd fAtilde = UX.leftCols(Nm).transpose()*snap_set.rightCols(Ns)*
//                                eig_vec_PODX.leftCols(Nm)*Sig_invX.block(0,0,Nm,Nm);
//
//    Eigen::MatrixXd bAtilde = UY.leftCols(Nm).transpose()*snap_set.leftCols(Ns)*
//                            eig_vec_PODY.leftCols(Nm)*Sig_invY.block(0,0,Nm,Nm);
//
//    Eigen::MatrixXd Atilde_sq = fAtilde*bAtilde.inverse();
//    std::cout << "Determinant bAtilde : " << bAtilde.determinant() << std::endl;
//
//    Eigen::MatrixXd Atilde = Atilde_sq.sqrt();
//    // std::cout << "Atilde :\n " << Atilde << std::endl;
//
//    Eigen::EigenSolver<Eigen::MatrixXd> es(Atilde);
//    lam = es.eigenvalues();
//    eig_vec = es.eigenvectors();
//
//    Eigen::MatrixXcd appo = snap_set.rightCols(Ns)*
//                            eig_vec_PODX.leftCols(Nm)*Sig_invX.block(0,0,Nm,Nm);
//
//    // Divido per lambda
//
//    Eigen::MatrixXcd phi(snap_set.rows(), Nm);
//    for (int i = 0; i < Nm; i++)
//        phi.col(i) = 1.0/lam(i)*appo*eig_vec.col(i);
//
//    return phi;
//
//    //Standard DMD
//    // return U*eig_vec;
//
//}


Eigen::MatrixXcd HODMD_basis( const Eigen::MatrixXd &snap_set,
                            Eigen::VectorXcd &lambda,
                            Eigen::MatrixXcd &eig_vec,
                            Eigen::VectorXcd &Coefs,
                            const double tol,
                            const int d )
{

    int Ns = snap_set.cols();
    Eigen::BDCSVD<Eigen::MatrixXd> svd( snap_set, 
                                        Eigen::ComputeThinU | Eigen::ComputeThinV );
    Eigen::VectorXd s = svd.singularValues();

    int N_red1 = 0;
    double eps1 = 0.0;
     
    while ( eps1 < tol && N_red1 < s.size() )
    {
        eps1 = s.tail(N_red1 + 1).sum()/s.sum();
        N_red1++;
     
    }

    int N_svd = s.size() - (N_red1 - 1);

    std::cout << std::endl;
    std::cout << " First reduction size : " << N_red1 - 1 << " out of " << Ns << std::endl; 

    Eigen::MatrixXd T = svd.matrixV();
    Eigen::MatrixXd U = svd.matrixU();
    Eigen::MatrixXd Sigma = Eigen::MatrixXd::Zero(N_svd,N_svd);

    for ( int i = 0; i < N_svd; i++ )
        Sigma(i,i) = s(i);

    Eigen::MatrixXd sn_set_hat = Sigma*T.leftCols(N_svd).transpose();
    Eigen::MatrixXd sn_set_hat_ho(d*N_svd, Ns - d);

    for ( int i = 0; i < (Ns - d); i++ )
    {
        for( int j = 0; j < d; j++ )
            sn_set_hat_ho.col(i).segment(j*N_svd, N_svd) = sn_set_hat.col(j+i);
    } 

    Eigen::BDCSVD<Eigen::MatrixXd> svd2( sn_set_hat_ho.leftCols(sn_set_hat_ho.cols()-1), 
                                        Eigen::ComputeThinU | Eigen::ComputeThinV );
    Eigen::VectorXd s2 = svd2.singularValues();
    Eigen::MatrixXd T_hat_ho = svd2.matrixV();
    Eigen::MatrixXd U_hat_ho = svd2.matrixU();

    eps1 = 0.0;
    N_red1 = 0;

    while ( eps1 < tol && N_red1 < s2.size() )
    {
        eps1 = s2.tail(N_red1 + 1).sum()/s2.sum();
        N_red1++;
     
    }

    std::cout << " Second reduction size : " << N_red1 - 1 << " out of " << s2.size() << std::endl;
    int N_svd2 = s2.size() - (N_red1 - 1);
    Eigen::MatrixXd Sigma2 = Eigen::MatrixXd::Zero(N_svd2, N_svd2);

    for ( int i = 0; i < N_svd2; i++ )
        Sigma2(i,i) = s2(i);

    Eigen::MatrixXd R_hat = sn_set_hat_ho.rightCols(sn_set_hat_ho.cols()-1)*
                            T_hat_ho.leftCols(N_svd2)*Sigma2.inverse()*U_hat_ho.leftCols(N_svd2).transpose();
    Eigen::EigenSolver<Eigen::MatrixXd> es(R_hat); 
    lambda = es.eigenvalues(); 
    Eigen::MatrixXcd Full_eig_vec = es.eigenvectors();
    eig_vec = Full_eig_vec.block(0,0,N_svd2,Full_eig_vec.cols());

    Coefs = Calculate_Coefs_DMD_exact ( sn_set_hat.leftCols(sn_set_hat.cols() - 1),  
                                        lambda,  
                                        eig_vec );

    return U.leftCols(N_svd2)*eig_vec;

}


Eigen::MatrixXd GPOD_basis( const double Dt,
                                const Eigen::MatrixXd &snap_set,
                                Eigen::VectorXd &lam,
                                Eigen::MatrixXd &Coeffs,
                                const int &r)
{
    int Ns = snap_set.cols();
    int Nr = snap_set.rows();
    Eigen::MatrixXd Gradients_T(Nr, Ns-2);

    for ( int i = 2; i < Ns; i++ )
        Gradients_T.col(i-2) = 0.5*(3.0*snap_set.col(i) - 4.0*snap_set.col(i-1) + snap_set.col(i-2))/Dt;


    Eigen::BDCSVD<Eigen::MatrixXd> svd( Gradients_T,
                                        Eigen::ComputeThinU | Eigen::ComputeThinV );

    lam = svd.singularValues();
    Eigen::MatrixXd eig_vec = svd.matrixV();
    eig_sort(lam, eig_vec);
    Eigen::MatrixXd U = svd.matrixU();
    int Nm;

    if ( r == 0) {
        Nm = SVHT ( lam, Ns-2, Nr );
        std::cout << "GPOD-rank from SVHT : " << Nm << std::endl;
    } else if ( r > 0 ) {
        Nm = std::min(r, Not_zero ( lam ));
        std::cout << "GPOD user-defined rank : " << Nm << std::endl;
    } else {
        Nm = Not_zero ( lam );
        std::cout << "GPOD rank based on non zero singular values : " << Nm << std::endl;
    }

    Coeffs = U.transpose()*snap_set;
    return U.leftCols(Nm);

}



Eigen::MatrixXcd DMD_Adaptive_basis ( const Eigen::MatrixXd &snap_set,
                             Eigen::VectorXcd &lam,
                             Eigen::MatrixXcd &eig_vec,
                             Eigen::VectorXd &lam_POD,
                             Eigen::MatrixXd &eig_vec_POD,
                             Eigen::VectorXi &tpos )
{

    int Nsamp = tpos.size();
    int Np = snap_set.rows();
    int Nm;
    //                                Eigen::MatrixXd U = SPOD_basis(snap_set.leftCols(Ns), lam_POD, K_pc, eig_vec_POD );
    Eigen::VectorXi tpos_shift = tpos + Eigen::VectorXi::Ones(Nsamp);
    Eigen::MatrixXd sub_sn_set = indexing(snap_set, Eigen::ArrayXi::LinSpaced(Np,0,Np-1),tpos);
    Eigen::MatrixXd sub_sn_set_shift = indexing(snap_set, Eigen::ArrayXi::LinSpaced(Np,0,Np-1),tpos_shift);

    Eigen::BDCSVD<Eigen::MatrixXd> svd( sub_sn_set,Eigen::ComputeThinU | Eigen::ComputeThinV );
    lam_POD = svd.singularValues();
    int N_notZero = Not_zero(lam_POD);

    Nm = std::min(N_notZero,Nsamp);
    eig_vec_POD = svd.matrixV();
    eig_sort(lam_POD, eig_vec_POD);

    Eigen::MatrixXd U = svd.matrixU();
    Eigen::MatrixXd Sig_inv = Eigen::MatrixXd::Zero(Nm, Nm);

    for ( int i = 0; i < Nm; i++ )
        Sig_inv(i, i) = 1.0/lam_POD(i);

    Eigen::MatrixXcd phi = Eigen::MatrixXd::Zero(Np, Nm);

    Eigen::MatrixXd Atilde = U.leftCols(Nm).transpose()*sub_sn_set_shift*
                             eig_vec_POD.leftCols(Nm)*Sig_inv;

    if ( Atilde.size() == 1 && Atilde(0,0) == 0.0 ) {
        lam = Eigen::VectorXcd::Zero(1);
        eig_vec = Eigen::MatrixXcd::Ones(1,1);
    } else {
        Eigen::EigenSolver<Eigen::MatrixXd> es(Atilde);
        // Full_lam = es.eigenvalues();
        lam = es.eigenvalues();
        eig_vec = es.eigenvectors();
    }

    Eigen::MatrixXcd appo = sub_sn_set_shift* eig_vec_POD.leftCols(Nm) * Sig_inv;

    //Final computation and normalization of modes
    for (int i = 0; i < Nm; i++) {
        phi.col(i) = 1.0/lam(i)*appo*eig_vec.col(i);
        phi.col(i) = phi.col(i)/phi.col(i).norm();
    }

    return phi;

    //Standard DMD
    // return U*eig_vec;
}


Eigen::MatrixXd RDMD_Adaptive_basis ( const Eigen::MatrixXd &sn_set,
                                      Eigen::MatrixXd &Coefs,     //Define N_mod RDMD through the dimension of matrix Coefs
                                      Eigen::VectorXd &K_pc,
                                      Eigen::VectorXi &tpos)
{
    int Np = sn_set.rows();
    int Ns = sn_set.cols();
    int rdmd = tpos.size();

    Eigen::MatrixXd Phi_RDMD = Eigen::MatrixXd::Zero(Np, Ns);
    Eigen::MatrixXd res_set = sn_set;
    Eigen::VectorXd lam_POD;
    Eigen::VectorXcd lam_DMD;
    Eigen::MatrixXd eig_vec_POD;
    Eigen::MatrixXcd eig_vec_DMD;

    Eigen::VectorXd residual_time(Np);
    Eigen::VectorXd residual_time_norm(Ns);

    double eps = 0.0;
    Eigen::VectorXd svd_new = Eigen::VectorXd::Zero(Ns);
    Eigen::VectorXd svd_old = Eigen::VectorXd::Zero(Ns);
    double count = 0;

    if ( rdmd > Ns ) {
        std::cout << "Number of adaptive samples greater than available samples\n Exiting ..." << std::endl;
        exit(EXIT_FAILURE);
    }

    for ( int i = 0; i <= rdmd; i++ ) { //Considering also i = rdmd only for computing energy at the last iteration
        //Perform pure DMD with equal time shifts
//        Eigen::MatrixXcd Phi = DMD_Adaptive_basis ( res_set,
//                                                    lam_DMD,
//                                                    eig_vec_DMD,
//                                                    lam_POD,
//                                                    eig_vec_POD,
//                                                    tpos );

        //Perform pure DMD with non equal time shifts
        Eigen::VectorXi Ipos = Inverse_POS(res_set, Ns);
        Eigen::MatrixXd sub_res_set = indexing(res_set, Eigen::ArrayXi::LinSpaced(Np,0,Np-1),Ipos);
        Eigen::MatrixXcd Phi = DMD_basis(sub_res_set,
                                         lam_DMD,
                                         eig_vec_DMD,
                                         lam_POD,
                                         eig_vec_POD,
                                         -1);

        if ( Phi.cols() == 0 )  break;
        if ( count == 0)    svd_old = lam_POD.cwiseProduct(lam_POD);
        svd_new = lam_POD.cwiseProduct(lam_POD);
        eps = (svd_old.sum() - svd_new.sum())/svd_old.sum();
        if ( i > 0 )    K_pc(i-1) = eps;

        count ++;
        std::cout << "Energy at Iteration " << i << " : " << std::setprecision(12) << eps*100 << "%" << std::endl;
        std::cout << "DMD-rank at current iteration :" << Phi.cols() << std::endl;

        if ( i == rdmd )    break;

        Eigen::MatrixXd Phi_r = Phi.real();
        Eigen::MatrixXd coef_mod(Phi.cols(),Ns);
        // Real part Modes normalization
        for ( int j = 0; j < Phi.cols(); j++ ) {
            double sum = 0.0;
            for ( int k = 0; k < Np; k++ )
                sum += Phi_r(k,j)*Phi_r(k,j);

            if ( std::abs(sum) < 1e-15 ) {
                Phi_r.col(j) = Eigen::MatrixXd::Zero(Np,1);
            } else {
                Phi_r.col(j) = Phi_r.col(j)/std::sqrt(sum);
            }
        }

        Eigen::VectorXd residual_average(Phi.cols());
        int min_idx;

        for ( int r_dmd = 0; r_dmd < Phi.cols(); r_dmd++ ) {
            for ( int nt = 0; nt < Ns; nt++ ) {
                coef_mod(r_dmd, nt) = res_set.col(nt).transpose()*Phi_r.col(r_dmd);
                residual_time = res_set.col(nt) - coef_mod(r_dmd, nt)*Phi_r.col(r_dmd);
                residual_time_norm(nt) = residual_time.norm();
            }
            double mean = 0.0;
            for ( int m = 0; m < Ns; m++ )
                mean += residual_time_norm(m);

            residual_average(r_dmd) = mean/Ns;
        }
        double min_Val = residual_average.minCoeff( &min_idx );
        Phi_RDMD.col(i) = Phi_r.col(min_idx);
        Coefs.row(i) = coef_mod.row(min_idx);
        res_set = res_set - Phi_RDMD.col(i)*Coefs.row(i);
    }

    return Phi_RDMD;
}
