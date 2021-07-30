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
    Read_cfg( filecfg, settings );

    // Pre-processing common vars
    int nC = settings.Cols.size();
    int Nr;
    Eigen::MatrixXd Coords;
    std::vector<double> t_vec(settings.Ns);
    common_vars( Nr, Coords, t_vec, settings);

    std::cout << "Storing snapshot Matrix ... \n ";
    Eigen::MatrixXd sn_set = generate_snap_matrix( Nr, settings);
    std::cout << std::endl;

    std::cout << "Computing mean/Initial Condition of CFD solution ... " << std::endl;
    Eigen::VectorXd Ic = Eigen::VectorXd::Zero(nC*Nr);

    std::string binary = "YES";
    Eigen::VectorXd svd_cum_sum(settings.Ns);

    if ( settings.flag_mean == "IC" ) {
        std::cout << "Subtracting Initial condition" << std::endl << std::endl;
        Ic = IC(sn_set,settings,nC,Nr);
    } else {
        std::cout << "Using data without subtracting any reference state" << std::endl << std::endl;
    }



    std::cout << "--------------------------------------" << std::endl ;
    std::cout << "-------Performin RDMD ResEval---------" << std::endl ;
    std::cout << "--------------------------------------" << std::endl ;
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

        COEF[ncons] = Coefs.transpose();
        surr_coefs[ncons] = getSurrCoefs(t_vec,
                                         Coefs.transpose(),
                                         settings.flag_interp);

        N_notZero = Phi[ncons].cols();
        if (settings.r == 0) Nm = Nmod(settings.En, K_pc);
        else Nm = std::min(settings.r, N_notZero);
        std::cout << "Number of modes used in reconstruction " << Nm << std::endl;

        }

    ModesDB_Write( Phi, COEF, settings);

    Eigen::VectorXd t_interp = Eigen::VectorXd::LinSpaced(t_vec.size()*100,0.,t_vec[t_vec.size()-1]);


    for (int ncons = 0; ncons < nC; ncons++) {

        Eigen::MatrixXd coef_t(t_interp.size(), Nm);
        std::vector<double> tr(1);

        for (int j = 0; j < t_interp.size(); j++) {
            tr[0] = t_interp[j];
            for (int i = 0; i < Nm; i++)
                surr_coefs[ncons][i].evaluate(tr, coef_t(j, i));
        }

        std::ofstream coef_data;
        std::string filetemp2 = "CoefsInterpConsVar" + std::to_string(ncons+1) + ".dat";
        coef_data.open(filetemp2.c_str());
        //Writting row of headers
        for ( int iMode = 0; iMode < Nm; iMode++ )
            coef_data << "\"CoefMode_" + std::to_string(iMode + 1) + "\"";
        coef_data << std::endl;
        //Writing Coefs field
        for ( int iTime = 0; iTime < t_interp.size(); iTime++ ) {
            for (int iMode = 0; iMode < Nm; iMode++){
                coef_data << std::setprecision(12) << std::scientific << coef_t(iTime,iMode) <<  " ";
            }
            coef_data << std::endl;
        }
        // Close file
        coef_data.close();
    }



    std::cout << "MODES Adaptive Residual Evaluation ends" << std::endl << std::endl;
    return 0;

}
