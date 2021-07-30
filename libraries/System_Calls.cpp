#include "System_Calls.hpp"

void SU2_DTR(prob_settings settings, std::string su2_conf,  std::string method, int it1, int it2 ){

    bool direct_error = settings.direct_error;
    bool surf_res = settings.surf_res;
    /*------------Defining all necessary strings----------------*/

    //New cfg file SU2 for residual evaluation and modification for SU2_DTR
    std::string root_conf;
    root_conf.assign ( su2_conf, 0, su2_conf.size() - 4);
    std::string su2_conf_new = root_conf + "-reseval.cfg";

        //If gust is active we need free stream velocity
    double M = settings.Mach;
    double T = settings.T;
    double alpha = M_PI*settings.alpha/double(180);
    double R = 287.058;
    double gamma = 1.4;

    double V_inf = M*std::sqrt(gamma*R*T)*std::cos(alpha);

    //Check if you can compute also direct error
    int iter = std::round(settings.t_res[it2]/settings.Dt_cfd);
    double temp;
    if( std::abs((double(iter) - settings.t_res[it2]/settings.Dt_cfd) > 1e-5) && settings.direct_error ) {
        std::cout << "lack = " << std::abs(double(iter) - settings.t_res[it2]/settings.Dt_cfd) << " WARNING! Not computing direct error for this time instant" << std::endl;
        direct_error = false;
    }

    Modify_su2_cfg ( su2_conf, su2_conf_new, settings, it1, it2, V_inf );

    //String to launch SU2_DTR from terminal
    std::string su2dtr_string = "mpirun -np 6 ./SU2_DTR " + su2_conf_new + " > SU2.log"; // + " > resEval_su2.log";
    int len_s = su2dtr_string.length();
    char su2_sys_call[len_s + 1];
    strcpy(su2_sys_call, su2dtr_string.c_str());

    //String for removing useless solutions
    std::string root_outputfile;
    root_outputfile.assign ( settings.out_file, 0, settings.out_file.size() - 4);
    std::string root_inputfile;
    std::string rmf_string = "rm -f " + root_outputfile + "_*";
    len_s = rmf_string.length();
    char rmf_sys_call[len_s + 20];
    strcpy(rmf_sys_call, rmf_string.c_str());


    //String for saving reconstructions when flag rec is set to yes
    std::stringstream buffer;
    buffer << std::setfill('0') << std::setw(5) << std::to_string(iter);
    std::string mv_resrec_string = "mv " + root_outputfile + "_" + buffer.str() + ".dat  RecRes" + method + "_" + buffer.str() + ".dat";
    len_s = mv_resrec_string.length();
    char mv_resrec_sys_call[len_s + 1];
    strcpy(mv_resrec_sys_call, mv_resrec_string.c_str());


    //Performing system calls: SU2 for residual and direct error evaluation, saving reconstructions, deleting useless solutions
    std::cout << "Calling SU2 for residual evaluation and writing file to history " << std::endl;
    auto opt = std::system(su2_sys_call);
    if ( settings.flag_rec == "YES" ) {
        std::cout << "Saving Reconstructions " << std::endl;
        opt = std::system(mv_resrec_sys_call);
    }
    opt = std::system(rmf_sys_call);

    Write_History_ResError(settings, method, it1, it2, direct_error);

    std::cout << std::endl;
}