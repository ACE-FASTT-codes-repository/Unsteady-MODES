//
// Created by haitan on 01/06/2020.
//
#include "Post-Process.hpp"

void Write_History_ResError(prob_settings settings, std::string method, int it1, int it2, bool d_err){

    //Defining Common Variables
    std::string header_history, header_error, value_history, value_error; //to write headers in global files
    std::ifstream history_su2;
    std::ifstream error_su2;
    std::ofstream history_global;
    std::ofstream error_global;

    //Getting filenames
    int iter = std::round(settings.t_res[it2]/settings.Dt_cfd);
    std::stringstream buffer;
    buffer << std::setfill('0') << std::setw(5) << std::to_string(iter);
    std::string filename_history_su2 = "history_rbm_" + buffer.str() + ".csv";
    std::string filename_error_su2 = "error_rbm.csv";
    std::string filename_history_global = "history_" + method + "_" + std::to_string(settings.r)
                                            + "_" + std::to_string(settings.Dt_res[it1]) + ".csv";
    std::string filename_error_global = "Error_" + method + "_" + std::to_string(settings.r) + ".csv";


    /*----------------------First Reading outputs from SU2-------------------*/
    history_su2.open(filename_history_su2);
    error_su2.open(filename_error_su2);

    if ( d_err ){
        if ( error_su2.is_open() ){
            std::string linedata2, token;

            //Getting row of headers
            getline(error_su2,linedata2);
            header_error = linedata2;

            //Getting values
            getline(error_su2,linedata2);
            value_error = linedata2;

        } else {
            std::cout << "Unable to open SU2 error files. Exiting ..." << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    if ( history_su2.is_open() ){
        std::string linedata1, token;

        //Getting row of headers
        getline(history_su2,linedata1);
        header_history = linedata1;

        //Getting values
        getline(history_su2,linedata1);
        value_history = linedata1;

    } else {
        std::cout << "Unable to open SU2 history files. Exiting ..." << std::endl;
        exit(EXIT_FAILURE);
    }


    /*-----------------Then appending to history and error global files-----------------------*/

    if ( it2 == 0 ){
        history_global.open(filename_history_global);

        if ( history_global.is_open() ) {
            header_history = "\"T(s)\"," + header_history;
            history_global << header_history << std::endl;
            history_global << settings.t_res[it2] << "," << value_history << std::endl;

            history_global.close();

        } else {
            std::cout << "Unable to open history files for writing. Exiting... " << std::endl;
            exit(EXIT_FAILURE);
        }

        if ( d_err ) {
            error_global.open(filename_error_global);

            if ( error_global.is_open() ) {
                header_error = "\"T(s)\"," + header_error;
                error_global << header_error << std::endl;
                error_global << settings.t_res[it2] << "," << value_error << std::endl;

                error_global.close();
            } else {
                std::cout << "Unable to open error files for writing. Exiting... " << std::endl;
                exit(EXIT_FAILURE);
            }
        }
    } else {
        history_global.open(filename_history_global, std::ofstream::out | std::ofstream::app);

        if ( history_global.is_open()  ) {

            history_global << settings.t_res[it2] << "," << value_history << std::endl;
            history_global.close();

        } else {
            std::cout << "Unable to open history and error files for writing. Exiting... " << std::endl;
            exit(EXIT_FAILURE);
        }

        if ( d_err ){
            error_global.open(filename_error_global, std::ofstream::out | std::ofstream::app);

            if ( error_global.is_open() ) {

                error_global << settings.t_res[it2] << "," << value_error << std::endl;
                error_global.close();

            } else {
                std::cout << "Unable to open error files for writing. Exiting... " << std::endl;
                exit(EXIT_FAILURE);
            }
        }
    }

    //Removing files generated from SU2_DTR no longer needed
    auto opt = system("rm -f history_rbm* error_rbm.csv forces_breakdown.dat surface_flow*");

}