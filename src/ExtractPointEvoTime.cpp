/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------ Copyright (C) 2021 University of Strathclyde and Authors ------
-------------------- e-mail: gaetano.pascarella.ac.uk ----------------
----------------------- Author: Gaetano Pascarella -----------------------
*/
/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------ Copyright (C) 2021 University of Strathclyde and Authors ------
-------------------- e-mail: gaetano.pascarella.ac.uk ----------------
----------------------- Author: Gaetano Pascarella -----------------------

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
    std::cout << "Extracting Single Point feature start" << std::endl;

    prob_settings settings;
    std::string filecfg = argv[1];

    Read_cfg( filecfg, settings );

    int iD;
    int Ns = settings.Ns;
    int ds = settings.Ds;
    int init = settings.nstart;
    std::vector<int> Cols = settings.Cols;

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
    Eigen::MatrixXd Coords = read_colnew( file_1, Nr, settings.Cols_coords );
    std::cout << "Done " << std::endl;
    Eigen::VectorXd distSource(Nr);

    std::cout << "Point source :" << std::endl;
    for ( auto x : settings.pointsource ) std::cout << x << "\t";
    std::cout << std::endl;
    Eigen::RowVectorXd psource(settings.pointsource.size());
    for ( int i = 0; i< settings.pointsource.size();i++ ) psource[i] = settings.pointsource[i];

    for (int iPoint = 0; iPoint<Nr; iPoint++ ) distSource(iPoint) = (Coords.row(iPoint)-psource).norm();
    double min_Val = distSource.minCoeff( &iD );
    std::cout << "ID of minimum " << iD << std::endl;

    Eigen::MatrixXd field(Nr,Cols.size());
    Eigen::MatrixXd EvoQuantities(settings.Ns,Cols.size());

    std::string file_temp;
    int count = 0;
    for( int i = init; i < (Ns*ds + init); i += ds ){

        std::stringstream buffer;
        buffer << std::setfill('0') << std::setw(5) << std::to_string(i);
        file_temp = root_inputfile + "_" + buffer.str() + "." + input_format;
        std::cout << "Reading fields from : " << file_temp << "\t";
        field = read_colnew(file_temp, Nr, Cols);
        std::cout << "Complete!" << std::endl;
        for (int iVar = 0; iVar<Cols.size(); iVar++)
            EvoQuantities(count,iVar) = field(iD,iVar);
        count++;
    }

    //Writing to file
    std::string filename = "ConservativeEvo_ID.dat";
    std::ofstream flow_data;
    flow_data.open(filename.c_str());

    // Write row of Headers
    std::string coef;
    for ( int i = 0; i < Cols.size(); i++ ) {
        coef = "\"Var" + std::to_string(i+1) + "\"";
        flow_data << coef << " ";
    }
    flow_data << std::endl;
    //Write fields
    for ( int i = 0; i < Ns; i++ ) {
        for (int j = 0; j < Cols.size(); j++)
            flow_data << std::setprecision(12) << std::scientific << EvoQuantities(i,j) <<  " ";

        flow_data << std::endl;
    }
    // Close file
    flow_data.close();

    std::cout << "Extracting Single Point Feature ends " << std::endl;

    return 0;
}
