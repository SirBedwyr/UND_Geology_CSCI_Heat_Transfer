/**
 * Performs a finite difference heat flow 
 * simulation using conduction and convection.
 */

#ifdef _WIN32
#define _USE_MATH_DEFINES
#define NOMINMAX //FYI need to disable min/max macro in windows.h
#include <windows.h>
#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <limits>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * Time_step (seconds/yr) divided by the product of density and heat capacity.
 * The value for the density is 2200 kg/m^3 and heat capacity is 1000 kJ/kg K.
 */
#define QFAC 14.33              //Description is defined in the previous comment
#define DTC 0.25                //
#define OUT_PRECISION 10        //Number of digits to print after the decimal place for floating point values
#define INDEX_WIDTH 2           //The number of characters to print and read for each conduction and convection code
#define REAL double             //The precision of the model.

using std::cerr;
using std::cin;
using std::cout;
using std::endl;
using std::string;
using std::ofstream;
using std::ifstream;
using std::ostringstream;
using std::setw;
using std::right;
using std::left;
using std::fixed;
using std::scientific;
using std::setprecision;
using std::setfill;
using std::ios;
using std::numeric_limits;
using std::streamsize;
using std::max;
using std::flush;

void save_surfer();
void save_model_state();
void conduction();
void convection();
void PressEnterToContinue();
REAL find_max_temp_diff();
void update_moving_sources();
void find_loc_index(REAL x_loc, REAL y_loc, REAL z_loc, int *index);

//Conduction code specific variables
int ***cond_codes;              //The unmodified conduction codes as read from the input file

//Convection code specific variables
int ***conv_codes;              //The convection codes as read from the input file

//File names
string source_filename;         //The input files name with extension
string output_filename;         //The output state files name with extension

//Input file variables
string title;                   //The title of the input file
int using_convection = -1;      //Indicates if convection is being used
REAL ***temp;                   //The current temperature array
int num_rows;                   //The number of rows for the simulation
int num_cols;                   //The number of columns for the simulation
int num_slices;                 //Total number of slices to form the 3d simulation (one 'slice' has dimension rows x columns)
int new_num_rows = -1;
int new_num_cols = -1;
int new_num_slices = -1;
REAL *dim_x;                    //The dimensions of each column in the x direction
REAL *dim_y;                    //The dimensions of each row in the y direction
REAL *dim_z;                    //The dimensions of each row in the z direction
REAL chf;                       //Constant Heat flow at base of model in mW M^2
REAL initial_time;              //The initial starting time of the model
int num_hp;                     //The number of heat production values
int num_tcd;                    //The number of thermal conductivity difference values
int num_hcf;                    //The number of fluid heat capacity values
int num_hcr;                    //The number of rock heat capacity values
int num_mtc;                    //The number of minimum convection temperature values
int num_vel;                    //The number of convection velocities
REAL *heat_production_values;   //The radioactive heat production values array used in conduction calculations
REAL *thermal_conduct_diff;     //The thermal conductivity difference array used in conduction calculations
REAL *heat_capac_fluid;         //The fluid heat capacity array used in convection calculations
REAL *heat_capac_rock;          //The rock heat capacity array used in convection calculations
REAL *min_temp_conv;            //The minimum temperature required for convection
REAL *vel;                      //The velocity array used for convection calculations


/**
 * Deallocates all allocated memory used by the program
 */
void deallocate_memory() {
     //Deletes allocated memory
    for(int i = 0; i < num_rows; i++) {
        for(int j = 0; j < num_cols; j++) {
            delete[] temp[i][j];
            delete[] cond_codes[i][j];
            if(using_convection) {
                delete[] conv_codes[i][j];
            }
        }
        delete[] temp[i];
        delete[] cond_codes[i];
        if(using_convection) {
            delete[] conv_codes[i];
        }
    }
    delete[] dim_x;
    delete[] dim_y;
    delete[] dim_z;
    delete[] temp;
    delete[] cond_codes;
    delete[] heat_production_values;
    delete[] thermal_conduct_diff;
    if(using_convection == 1) {
        delete[] conv_codes;
        delete[] heat_capac_fluid;
        delete[] heat_capac_rock;
        delete[] min_temp_conv;
        delete[] vel;
    }
}

/*
 * Clears the cin buffer
 */
void clear_cin() {
	cin.clear();
	cin.ignore(numeric_limits <streamsize> ::max(), '\n' );
}

/*
 * This function waits for the user to hit enter before continuing
 */
void PressEnterToContinue() {
    cout << "Press ENTER to continue... " << flush;
    clear_cin();
}

/*
 * Loads the input file into program memory and allocates
 * necessary memory to store the input variables
 */
void load_file() {
    ifstream source_file;       //Input file stream
    string temp_str;
    ostringstream str_conv;
 
	source_file.open(source_filename.c_str(),ios::in);
	if(!source_file.is_open()) {
		cerr << "Soruce File Not found!" << endl;
		exit(1);
	}
	
    //Loads the input file
    cout << endl << endl << "Loading Input File";
    
    //Retrieves the simulation parameters from the input file
    source_file >> num_rows >> num_cols >> num_slices >> using_convection;
    source_file >> chf >> initial_time;
    getline(source_file,title);
    getline(source_file,title);
    
    //displays parameters of the input file
    cout << endl << endl << "Number of rows   = " << num_rows << endl;
    cout << "Number of cols   = " << num_cols << endl;
    cout << "Number of slices = " << num_slices << endl;
    if(using_convection == 1) {
        cout << "Using convection" << endl;
    }
    else {
        cout << "No Convection" << endl;
    }
    cout << endl << "Constant Heat Flow at Base of Model = " << chf << "mW M^2" << endl;
    cout << "Model time elapsed = " << initial_time << " Years" << endl << endl;
    
    //Allocates memory for the conduction variables based on the previously read in simulation
    //parameters
    dim_x = new REAL[num_cols];
    dim_y = new REAL[num_rows];
    dim_z = new REAL[num_slices];
    temp = new REAL**[num_rows];
    cond_codes = new int**[num_rows];
    for(int i = 0; i < num_rows; i++) {
        temp[i] = new REAL*[num_cols];
        cond_codes[i] = new int*[num_cols];
        for (int j = 0; j < num_cols; j++) {
            temp[i][j] = new REAL[num_slices];
            cond_codes[i][j] = new int[num_slices];
        }
    }
    
    //Reads in the Y (column) dimensions and finds the minimum column distance
    for(int i = 0; i < num_cols; i++) {
        source_file >> dim_x[i];
    }

    //Reads in the X (row) dimensions and finds the minimum row distance
    for(int i = 0; i < num_rows; i++) {
        source_file >> dim_y[i];
    }

    //Reads in the Z (slice depth) dimension and finds the minimum row distance
    for (int i = 0; i < num_slices; i++) {
        source_file >> dim_z[i];
    }

    //Reads in the conduction heat production values
    source_file >> num_hp;
    heat_production_values = new REAL[num_hp];
    for(int i = 0; i < num_hp; i++) {
        source_file >> heat_production_values[i];
    }
    cout << "Read "<< num_hp << " heat production values" << endl;
    
    //Reads in the thermal conduction difference values
    //Finds the minimum and maximum thermal conductivity differences and
    //performs some scaling of the conduction associated variables
    source_file >> num_tcd;
    thermal_conduct_diff = new REAL[num_tcd];
    for(int i = 0; i < num_tcd; i++) {
        source_file >> thermal_conduct_diff[i];
    }

    //Reads in the convection specific variables if convection
    //is used by the user specified input file
    if(using_convection) {
        //Reads in the fluid heat capacity values
        source_file >> num_hcf;
        heat_capac_fluid = new REAL[num_hcf];
        for(int i = 0; i < num_hcf; i++) {
            source_file >> heat_capac_fluid[i];
        }
        
        //Reads in the rock heat capacity values
        source_file >> num_hcr;
        heat_capac_rock = new REAL[num_hcr];
        for(int i = 0; i < num_hcr; i++) {
            source_file >> heat_capac_rock[i];
        }
        
        //Reads in the minimum convection temperatures
        source_file >> num_mtc;
        min_temp_conv = new REAL[num_mtc];
        for(int i = 0; i < num_mtc; i++) {
            source_file >> min_temp_conv[i];
        }
        //Reads in the convection velocities
        source_file >> num_vel;
        vel = new REAL[num_vel];
        for(int i = 0; i < num_vel; i++) {
            source_file >> vel[i];
        }
        cout << endl << "Read " << num_vel << " Velocities in m/yr" << endl;

    }
    
    //Reads in the starting temperatures of the simulation from the input file
    for (int k = 0; k < num_slices; k++) {
        for(int i = 0; i < num_rows; i++) {
            for(int j = 0; j < num_cols; j++) {
                source_file >> temp[i][j][k];
            }
        }
    }
    cout << "Read " << num_rows << " X " << num_cols << " X " << num_slices << " temps" << endl;

    //Reads in the conduction codes for each cell of the simulation and parses
    //the array indexs from the codes
    //Unlike, the Fortran version of the program, the conduction direction codes
    //are ignored since the simulation accounts for them internally
    for (int k = 0; k < num_slices; k++) {
        for(int i = 0; i < num_rows; i++) {
            for(int j = 0; j < num_cols; j++) {
                source_file >> temp_str;
                cond_codes[i][j][k] = atoi(temp_str.c_str());
            }
        }
    }
    cout << "Read " << num_rows << " X " << num_cols << " X " << num_slices << " conduction codes" << endl;

    //If convection is used for the user specified input file, memory is allocated for its
    //variables and they are read in from the input file
    if(using_convection) {     
        //Allocates memory for the convection variables based on the previously read in simulation
        //parameters
        conv_codes = new int**[num_rows];
        for(int i = 0; i < num_rows; i++) {
            conv_codes[i] = new int*[num_cols];
            for (int j = 0; j < num_cols; j++) {
                conv_codes[i][j] = new int[num_slices];
            }
        }
        
        //Reads in the convection codes for each cell of the simulation and parses the array
        //indexs from the ocdes
        for (int k = 0; k < num_slices; k++) {
            for(int i = 0; i < num_rows; i++) {
                for(int j = 0; j < num_cols; j++) {
                    source_file >> temp_str;
                    conv_codes[i][j][k] = atoi(temp_str.c_str());
                }
            }
        }
        cout << "Read " << num_rows << " X " << num_cols << " X " << num_slices << " convection codes" << endl;
    }
    
    //Closes the input file
    source_file.close();

    cout << endl << "Done Loading Input File" << endl;
}

/*
 * Saves the current state of the simulation, using the same format
 * as the input file
 */
void save_model_state() {
    ofstream output_file;    //Output file stream
    
    //Opens the output file for writing
    output_file.open(output_filename.c_str(),ios::out);
    if(!output_file.is_open()) {
        cerr << "Failed to write state to file" << endl;
        exit(1);
    }
    else {
        //Prints the simulation parameters to the output file
        output_file << setw(20) << new_num_rows << " " << setw(20) << new_num_cols << " " << setw(20) << new_num_slices << setw(20) << using_convection << endl;
        output_file << setw(20) << fixed << setprecision(OUT_PRECISION) << chf << " " << setw(20) << initial_time << endl;
        output_file << title << endl;
        
        output_file << setfill(' ');
        output_file << setprecision(3);
        //Prints the column (X) dimensions of the simulation to the output file
        for(int i = 0; i < new_num_cols; i++) {
            output_file << " " << dim_x[i%num_cols];
        }
        output_file << endl;
        
        //Prints the row (Y) dimensions of the simulation to the output file
        for(int i = 0; i < new_num_rows; i++) {
            output_file << " " << dim_y[i%num_rows];
        }
        output_file << endl;
        
        // Prints the slice (Z) dimensions of the simulation to the oputput file

        for (int i = 0; i < new_num_slices; i++) {
            output_file << " " << dim_z[i%num_slices];
        }
        output_file << endl;

        //Prints the heat production values of the simulation to the output file
        output_file << " " << num_hp;
        for(int i = 0; i < num_hp; i++) {
            output_file << " " << scientific << heat_production_values[i];
        }
        output_file << endl;
        
        //Prints the thermal conductivity difference values to the output file
        output_file << " " << num_tcd;
        for(int i = 0; i < num_tcd; i++) {
            output_file << " " << thermal_conduct_diff[i];
        }
        output_file << endl;
        
        //Prints the convection specific variables to the output file if convection is used
        if(using_convection) {
            //Prints the fluid heat capacity values to the output file
            output_file << " " << num_hcf;
            for(int i = 0; i < num_hcf; i++) {
                output_file << " " << heat_capac_fluid[i];
            }
            output_file << endl;
            
            //Prints the rock heat capacity values to the output file
            output_file << " " << num_hcr;
            for(int i = 0; i < num_hcr; i++) {
                output_file << " " << heat_capac_rock[i];
            }
            output_file << endl;
            
            //Prints the minimum convection temps to the output file
            output_file << " " << num_mtc;
            for(int i = 0; i < num_mtc; i++) {
                output_file << " " << min_temp_conv[i];
            }
            output_file << endl;
            
            //Prints the convection velocities to the output file
            output_file << " " << num_vel;
            for(int i = 0; i < num_vel; i++) {
                output_file << " " << vel[i];
            }
            output_file << endl;
        }
        
        output_file << setprecision(OUT_PRECISION) << fixed;
        output_file << endl;
        //Prints the current temperature array of the simulation
        for (int k = 0; k < new_num_slices; k++) {
            for(int i = 0; i < new_num_rows; i++) {
                for(int j = 0; j < new_num_cols; j++) {
                    output_file << " " << setw(OUT_PRECISION+5) << temp[i%num_rows][j%num_cols][k%num_slices];
                }
                output_file << endl;
            }
            output_file << endl;
        }
		cout << "Wrote " << new_num_rows << " X " << new_num_cols << " X " << new_num_slices << " temps" << endl; 
        
        //Prints the conduction codes of the simulation to the output file
        output_file << setfill('0');
        for (int k = 0; k < new_num_slices; k++) {
            for(int i = 0; i < new_num_rows; i++) {
                for(int j = 0; j < new_num_cols; j++) {
                    output_file << " " << setw(2*INDEX_WIDTH+1) << cond_codes[i%num_rows][j%num_cols][k%num_slices];
                }
                output_file << endl;
            }
            output_file << endl;
        }
		cout << "Wrote " << new_num_rows << " X " << new_num_cols << " X " << new_num_slices << " conduction codes" << endl; 
        
        //Prints the convection codes to the output file if convection is being used
        if(using_convection) {
            for (int k = 0; k < new_num_slices; k++) {
				for(int i = 0; i < new_num_rows; i++) {
					for(int j = 0; j < new_num_cols; j++) {
                        output_file << " " << setw(4*INDEX_WIDTH+2) << conv_codes[i%num_rows][j%num_cols][k%num_slices];
                    }
                    output_file << endl;
                }
                output_file << endl;
            }
        }
		cout << "Wrote " << new_num_rows << " X " << new_num_cols << " X " << new_num_slices << " convection codes" << endl; 
        
        //Closes the output file
        output_file.close();
    }
}

void usage() {
    cerr << "Arguments: " << endl;
    cerr << "   --num_row s         :       The number of rows in the output file" << endl;
    cerr << "   --num_cols          :       The number of columns in the output file" << endl;
    cerr << "   --num_slices        :       The number of slices in the output file" << endl;
    cerr << "   --input_filename    :       The filename of the input file" << endl;
    cerr << "   --source_filename   :       The filename of the output file" << endl;
}

/*
 * Performs a finite heat flow simulation using
 * conduction and convection.
 */
int main(int argc, char **argv) {

	for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--num_rows") == 0) {
            new_num_rows = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--num_cols") == 0) {
            new_num_cols = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--num_slices") == 0) {
            new_num_slices = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--input_filename") == 0) {
            source_filename = argv[++i];
        } else if (strcmp(argv[i], "--output_filename") == 0) {
            output_filename = argv[++i];
        } else {
            cerr << "Unknown argument " << argv[i] << endl;
            usage();
            exit(0);
        }
    }
    
	if(new_num_rows <= 0) {
		cerr << "New number of rows not specified" << endl;
		exit(1);
	}
	if(new_num_cols <= 0) {
		cerr << "New number of columns not specified" << endl;
		exit(1);
	}
	if(new_num_slices <= 0) {
		cerr << "New number of slices not specified" << endl;
		exit(1);
	}
	if(source_filename.size() <= 0) {
		cerr << "Source Filename not specified" << endl;
		exit(1);
	}
	if(output_filename.size() <= 0) {
		cerr << "Source Filename not specified" << endl;
		exit(1);
	}
	
    //Loads the input file for the simulation
    load_file();

	save_model_state();

	deallocate_memory();
      
}
