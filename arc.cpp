/**
 * Performs a finite difference heat flow 
 * simulation using conduction and convection.
 */

#ifdef _WIN32
#include <windows.h>
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <limits>
#include <math.h>

/**
 * Time_step (seconds/yr) divided by the product of density and heat capacity.
 * The value for the density is 2200 kg/m^3 and heat capacity is 1000 kJ/kg K.
 */
#define QFAC 14.33              //Description is defined in the previous comment
#define DTC 0.25                //
#define OUT_PRECISION 10        //Number of digits to print after the decimal place for floating point values
#define REAL float              //The precision of the model.

using std::cerr;
using std::cin;
using std::cout;
using std::endl;
using std::string;
using std::ofstream;
using std::ifstream;
using std::setw;
using std::right;
using std::left;
using std::fixed;
using std::scientific;
using std::setprecision;
using std::ios;
using std::numeric_limits;
using std::streamsize;
using std::max;
using std::flush;


//Conduction code specific variables
int **cond_codes;               //The unmodified conduction codes as read from the input file
int **cond_hp_index;            //The conduction index for the radioactive heat production array
int **cond_tc_index;            //The conduction index for the thermal conductivity array
REAL DHF;                       //

/**
 * Convection Directions
 * 1  2  3
 * 4  5  6
 * 7  8  9
 */
 //Convection code specific variables
int **conv_codes;               //The convection codes as read from the input file
int **conv_min_temp_index;      //The convection index for the minimum temp for convection array
int **conv_direction;           //The direction of convection following the direction matrix in the previous comment
int **conv_vel_index;           //The convection index for the velocity array
int **conv_fluid_index;         //The convection index for the fluid heat capacity array
int **conv_rock_index;          //The convection index for the rock heat capacity array
int num_conv_loops;             //The number of convection updates to perform per time step
REAL time_inc;                  //The amount of time increment per convection loop

//File names
string source_filename;         //The input files name with extension
string output_filename;         //The output state files name with extension
string output_su_filename;      //The output surfer files name with extension

//Input file variables
string title;                   //The title of the input file
int using_convection;           //Indicates if convection is being used
REAL **temp;                    //The current temperature array
int num_rows;                   //The number of rows for the simulation
int num_cols;                   //The number of columns for the simulation
REAL *dim_x;                    //The dimensions of each column in the x direction
REAL *dim_y;                    //The dimensions of each row in the y direction
REAL chf;                       //Constant Heat flow at base of model in mW M^2
REAL initial_time;              //The initial starting time of the model
REAL heat_production_values[8]; //The radioactive heat production values array used in conduction calculations
REAL thermal_conduct_diff[8];   //The thermal conductivity difference array used in conduction calculations
REAL heat_capac_fluid[8];       //The fluid heat capacity array used in convection calculations
REAL heat_capac_rock[8];        //The rock heat capacity array used in convection calculations
REAL min_temp_conv[8];          //The minimum temperature required for convection
REAL vel[8];                    //The velocity array used for convection calculations

//Moving source variables
int using_moving_source;        //Indicates if a moving source is being used
int mvsrc_start_col;            //Starting column of the moving source
int mvsrc_end_col;              //Ending column of the moving source
int mvsrc_depth;                //The depth of the moving source (defines the row in the simulation)
REAL mvsrc_temp;                //The temperature of the moving source

//Time specific variables
REAL sim_time;                  //The current time of the simulation
REAL tic;                       //Time variable used in convection calculations
REAL time_step;                 //The amount of time that passes between each update of the simulation, time step
REAL run_time;                  //The run time of the simulation

//Global variables
int save_result;                //Indicates if the model should save the current state at each screen update
REAL max_vel;                   //The maximum convection velocity of the velocity array
REAL min_row_dim;               //The minimum y dimension of each cell of the simulation
REAL min_col_dim;               //The minimum x dimension of each cell of the simulation
REAL thermal_time_constant;     //The thermal time constant of the model, used in the selection of the run time
REAL **next_temp;               //The next temperature array
REAL max_thermal_conduct_diff;  //The maximum thermal conductivity difference
REAL min_thermal_conduct_diff;  //The minimum thermal conductivity difference
int num_loops;                  //The number of loops between screen updates

/**
 * This function waits for the user to hit enter before continuing
 */
void PressEnterToContinue() {
    cout << "Press ENTER to continue... " << flush;
    cin.ignore(numeric_limits <streamsize> ::max(), '\n' );
}

/**
 * Swaps the temp arrays
 */
void swap_temp_array() {
    REAL **tmp;
    
    tmp = temp;
    temp = next_temp;
    next_temp = tmp;
}

/**
 * Loads the input file into program memory and allocates
 * necessary memory to store the input variables
 */
void load_file() {
    ifstream source_file;       //Input file stream
    
    //Ask for the input file names and displays an error message
    //if the file does not exist
    do {
        cout << "Input File Name: ";
        cin >> source_filename;
        source_file.open(source_filename.c_str(),ios::in);
        if(!source_file.is_open()) {
            cout << "File Not found!" << endl;
        }
    } while(!source_file.is_open());
    
    //Asks the user if the state of the model should be saved every screen update
    cout << endl << "To save the state of the model every screen update enter 1, otherwise 0: ";
    cin >> save_result;
    while(save_result < 0 || save_result > 1) {
        cout << "Incorrect input, to save the state of the model enter 1, else 0: ";
        cin >> save_result;
    }
    
    //Ask for the state filename if the user specified that the file should be saved
    if(save_result == 1) {
        cout << "Output filename: ";
        cin >> output_filename;
        while(output_filename.size() <= 0) {
            cout << "Incorrect input, please enter a valid filename: ";
            cin >> output_filename;
        }
    }
    
    //Asks for the DSAA surfer grid filenmae
    cout << "Surfer filename: ";
    cin >> output_su_filename;
    while(output_su_filename.size() <= 0) {
        cout << "Incorrect input, please enter a valid filename: ";
    }
    
    //Loads the input file
    cout << endl << endl << "Loading Input File";
    
    //Retrieves the simulation parameters from the input file
    source_file >> num_rows >> num_cols >> using_convection;
    source_file >> chf >> initial_time;
    source_file >> title;
    
    //displays parameters of the input file
    cout << endl << endl << "Number of rows  = " << num_rows << endl;
    cout << "Number of cols  = " << num_cols << endl;
    if(using_convection == 1) {
        cout << "Using convection" << endl;
    }
    else {
        cout << "No Convection" << endl;
    }
    cout << endl << "Constant Heat Flow at Base of Model = " << chf << "mW M^2" << endl;
    chf *= 0.001;
    cout << "Model time elapsed = " << initial_time << " Years" << endl << endl;
    
    //Allocates memory for the conduction variables based on the previously read in simulation
    //parameters
    dim_x = new REAL[num_cols];
    dim_y = new REAL[num_rows];
    temp = new REAL*[num_rows];
    next_temp = new REAL*[num_rows];
    cond_codes = new int*[num_rows];
    cond_hp_index = new int*[num_rows];
    cond_tc_index = new int*[num_rows];
    for(int i = 0; i < num_rows; i++) {
        temp[i] = new REAL[num_cols];
        next_temp[i] = new REAL[num_cols];
        cond_codes[i] = new int[num_cols];
        cond_hp_index[i] = new int[num_cols];
        cond_tc_index[i] = new int[num_cols];
    }
    //Reads in the starting temperatures of the simulation from the input file
    for(int i = 0; i < num_rows; i++) {
        for(int j = 0; j < num_cols; j++) {
            source_file >> temp[i][j];
        }
    }
    cout << "Read " << num_rows << " X " << num_cols << " temps" << endl;

    //Reads in the conduction codes for each cell of the simulation and parses
    //the array indexs from the codes
    //Unlike, the Fortran version of the program, the conduction direction codes
    //are ignored since the simulation accounts for them internally
    for(int i = 0; i < num_rows; i++) {
        for(int j = 0; j < num_cols; j++) {
            source_file >> cond_codes[i][j];
            cond_tc_index[i][j] = cond_codes[i][j]/1000 - 1;
            cond_hp_index[i][j] = (cond_codes[i][j] - (cond_tc_index[i][j]+1)*1000)/100 - 1;
        }
    }
    cout << "Read " << num_rows << " X " << num_cols << " conduction codes" << endl;

    //If convection is used for the user specified input file, memory is allocated for its
    //variables and they are read in from the input file
    if(using_convection) {
        int tmp_val = 0.0;
        
        //Allocates memory for the convection variables based on the previously read in simulation
        //parameters
        conv_codes = new int*[num_rows];
        conv_min_temp_index = new int*[num_rows];
        conv_direction = new int*[num_rows];
        conv_vel_index = new int*[num_rows];
        conv_fluid_index = new int*[num_rows];
        conv_rock_index = new int*[num_rows];
        for(int i = 0; i < num_rows; i++) {
            conv_codes[i] = new int[num_cols];
            conv_min_temp_index[i] = new int[num_cols];
            conv_direction[i] = new int[num_cols];
            conv_vel_index[i] = new int[num_cols];
            conv_fluid_index[i] = new int[num_cols];
            conv_rock_index[i] = new int[num_cols];
        }
        
        //Reads in the convection codes for each cell of the simulation and parses the array
        //indexs from the ocdes
        for(int i = 0; i < num_rows; i++) {
            for(int j = 0; j < num_cols; j++) {
                source_file >> conv_codes[i][j];
                tmp_val = conv_codes[i][j];
                conv_min_temp_index[i][j] = conv_codes[i][j] / 10000 - 1;
                tmp_val -= 10000*(conv_min_temp_index[i][j]+1);
                conv_vel_index[i][j] = tmp_val/1000 - 1;
                tmp_val -= 1000*(conv_vel_index[i][j]+1);
                conv_fluid_index[i][j] = tmp_val/100 - 1;
                tmp_val -= 100*(conv_fluid_index[i][j]+1);
                conv_rock_index[i][j] = tmp_val/10 - 1;
                conv_direction[i][j] = tmp_val - 10*(conv_rock_index[i][j]+1);
            }
        }
        cout << "Read " << num_rows << " X " << num_cols << " convection codes" << endl;
    }
    //Reads in the Y (column) dimensions and finds the minimum column distance
    for(int i = 0; i < num_cols; i++) {
        source_file >> dim_x[i];
        if(i == 0) {
            min_col_dim = dim_x[i];
        }
        else if(dim_x[i] < min_col_dim) {
            min_col_dim = dim_x[i];
        }
    }

    //Reads in the X (row) dimensions and finds the minimum row distance
    for(int i = 0; i < num_rows; i++) {
        source_file >> dim_y[i];
        if(i == 0) {
            min_row_dim = dim_y[i];
        }
        else if(dim_y[i] < min_row_dim) {
            min_row_dim = dim_y[i];
        }
    }
    
    //Reads in the conduction heat production values
    for(int i = 0; i < 8; i++) {
        source_file >> heat_production_values[i];
    }
    cout << "Read 8 heat production values" << endl;
    
    //Reads in the thermal conduction difference values
    for(int i = 0; i < 8; i++) {
        source_file >> thermal_conduct_diff[i];
    }
    cout << "Converted 8 Thermal Conductivities to Diff. in m^2/y" << endl;
    
    //Finds the minimum and maximum thermal conductivity differences and
    //performs some scaling of the conduction associated variables
    max_thermal_conduct_diff = thermal_conduct_diff[0];
    min_thermal_conduct_diff = thermal_conduct_diff[0];
    for(int i = 0; i < 8; i++) {
        heat_production_values[i] /= 1E6;
        thermal_conduct_diff[i] *= 14.33;

        if(thermal_conduct_diff[i] > max_thermal_conduct_diff) {
            max_thermal_conduct_diff = thermal_conduct_diff[i];
        }
        if(thermal_conduct_diff[i] < min_thermal_conduct_diff) {
            min_thermal_conduct_diff = thermal_conduct_diff[i];
        }
        cout << "  " << thermal_conduct_diff[i];
    }
    
    //Reads in the convection specific variables if convection
    //is used by the user specified input file
    if(using_convection) {
    
        //Reads in the fluid heat capacity values
        for(int i = 0; i < 8; i++) {
            source_file >> heat_capac_fluid[i];
        }
        
        //Reads in the rock heat capacity values
        for(int i = 0; i < 8; i++) {
            source_file >> heat_capac_rock[i];
        }
        
        //Reads in the minimum convection temperatures
        for(int i = 0; i < 8; i++) {
            
            source_file >> min_temp_conv[i];
        }
        //Reads in the convection velocities
        for(int i = 0; i < 8; i++) {
            source_file >> vel[i];
        }
        cout << endl << "Read 8 Velocities in m/yr" << endl;
        
        //Finds the maximum convection velocity
        max_vel = vel[0];
        cout << " " << vel[0];
        for(int i = 1 ; i < 8; i++) {
            if(vel[i] > max_vel) {
                max_vel = vel[i];
            }
            cout << " " << vel[i];
        }
        cout << endl;
    }
    
    //Closes the input file
    source_file.close();
    
    /*
    T1 = max_thermal_conduct_diff;
    T2 = min_col_dim;
    T3 = min_row_dim;
    */
    //Finds the convection time increment
    if(using_convection) {
        if(min_col_dim > min_row_dim) {
            tic = min_row_dim/max_vel;
        }
        else {
            tic = min_col_dim/max_vel;
        }
    }
    //Calculates the maximum time step of the simulation
    if(min_col_dim < min_row_dim) {
        time_step = min_col_dim*min_col_dim/(5*max_thermal_conduct_diff);
    }
    else {
        time_step = min_row_dim*min_row_dim/(5*max_thermal_conduct_diff);
    }
    cout << endl << "Done Loading Input File" << endl;
}


/**
 * Saves the current state of the simulation, using the same format
 * as the input file
 */
void save_state() {
    ofstream output_file;    //Output file stream
    
    //Opens the output file for writing
    output_file.open(output_filename.c_str(),ios::out);
    if(!output_file.is_open()) {
        cerr << "Failed to write state to file" << endl;
        exit(1);
    }
    else {
        //Prints the simulation parameters to the output file
        output_file << setw(20) << num_rows << " " << setw(20) << num_cols << " " << setw(20) << using_convection << endl;
        output_file << setw(20) << fixed << setprecision(2) << chf*1000.0 << " " << setw(20) << initial_time + sim_time << endl;
        output_file << title << endl;
        
        output_file << setprecision(OUT_PRECISION);
        //Prints the current temperature array of the simulation
        for(int i = 0; i < num_rows; i++) {
            for(int j = 0; j < num_cols; j++) {
                output_file << " " << setw(OUT_PRECISION+5) << temp[i][j];
            }
            output_file << endl;
        }
        
        //Prints the conduction codes of the simulation to the output file
        for(int i = 0; i < num_rows; i++) {
            for(int j = 0; j < num_cols; j++) {
                output_file << " " << cond_codes[i][j];
            }
            output_file << endl;
        }
        
        //Prints the convection codes to the output file if convection is being used
        if(using_convection) {
            for(int i = 0; i < num_rows; i++) {
                for(int j = 0; j < num_cols; j++) {
                    output_file << " " << conv_codes[i][j];
                }
                output_file << endl;
            }
        }
        
        output_file << setprecision(3);
        //Prints the column (X) dimensions of the simulation to the output file
        for(int i = 0; i < num_cols; i++) {
            output_file << " " << dim_x[i];
        }
        output_file << endl;
        
        //Prints the row (Y) dimensions of the simulation to the output file
        for(int i = 0; i < num_rows; i++) {
            output_file << " " << dim_y[i];
        }
        output_file << endl;
        //Prints the heat production values of the simulation to the output file
        for(int i = 0; i < 8 ; i++) {
            output_file << " " << scientific << heat_production_values[i]*1E6;
        }
        output_file << endl;
        
        //Prints the thermal conductivity difference values to the output file
        for(int i = 0; i < 8 ; i++) {
            output_file << " " << thermal_conduct_diff[i]/14.33;
        }
        output_file << endl;
        
        //Prints the convection specific variables to the output file if convection is used
        if(using_convection) {
            //Prints the fluid heat capacity values to the output file
            for(int i = 0; i < 8 ; i++) {
                output_file << " " << heat_capac_fluid[i];
            }
            output_file << endl;
            
            //Prints the rock heat capacity values to the output file
            for(int i = 0; i < 8 ; i++) {
                output_file << " " << heat_capac_rock[i];
            }
            output_file << endl;
            
            //Prints the minimum convection temps to the output file
            for(int i = 0; i < 8 ; i++) {
                output_file << " " << min_temp_conv[i];
            }
            output_file << endl;
            
            //Prints the convection velocities to the output file
            for(int i = 0; i < 8 ; i++) {
                output_file << " " << vel[i];
            }
            output_file << endl;
        }
        
        //Closes the output file
        output_file.close();
    }
}

/**
 * Saves the current temperatures of the simulation to a DSAA surfer grid file
 */
void save_surfer() {
    ofstream output_file;            //Output file stream
    
    //Opens the output file for writting
    output_file.open(output_su_filename.c_str(),ios::out);
    if(!output_file.is_open()) {
        cerr << "Failed to write surfer file" << endl;
        exit(1);
    }
    else {
        REAL min_temp, max_temp, temp_range;    //Minimum and maximum temperatures.
        REAL xmax,ymin;                         //Maximum x and minimum y distances
        
        //Finds the minimum and maximum temps in the temperature array
        min_temp = max_temp = temp[0][0];
        for(int i = 0; i < num_rows; i++) {
            for(int j = 0; j < num_cols; j++) {
                if(temp[i][j] > max_temp) {
                    max_temp = temp[i][j];
                }
                if(temp[i][j] < min_temp) {
                    min_temp = temp[i][j];
                }
            }
        }
        
        //Calculates the temperature range.
        temp_range = max_temp - min_temp;
        if(temp_range == 0) {
            temp_range = 1.0;
        }
        
        //Calculates the maximum x distance and the
        //minimum y distance
        xmax = dim_x[0]*num_cols;
        ymin = dim_y[0]*num_rows;
        if(dim_x[0] < 0.01) {
            xmax *= 1000;
        }
        else if(dim_x[0] < 0.1) {
            xmax *= 100;
        }
        else if(dim_x[0] < 1) {
            xmax *= 10;
        }
        if(dim_y[0] < 0.01) {
            ymin *= 1000;
        }
        else if(dim_y[0] < 0.1) {
            ymin *= 100;
        }
        else if(dim_y[0] < 1) {
            ymin *= 10;
        }
        
        //Prints the DSAA surfer grid parameters to the output file
        output_file << "DSAA" << endl;
        output_file << setw(20) << num_cols << " " << setw(20) << num_rows << endl;
        output_file << fixed << setprecision(3) << setw(20) << 0.0 << " " << setw(20) << xmax << endl;
        output_file << setw(20) << -ymin << " " << setw(20) << 0.0 << endl;
        output_file << setw(20) << setprecision(OUT_PRECISION) << min_temp << " " << setw(20) << max_temp << endl;
        
        //Prints the temperature array to the output file
        for(int i = num_rows-1; i >= 0; i--) {
            for(int j = 0; j < num_cols; j++) {
                output_file << " " << setw(OUT_PRECISION+5) << temp[i][j];
            }
            output_file << endl;
        }
        
        //Closes the output file
        output_file.close();
    }
}

/**
 * Calculates and returns the heat flow per year between two cells in the X direction
 * based on the provided indexes
 */
REAL cond_add_x_2D(int x1, int y1, int x2, int y2) {
    REAL temp_diff;    //Temperature difference between the two cells
    REAL ad;           //
    
    temp_diff = temp[x1][y2] - temp[x1][y1];
    ad = dim_x[y2]/thermal_conduct_diff[cond_tc_index[x1][y2]] + dim_x[y1]/thermal_conduct_diff[cond_tc_index[x1][y1]];
    
    return 2*temp_diff/(ad*dim_x[y1]);
}

/**
 * Calculates and returns the heat flow per year between two cells in the Y direction
 * based on the provided indexes
 */
REAL cond_add_y_2D(int x1, int y1, int x2, int y2) {
    REAL temp_diff;    //Temperature difference between the two cells
    REAL ad;           //
    
    temp_diff = temp[x2][y1] - temp[x1][y1];
    ad = dim_y[x2]/thermal_conduct_diff[cond_tc_index[x2][y1]] + dim_y[x1]/thermal_conduct_diff[cond_tc_index[x1][y1]];
    
    return 2*temp_diff/(ad*dim_y[x1]);
}

/**
 * Updates the temperature array using 2D conduction with finite
 * difference heat flow.
 */
void cond_2D(){
    REAL heat_flow_x;    //Heat flow into the cell in the x direction per year
    REAL heat_flow_y;    //Heat flow into the cell in the y direction per year
    
    for(int i = 0; i < num_rows; i++) {
        for(int j = 0; j < num_cols; j++) {
            //The top of the model, i.e. row 0's values do not change due to conduction during the simulation
            if(i == 0) {
                next_temp[i][j] = temp[i][j];
            }
            else {
                //Calculates heat flow in the X and Y direction into the current
                //cell based on its location within the model
                if(i == num_rows-1 && j == 0) { //Bottom left corner
                    heat_flow_x = cond_add_x_2D(i,j,i,j+1);
                    heat_flow_y = cond_add_y_2D(i,j,i-1,j);
                    next_temp[i][j] = DHF/dim_y[i];    //Constant heat flow at the bottom of the model
                }
                else if(i == num_rows-1 && j == num_cols-1) { //Bottom right corner
                    heat_flow_x = cond_add_x_2D(i,j,i,j-1);
                    heat_flow_y = cond_add_y_2D(i,j,i-1,j);
                    next_temp[i][j] = DHF/dim_y[i];    //Constant heat flow at the bottom of the model
                }
                else if(i == num_rows-1) { //Bottom
                    heat_flow_x = cond_add_x_2D(i,j,i,j+1) + cond_add_x_2D(i,j,i,j-1);
                    heat_flow_y = cond_add_y_2D(i,j,i-1,j);
                    next_temp[i][j] = DHF/dim_y[i];    //Constant heat flow at the bottom of the model
                }
                else if(j == 0) { //Left side
                    heat_flow_x = cond_add_x_2D(i,j,i,j+1);
                    heat_flow_y = cond_add_y_2D(i,j,i-1,j) + cond_add_y_2D(i,j,i+1,j);
                    next_temp[i][j] = 0.0;
                }
                else if(j == num_cols-1) { //Right side
                    heat_flow_x = cond_add_x_2D(i,j,i,j-1);
                    heat_flow_y = cond_add_y_2D(i,j,i-1,j) + cond_add_y_2D(i,j,i+1,j);
                    next_temp[i][j] = 0.0;
                }
                else { //Middle
                    heat_flow_x = cond_add_x_2D(i,j,i,j-1) + cond_add_x_2D(i,j,i,j+1);
                    heat_flow_y = cond_add_y_2D(i,j,i-1,j) + cond_add_y_2D(i,j,i+1,j);
                    next_temp[i][j] = 0.0;
                }
                
                //Heat flow from the adjacent cells
                next_temp[i][j] += temp[i][j] + time_step*(heat_flow_x + heat_flow_y);
                //Heat flow due to radioactive heat production
                next_temp[i][j] += heat_production_values[cond_hp_index[i][j]]*time_step/DTC;
            }
        }
    }
    swap_temp_array();    //Swaps the current and next temperature arrays
}

/**
 * Performs convection between two specified cells
 */
void perform_conv_2D(REAL time_inc, int x1, int y1, int x2, int y2) {
    REAL avg_x_dim;    //Average X dimension for the two cells
    REAL avg_y_dim;    //Average Y dimension for the two cells
    REAL amt;          //
    REAL dist;         //Distance between the two cells
    REAL ratio;        //Ratio of amt to distance
    
    //Checks if the specified cell is within the bounds of the simulation and if it has a high enough
    //temperature to perform convection
    if((x2 >= 0) && (x2 < num_rows) && (y2 >= 0) && (y2 < num_cols) && (temp[x2][y2] - min_temp_conv[conv_min_temp_index[x1][y1]] >= 0)) {
        avg_x_dim = (dim_x[y2] + dim_x[y1])/2.0;
        avg_y_dim = (dim_y[x2] + dim_y[x1])/2.0;
        amt = (vel[conv_vel_index[x1][y1]]*heat_capac_fluid[conv_fluid_index[x1][y1]]/heat_capac_rock[conv_rock_index[x1][y1]])*time_inc;
        dist = sqrt(avg_x_dim*avg_x_dim + avg_y_dim*avg_y_dim);
        ratio = amt/dist;
        if(ratio > 1) {
            ratio = 0.999999;
        }
        next_temp[x1][y1] = temp[x1][y1] + ratio *(temp[x2][y2]-temp[x1][y1]);
    }
    else {
        next_temp[x1][y1] = temp[x1][y1];
    }
}

/**
 * Updates the temperature array using convection
 */
void conv_2D() {
    //Performs the calculated number of convection updates per time step
    for(int k = 0; k < num_conv_loops; k++) {
        for(int i = 0; i < num_rows; i++) {
            for(int j = 0; j < num_cols; j++) {
                //Checks if convection can occur for the specified cell
                if((conv_codes[i][j] <= 0) || (i == 0) || (conv_direction[i][j] == 5) || (conv_direction[i][j] < 1) || (conv_direction[i][j] > 9)) {
                    next_temp[i][j] = temp[i][j];
                }
                else {
                    //Performs convection based on the convection direction code
                    switch(conv_direction[i][j]) {
                        case 1:
                            perform_conv_2D(time_inc,i,j,i-1,j-1);
                            break;
                        case 2:
                            perform_conv_2D(time_inc,i,j,i-1,j);
                            break;
                        case 3:
                            perform_conv_2D(time_inc,i,j,i-1,j+1);
                            break;
                        case 4:
                            perform_conv_2D(time_inc,i,j,i,j-1);
                            break;
                        case 6:
                            perform_conv_2D(time_inc,i,j,i,j+1);
                            break;
                        case 7:
                            perform_conv_2D(time_inc,i,j,i+1,j-1);
                            break;
                        case 8:
                            perform_conv_2D(time_inc,i,j,i+1,j);
                            break;
                        case 9:
                            perform_conv_2D(time_inc,i,j,i+1,j+1);
                            break;
                    }
                }
            }
        }
        swap_temp_array();    //Swaps the current and next temperature arrays
    }
}

/**
 * Performs a finite heat flow simulation using
 * conduction and convection.
 */
int main(int argc, char **argv) {
    int input_val;    //Temporary int value
    REAL temp_val;    //Temporary REAL value
    
    cout << "\t\t Finite Difference Heat Flow Simulation" << endl;
    
    //Loads the input file for the simulation
    load_file();

    /**
     * Allows the user to change multiple rectangular blocks of temperatures
     * within the model
     */
    cout << endl << endl << "To Change the Temp. on a Block, Enter 1, Else 0: ";
    cin >> input_val;
    while(!(input_val == 0 || input_val == 1)) {
        cout << "Incorrect Input, Enter 1 to Change, Else 0: ";
        cin >> input_val;
    }
    
    //Warning, the row column pairs need to be space seperated not comma seperated
    if(input_val == 1) {
        int num_block, row1, row2, col1, col2;
        REAL new_temp;
        cout << "Enter the Number of Blocks to Change: ";
        cin >> num_block;
        for(int i = 0; i < num_block; i++) {
            cout << endl << "Block " << i << endl;
            cout << "Enter the Coordinates of the Upper Left Corner <row> <column>: ";
            cin >> row1 >> col1;
            cout << "Enter the Coordinates of the Lower Right Corner <row> <column>: ";
            cin >> row2 >> col2;
            cout << endl << "Current Block Temps" << endl;
            cout << setw(10) << "row" << " " << setw(10) << "col" << " " << setw(OUT_PRECISION+5) << "temp" << endl;
            cout << setw(10) << row1 << " " << setw(10) << col1 << " " << setw(OUT_PRECISION+5) << fixed << setprecision(OUT_PRECISION) << temp[row1][col1] << endl;
            cout << setw(10) << row2 << " " << setw(10) << col2 << " " << setw(OUT_PRECISION+5) << temp[row2][col2] << endl;
            cout << "Enter a New Temperature For the Block: ";
            cin >> new_temp;
            for(int i = row1; i < row2; i++) {
                for(int j = col1; j < col2; j++) {
                    temp[i][j] = new_temp;
                }
            }
        }
    }
    
    /**
     * Allows the user to start a moving source.
     * the only direction currently implemented is from left to right and
     * the move is made per iteration. Needs to be updated to
     * allow for multidirectional movement and better velocity control
     */
    cout << endl << endl << "To Start a Moving Source Enter 1, Else Enter 0: ";
    cin >> using_moving_source;
    while(!(using_moving_source == 0 || using_moving_source == 1)) {
        cout << "Incorrect Input, Enter 1 to Change, Else 0: ";
        cin >> using_moving_source;
    }
    if(using_moving_source == 1) {
        cout << "Enter Starting and Ending Columns: ";
        cin >> mvsrc_start_col >> mvsrc_end_col;
        cout << "Maximum Depth is " << (int)num_rows*dim_y[0] << " Meters" << endl;
        cout << "Enter the Depth to Source in Meters: ";
        cin >> mvsrc_depth;
        mvsrc_depth  = (int)mvsrc_depth/dim_y[0];
        cout << "Enter the Temperature of the Moving Source: ";
        cin >> mvsrc_temp;
    }
    
    //Allows the user to decrease the size of the time step
    cout << endl << endl << "Each Iteration in Time Spans " << scientific << time_step << " Years" << endl;
    cout << "Enter a Shorter Iteration Time in Years if Desired: ";
    cin >> temp_val;
    if(temp_val < time_step) {
        time_step = temp_val;
    }

    DHF = chf * QFAC * time_step;

    //Calculates the number of convection loops to perform per time step
    num_conv_loops = (int)(time_step/(10*tic));
    if(num_conv_loops > 5) {
        num_conv_loops = 5;
    }
    else if(num_conv_loops <= 0) {
        num_conv_loops = 1;
    }
    
    //Calculates the time increment per convection loop
    time_inc = time_step/num_conv_loops;
    
    min_row_dim = 100.0;
    for(int i = 0; i < num_rows; i++) {
        if(dim_y[i] < min_row_dim) {
            min_row_dim = dim_y[i];
        }
    }
    
    //Asks the user for the runtime of the simulation
    thermal_time_constant = min_row_dim*min_row_dim/max_thermal_conduct_diff;
    cout << endl << endl << "The Thermal Time Constant for the Vertical Dimension is " << thermal_time_constant << " Years" << endl;
    cout << "Enter Time Duration for Calculation in Years: ";
    cin >> run_time;
    
    //Asks the user for the number of loops to perform between screen updates
    cout << endl << endl << "Enter the Number of Loops Between Screen Updates: ";
    cin >> num_loops;
    
    //Initializes the simulation time to 0.0
    sim_time = 0.0;
    
    //Waits for the user to hit enter before beginning the simulation
    cout << endl;
    cin.ignore(numeric_limits <streamsize> ::max(), '\n' );
    PressEnterToContinue();
    
    /**
     * The main loop of the simulation
     */
    unsigned long count = 0;    //Number of loops performed
    cout << endl << endl << num_loops << " loops between screen updates" << endl << endl;
    cout << setw(15) << "num loops" << setw(20) << "run time (years)" << setw(20) << "sim time (years)" << endl;
    while(sim_time <= run_time) {
        //Displays status information for the current loop
        if(count%num_loops == 0) {
            cout << setw(15) << count << setw(20) << fixed << setprecision(5) << sim_time << setw(20) << initial_time + sim_time << endl;
            
            //Saves the current state of the simulation if the save_result flag is set
            if(save_result) {
                save_state();
            }
        }
        
        //Updates the moving source
        if(using_moving_source) {
            if(mvsrc_start_col < mvsrc_end_col) {
                mvsrc_start_col++;
                temp[mvsrc_depth][mvsrc_start_col] = mvsrc_temp;
            }
        }
        
        //Performs convection updates if the current simulation is using convection
        if(using_convection) {
            conv_2D();
        }
        
        //Performs conduction calculations
        cond_2D();
        
        //Increments the simulation time and loop count
        sim_time += time_step;
        count++;
    }
    
    //Saves the final result of the simulation
    if(save_result) {
        save_state();
    }
    save_surfer();
    
    //Deletes allocated memory
    for(int i = 0; i < num_rows; i++) {
        delete[] temp[i];
        delete[] next_temp[i];
        delete[] cond_codes[i];
        delete[] cond_hp_index[i];
        delete[] cond_tc_index[i];
        delete[] conv_codes[i];
        delete[] conv_min_temp_index[i];
        delete[] conv_direction[i];
        delete[] conv_vel_index[i];
        delete[] conv_fluid_index[i];
        delete[] conv_rock_index[i];
    }
    delete[] dim_x;
    delete[] dim_y;
    delete[] temp;
    delete[] next_temp;
    delete[] cond_codes;
    delete[] cond_hp_index;
    delete[] cond_tc_index;
    delete[] conv_codes;
    delete[] conv_min_temp_index;
    delete[] conv_direction;
    delete[] conv_vel_index;
    delete[] conv_fluid_index;
    delete[] conv_rock_index;
        
    //Waits for the user to hit enter before ending the simulation
    cout << endl << "Simulation Complete" << endl;
    PressEnterToContinue();
}
