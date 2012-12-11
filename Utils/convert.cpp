/**
 * Converts old 2D simulation files to
 * the new 3D format.
 */

#ifdef _WIN32
#include <windows.h>
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <limits>

#define OUT_PRECISION 2         //Number of digits to print after the decimal place for floating point values
#define INDEX_WIDTH 2           //The number of characters to print and read for each conduction and convection code
#define REAL double             //The precision of the model.

using std::cerr;
using std::cout;
using std::cin;
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
using std::setfill;
using std::ios;

/**
 * Loads an old 2D simulation input file and converts
 * it to the new 3D format.
 * Takes two command line arguments: the input filename and
 * output filename
 */
int main(int argc, char **argv) {
    string in_filename;     //The input filename
    string out_filename;    //The output filename
    ifstream in_file;       //Input file stream
    ofstream out_file;      //Output file stream
    int temp_int;           //Temporary int
    REAL temp_real;         //Temporary REAL
    string temp_str, title; //Temporary string
    
    int num_rows;           //The number of rows for the simulation
    int num_cols;           //The number of columns for the simulation
    int num_slices = 1;     //Total number of slices to form the 3d simulation (one 'slice' has dimension rows x columns)
    int using_convection;   //Indicates if convection is being used
    
    //Opens the input file for reading
    //Exits the program if it fails to open
    do {
        cout << "Input File Name: ";
        cin >> in_filename;
        in_file.open(in_filename.c_str(),ios::in);
        if(!in_file.is_open()) {
            cout << "File Not found!" << endl;
        }
    } while(!in_file.is_open());
    
    //Opens the output file for reading
    //Exits the program if it fails to open
    do {
        cout << "Output File Name: ";
        cin >> out_filename;
        out_file.open(out_filename.c_str(),ios::out);
        if(!out_file.is_open()) {
            cout << "Unable to create output file" << endl;
        }
    } while(!out_file.is_open());
    
    //num_rows num_cols num_slices using_convection
    in_file >> num_rows >> num_cols >> using_convection;
    out_file << setw(20) << num_rows << setw(20) << num_cols << setw(20) << 1 << setw(20) << using_convection << endl;
    
    //chf
    in_file >> temp_real;
    out_file << setw(20) << fixed << setprecision(2) << temp_real << " ";
    
    //initial_time
    in_file >> temp_real;
    out_file << setw(20) << temp_real << endl;
    
    //title
    getline(in_file,title);
    getline(in_file,title);
    out_file << title << endl;
    
    //Temperatures
    out_file << setprecision(OUT_PRECISION);
    for(int i = 0; i < num_rows; i++) {
        for(int j = 0; j < num_cols; j++) {
            in_file >> temp_real;
            out_file << " " << setw(OUT_PRECISION+5) << temp_real;
        }
        out_file << endl;
    }
    out_file << endl;
    
    //Conduction Codes
    out_file << setfill('0');
    for(int i = 0; i < num_rows; i++) {
        for(int j = 0; j < num_cols; j++) {
            in_file >> temp_str;
            temp_int = atoi(temp_str.substr(2,2).c_str());
            out_file << " " << setw(INDEX_WIDTH) << atoi(temp_str.substr(0,1).c_str()) << setw(INDEX_WIDTH) << atoi(temp_str.substr(1,1).c_str());
            if(temp_int == 2) {
                out_file << setw(1) << 0;
            }
            else {
                out_file << setw(1) << 1;
            }
        }
        out_file << endl;
    }
    out_file << endl;
    
    if(using_convection) {
        //Convection Codes
        for(int i = 0; i < num_rows; i++) {
            for(int j = 0; j < num_cols; j++) {
                in_file >> temp_str;
                out_file << " " << setw(INDEX_WIDTH) << atoi(temp_str.substr(0,1).c_str()); 
                out_file << setw(INDEX_WIDTH) << atoi(temp_str.substr(1,1).c_str());
                out_file << setw(INDEX_WIDTH) << atoi(temp_str.substr(2,1).c_str());
                out_file << setw(INDEX_WIDTH) << atoi(temp_str.substr(3,1).c_str());
                out_file << setw(2) << atoi(temp_str.substr(4,1).c_str());
            }
            out_file << endl;
        }
        out_file << endl;
    }
    
    out_file << setprecision(3) << setfill(' ');
    //Prints the column (X) dimensions of the simulation to the output file
    for(int i = 0; i < num_cols; i++) {
        in_file >> temp_real;
        out_file << " " << temp_real;
    }
    out_file << endl;
    
    //Prints the row (Y) dimensions of the simulation to the output file
    for(int i = 0; i < num_rows; i++) {
        in_file >> temp_real;
        out_file << " " << temp_real;
    }
    out_file << endl;
    
    //Prints the slice (Z) dimensions of the simulation to the output file
    out_file << " " << 1.0 << endl;
    
    //Heat Production Values
    out_file << " " << 8;
    for(int i = 0; i < 8; i++) {
        in_file >> temp_real;
        out_file << " " << scientific << temp_real;
    }
    out_file << endl;
    
    //Thermal Conductivity Difference Values
    out_file << " " << 8;
    for(int i = 0; i < 8; i++) {
        in_file >> temp_real;
        out_file << " " << scientific << temp_real;
    }
    out_file << endl;
    
    if(using_convection) {
        //Fluid Heat Capacity values
        out_file << " " << 8;
        for(int i = 0; i < 8; i++) {
            in_file >> temp_real;
            out_file << " " << scientific << temp_real;
        }
        out_file << endl;
        
        //Rock Heat Capacity values
        out_file << " " << 8;
        for(int i = 0; i < 8; i++) {
            in_file >> temp_real;
            out_file << " " << scientific << temp_real;
        }
        out_file << endl;
        
        //Minimum Convection Temps
        out_file << " " << 8;
        for(int i = 0; i < 8; i++) {
            in_file >> temp_real;
            out_file << " " << scientific << temp_real;
        }
        out_file << endl;
        
        //Convection Velocity Values
        out_file << " " << 8;
        for(int i = 0; i < 8; i++) {
            in_file >> temp_real;
            out_file << " " << scientific << temp_real;
        }
        out_file << endl;
    }
    
    //Close the input files
    in_file.close();
    out_file.close();
}
