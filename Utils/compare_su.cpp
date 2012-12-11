/**
 * Compares two DSAA surfer grid files.
 * Reports any differences above a user specified tolerance and displays
 * a summary of the results of the comparision. All cli output is also
 * written to a log file.
 */
 
#ifdef _WIN32
#define NOMINMAX //FYI need to disable min/max macro in windows.h
#include <windows.h>
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <limits>
#include <math.h>

#define REAL double

using std::cerr;
using std::cin;
using std::cout;
using std::endl;
using std::string;
using std::ofstream;
using std::ifstream;
using std::setw;
using std::fixed;
using std::scientific;
using std::setprecision;
using std::ios;
using std::numeric_limits;
using std::streamsize;
using std::max;
using std::flush;

/*
 * Clears the cin buffer
 */
void clear_cin() {
	cin.clear();
	cin.ignore(numeric_limits <streamsize> ::max(), '\n' );
}

int main(int argc, char **argv) {
    string input;                                     //Input buffer
    string in1_filename;                              //The first surfer filename
    string in2_filename;                              //The second surfer filename
    ifstream in1, in2;                                //Input file streams
    ofstream out;                                     //Output file stream
    int dim_x[2], dim_y[2];                           //X and Y dimensions of the two files
    REAL min_x[2],min_y[2],max_x[2],max_y[2];         //Minimum and Maximum X and Y values
    REAL min_temp[2],max_temp[2],temp[2];             //Minimum, Maximum, and temporary temperature variables
    REAL max_diff = 0.0;                              //Maximum difference found
    REAL sum = 0.0;                                   //Sum of the differences found
    REAL diff = 0.0;                                  //Absolute temperature difference
    REAL tolerance = 0.01;                            //Defines how much of a differnce in temperature is required for the temperatures to be considered different
    
    //Opens the log file for writting
    out.open("log.txt",ios::out);
    
    //Asks the user for the first surfer filename
    do {
        cout << "First Surfer File Name: ";
        cin >> in1_filename;
        in1.open(in1_filename.c_str(),ios::in);
        if(!in1.is_open()) {
            cout << "File Not found!" << endl;
        }
    } while(!in1.is_open());
    
    //Asks the user for the second surfer filename
    do {
        cout << "Second Surfer File Name: ";
        cin >> in2_filename;
        in2.open(in2_filename.c_str(),ios::in);
        if(!in2.is_open()) {
            cout << "File Not found!" << endl;
        }
    } while(!in2.is_open());
    
    //Asks the user for teh tolerance
    cout << "Enter the minimum value to be considered different: ";
    while(!(cin >> tolerance)) {
		clear_cin();
        cout << "Incorrect input, enter a number greater than 0: ";
    }
    
    
    cout << "Opening \"" << in1_filename << "\" and \"" << in2_filename << "\" for comparision" << endl;
    cout << "Tolerance = " << tolerance << endl;
    out << "Opening \"" << in1_filename << "\" and \"" << in2_filename << "\" for comparision" << endl;
    out << "Tolerance = " << tolerance << endl;
    
    //Checks if file 1 is a DSAA surfer grid file
    in1 >> input;
    if(input.compare("DSAA") != 0) {
        cout << "File 1 is not a DSAA surfer grid file" << endl;
        out << "File 1 is not a DSAA surfer grid file" << endl;
        in1.close();
        out.close();
        exit(1);
    }
    
    //Checks if file 2 is a DSAA surfer grid file
    in2 >> input;
    if(input.compare("DSAA") != 0) {
        cout << "File 2 is not a DSAA surfer grid file" << endl;
        out << "File 2 is not a DSAA surfer grid file" << endl;
        in1.close();
        in2.close();
        out.close();
        exit(1);
    }
    
    //Retrieves and displays the specifications of the two files
    in1 >> dim_y[0] >> dim_x[0] >> min_x[0] >> max_x[0] >> min_y[0] >> max_y[0] >> min_temp[0] >> max_temp[0];
    in2 >> dim_y[1] >> dim_x[1] >> min_x[1] >> max_x[1] >> min_y[1] >> max_y[1] >> min_temp[1] >> max_temp[1];
    
    cout << dim_x[0] << " " << dim_y[0] << " " << min_x[0] << " " << max_x[0] << " " << min_y[0] << " " << max_y[0] << " " << min_temp[0] << " " << max_temp[0] << endl;
    cout << dim_x[1] << " " << dim_y[1] << " " << min_x[1] << " " << max_x[1] << " " << min_y[1] << " " << max_y[1] << " " << min_temp[1] << " " << max_temp[1] << endl;
    out << dim_x[0] << " " << dim_y[0] << " " << min_x[0] << " " << max_x[0] << " " << min_y[0] << " " << max_y[0] << " " << min_temp[0] << " " << max_temp[0] << endl;
    out << dim_x[1] << " " << dim_y[1] << " " << min_x[1] << " " << max_x[1] << " " << min_y[1] << " " << max_y[1] << " " << min_temp[1] << " " << max_temp[1] << endl;
    
    //Checks if the dimensions of the two files match
    if((dim_x[0] != dim_x[1]) || (dim_y[0] != dim_y[1])) {
        cout << "file specifications do not match, files are different" << endl;
        out << "file specifications do not match, files are different" << endl;
    }
    else {
        //List any differences in the other file specifications
        if(max_x[0] != max_x[1]) {
            cout << "X dimensions do not match!" << endl;
            out << "X dimensions do not match!" << endl;
        }
        if(min_y[0] != min_y[1]) {
            cout << "Y dimensions do not match!" << endl;
            out << "Y dimensions do not match!" << endl;
        }
        if(min_temp[0] != min_temp[1]) {
            cout << "Minimum Temperatures do not match!" << endl;
            out << "Minimum Temperatures do not match!" << endl;
        }
        if(max_temp[0] != max_temp[1]) {
            cout << "Maximum Temperatures do not match!" << endl;
            out << "Maximum Temperatures do not match!" << endl;
        }
        cout <<endl;
        out << endl;
        
        //Compares the temperatures of the two files
        int count = 0;
        for(int i = 0; i < dim_x[0]; i++) {
            for(int j = 0; j < dim_y[0]; j++) {
                in1 >> temp[0];
                in2 >> temp[1];
                diff = fabs(temp[0]-temp[1]);
                if(diff >= tolerance) {
                    if(diff > max_diff) {
                        max_diff = diff;
                    }
                    cout << "Diff found at x=" << setw(3) << i << ", y=" << setw(3) << j << "; Temps = " << setw(8) << fixed << setprecision(2) << temp[0] << ":" << setw(8) << temp[1] << "; Diff = " << setw(15) << setprecision(10) << diff << endl;
                    out << "Diff found at x=" << setw(3) << i << ", y=" << setw(3) << j << "; Temps = " << setw(8) << fixed << setprecision(2) << temp[0] << ":" << setw(8) << temp[1] << "; Diff = " << setw(15) << setprecision(10) << diff << endl;
                    count++;
                    sum += diff;
                }
            }
        }
        
        //Displays a summary of the comparision of the two files
        if(count != 0) {
            cout << endl << endl << "Tolerance                   = " << tolerance << endl;
            cout << "Number of Differences Found = " << count << endl;
            cout << "Maximum Difference Found    = " << max_diff << endl;
            cout << "Average Difference Found    = " << sum/count << endl;
            
            out << endl << endl << "Tolerance                   = " << tolerance << endl;
            out << "Number of Differences Found = " << count << endl;
            out << "Maximum Difference Found    = " << max_diff << endl;
            out << "Average Difference Found    = " << sum/count << endl;
        }
        else {
            cout << endl << endl << "Files are identical, No differences found" << endl;
            out << endl << endl << "Files are identical, No differences found" << endl;
        }
    }
    
    //Closes the input and output files
    in1.close();
    in2.close();
    out.close();
}
