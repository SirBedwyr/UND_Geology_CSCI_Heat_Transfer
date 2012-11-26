/**
 * Compares two DSAA surfer grid files.
 * Reports any differences above a user specified tolerance and displays
 * a summary of the results of the comparision. All cli output is also
 * written to a log file.
 * The program takes in three command line arguments:
 * <filename 1> <filename 2> <tolerance>
 *         filename 1 - The name of the first surfer file
 *         filename 2 - The name of the second surfer file
 *        tolerance  - The minimum temperature difference to be considered different
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

#define REAL float

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

int main(int argc, char **argv) {
    string input;                                     //Input buffer
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
    
    //Checks if the correct number of command line arguments was provided
    if(argc != 4) {
        cout << "Incorrect command line arguments" << endl;
        cout << "compare <filename 1> <filename 2> <tolerance>" << endl;
        out << "Incorrect command line arguments" << endl;
        out << "compare <filename 1> <filename 2> <tolerance>" << endl;
        exit(1);
    }
    
    tolerance = atof(argv[3]);    //Retrieves the user specified tolerance
    
    cout << "Opening \"" << argv[1] << "\" and \"" << argv[2] << "\" for comparision" << endl;
    cout << "Tolerance = " << tolerance << endl;
    out << "Opening \"" << argv[1] << "\" and \"" << argv[2] << "\" for comparision" << endl;
    out << "Tolerance = " << tolerance << endl;
    
    //Checks if file 1 exists and is a DSAA surfer grid file
    in1.open(argv[1],ios::in);
    if(!in1.is_open()) {
        cout << "Unable to open surfer file 1" << endl;
        out << "Unable to open surfer file 1" << endl;
        out.close();
        exit(1);
    }
    in1 >> input;
    if(input.compare("DSAA") != 0) {
        cout << "File 1 is not a DSAA surfer grid file" << endl;
        out << "File 1 is not a DSAA surfer grid file" << endl;
        in1.close();
        out.close();
        exit(1);
    }
    
    //Checks if file 1 exists and is a DSAA surfer grid file
    in2.open(argv[2],ios::in);
    if(!in2.is_open()) {
        cout << "Unable to open surfer file 2" << endl;
        out << "Unable to open surfer file 2" << endl;
        in1.close();
        out.close();
        exit(1);
    }
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
