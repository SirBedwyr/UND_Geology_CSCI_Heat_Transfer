/**
 * Performs a finite difference heat flow 
 * simulation using conduction and convection.
 */

#ifdef _WIN32
#define _USE_MATH_DEFINES
#define NOMINMAX //FYI need to disable min/max macro in windows.h
#include <windows.h>
#endif

#ifdef DISPLAY
#ifdef __APPLE__
#	include <OpenGL/gl.h>
#	include <OpenGL/glu.h>
#	include <GLUT/glut.h>
#else
#	include <GL/GL.h>
#	include <GL/GLU.h>
#	include <GL/glut.h>
#endif
#	include "GL/freeglut.h"
#endif

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <limits>
#include <cmath>
#include <string.h>
#include <mpi.h>

#ifdef TIME
#ifdef _WIN32
#include <time.h>
SYSTEMTIME start, end;
#else
#include<time.h>		//Time library
#include<sys/time.h>	//Used for the gettimeofday function and timeval struct
struct timeval start, end;
#endif
#endif

/**
 * Time_step (seconds/yr) divided by the product of density and heat capacity.
 * The value for the density is 2200 kg/m^3 and heat capacity is 1000 kJ/kg K.
 */
#define QFAC 14.33                 //Description is defined in the previous comment
#define DTC 0.25                   //
#define OUT_PRECISION 10           //Number of digits to print after the decimal place for floating point values
#define INDEX_WIDTH 2              //The number of characters to print and read for each conduction and convection code
#define REAL double                //The precision of the model.
#define MPI_REAL_VALUE MPI_DOUBLE  //The data type of the MPI messages

using std::cerr;
using std::cin;
using std::cout;
using std::endl;
using std::string;
using std::ofstream;
using std::ifstream;
using std::ostringstream;
using std::stringstream;
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

//Function prototypes
void save_surfer();
void save_model_state();
void conduction();
void convection();
void PressEnterToContinue();
REAL find_max_temp_diff();
void update_moving_sources();
void find_loc_index(REAL x_loc, REAL y_loc, REAL z_loc, int *index);


//MPI Requests for asynchronous sends and receives
MPI_Request *request1;
MPI_Request *request2;
MPI_Request *request3;
MPI_Request *request4;
REAL temp_time_step;

//Conduction code specific variables
int ***cond_codes;              //The unmodified conduction codes as read from the input file
int ***cond_hp_index;           //The conduction index for the radioactive heat production array
int ***cond_tc_index;           //The conduction index for the thermal conductivity array
int ***use_cond;                //Flag to indicate if conduction occurs for a given cell
REAL DHF;                       //


 //Convection code specific variables
int ***conv_codes;              //The convection codes as read from the input file
int ***conv_min_temp_index;     //The convection index for the minimum temp for convection array
int ***conv_direction;          //The direction of convection following the direction matrix in the previous comment
int ***conv_vel_index;          //The convection index for the velocity array
int ***conv_fluid_index;        //The convection index for the fluid heat capacity array
int ***conv_rock_index;         //The convection index for the rock heat capacity array
int num_conv_loops;             //The number of convection updates to perform per time step
REAL time_inc;                  //The amount of time increment per convection loop

//File names
string source_filename;         //The input files name with extension
string output_filename;         //The output state files name with extension
string output_su_filename;      //The output surfer files name with extension

//Input file variables
string title;                   //The title of the input file
int using_convection = -1;      //Indicates if convection is being used
REAL ***temp;                   //The current temperature array
int num_rows;                   //The number of rows for the simulation
int num_cols;                   //The number of columns for the simulation
int num_slices;                 //Total number of slices to form the 3d simulation (one 'slice' has dimension rows x columns)
REAL *dim_x;                    //The dimensions of each column in the x direction
REAL *dim_y;                    //The dimensions of each row in the y direction
REAL *dim_z;                    //The dimensions of each row in the z direction
REAL *dist_x;                   //The distance from the origin to the center of a column for a given column index in the x direction
REAL *dist_y;                   //The distance from the origin to the center of a row for a given row index in the y direction
REAL *dist_z;                   //The distance from the origin to the center of a slice for a given slice index in the z direction
REAL max_dist_x;                //The maximum x distance
REAL max_dist_y;                //The maximum y distance
REAL max_dist_z;                //The maximum z distance
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

//Moving source variables
/*int using_moving_source = -1;   //Indicates if a moving source is being used
int mvsrc_start_col = -1;       //Starting column of the moving source
int mvsrc_end_col = -1;         //Ending column of the moving source
int mvsrc_depth = -1;           //The depth of the moving source (defines the row in the simulation)
*/
int using_moving_source = -1;   //Indicates if a moving source is being used
int num_mvsrc = -1;                  //The number of moving sources
REAL *mvsrc_x;                  //The x location of the moving sources
REAL *mvsrc_y;                  //The y location of the moving sources
REAL *mvsrc_z;                  //The z location of the moving sources
REAL *mvsrc_offset_x;           //The size of the moving sources in the x direction
REAL *mvsrc_offset_y;           //The size of the moving sources in the y direction
REAL *mvsrc_offset_z;           //The size of the moving sources in the z direction
REAL *mvsrc_vel_x;              //The x component of the moving sources velocity vectors
REAL *mvsrc_vel_y;              //The y component of the moving sources velocity vectors
REAL *mvsrc_vel_z;              //The z component of the moving sources velocity vectors
REAL *mvsrc_accel_x;            //The x component of the moving sources acceleration vectors
REAL *mvsrc_accel_y;            //The y component of the moving sources acceleration vectors
REAL *mvsrc_accel_z;            //The z component of the moving sources acceleration vectors
REAL *mvsrc_temp;               //The temperature of the moving sources
int *mvsrc_valid;               //Indicates if a moving source is valid

//Time specific variables
REAL sim_time;                  //The current time of the simulation
REAL tic;                       //Time variable used in convection calculations
REAL time_step = -1;            //The amount of time that passes between each update of the simulation, time step
REAL run_time = -1;             //The run time of the simulation

//Global variables
int save_state = 0;            //Indicates if the model should save the current state at each screen update
int save_result = 0;           //Indicates if the model should save the final result of the simulaiton
int use_tolerance = 0;         //Indicates if the model should stop once a user specified tolerance is met for temperature change
REAL max_vel;                   //The maximum convection velocity of the velocity array
REAL min_row_dim;               //The minimum y dimension of each cell of the simulation
REAL min_col_dim;               //The minimum x dimension of each cell of the simulation
REAL min_slice_dim;             // The minimum z dimension of each cell in the simulation
REAL thermal_time_constant;     //The thermal time constant of the model, used in the selection of the run time
REAL ***next_temp;              //The next temperature array
REAL max_thermal_conduct_diff;  //The maximum thermal conductivity difference
REAL min_thermal_conduct_diff;  //The minimum thermal conductivity difference
int num_loops = -1;             //The number of loops between screen updates
int su_num_width;               //The number of characters for the slice number in the output surfer filenames
unsigned long long count = 0;	//The current loop
REAL tolerance = -1;            //The maximum difference required for the model to stop
REAL max_temp_diff;             //The maximum temperature difference between the current and next temperature arrays
int my_rank;                    //The rank of the current process
int comm_sz;                    //The total number of processes in the simulation
int *chunk_size;                //The number of rows per MPI process
int *displacements;             //The row displacements for each process
int st_index = 0;               //Starting row index for a each mpi process
int en_index = 0;               //Ending row index for each mpi process

/**
 * Deallocates all allocated memory used by the program
 */
void deallocate_memory() {
     //Deletes allocated memory
    for(int i = 0; i < chunk_size[my_rank]+2; i++) {
        for(int j = 0; j < num_cols; j++) {
            delete[] temp[i][j];
            delete[] next_temp[i][j];
            delete[] cond_codes[i][j];
            delete[] cond_hp_index[i][j];
            delete[] cond_tc_index[i][j];
            delete[] use_cond[i][j];
            if(using_convection) {
                delete[] conv_codes[i][j];
                delete[] conv_min_temp_index[i][j];
                delete[] conv_direction[i][j];
                delete[] conv_vel_index[i][j];
                delete[] conv_fluid_index[i][j];
                delete[] conv_rock_index[i][j];
            }
        }
        delete[] temp[i];
        delete[] next_temp[i];
        delete[] cond_codes[i];
        delete[] cond_hp_index[i];
        delete[] cond_tc_index[i];
        delete[] use_cond[i];
        if(using_convection) {
            delete[] conv_codes[i];
            delete[] conv_min_temp_index[i];
            delete[] conv_direction[i];
            delete[] conv_vel_index[i];
            delete[] conv_fluid_index[i];
            delete[] conv_rock_index[i];
        }
    }
    delete[] dim_x;
    delete[] dim_y;
    delete[] dim_z;
    delete[] dist_x;
    delete[] dist_y;
    delete[] dist_z;
    delete[] temp;
    delete[] next_temp;
    delete[] cond_codes;
    delete[] cond_hp_index;
    delete[] cond_tc_index;
    delete[] use_cond;
    delete[] heat_production_values;
    delete[] thermal_conduct_diff;
    if(using_convection == 1) {
        delete[] conv_codes;
        delete[] conv_min_temp_index;
        delete[] conv_direction;
        delete[] conv_vel_index;
        delete[] conv_fluid_index;
        delete[] conv_rock_index;
        delete[] heat_capac_fluid;
        delete[] heat_capac_rock;
        delete[] min_temp_conv;
        delete[] vel;
    }
    if(using_moving_source == 1) {
        delete[] mvsrc_x;
        delete[] mvsrc_y;
        delete[] mvsrc_z;
        delete[] mvsrc_offset_x;
        delete[] mvsrc_offset_y;
        delete[] mvsrc_offset_z;
        delete[] mvsrc_vel_x;
        delete[] mvsrc_vel_y;
        delete[] mvsrc_vel_z;
        delete[] mvsrc_accel_x;
        delete[] mvsrc_accel_y;
        delete[] mvsrc_accel_z;
        delete[] mvsrc_temp;
        delete[] mvsrc_valid;
    }
}

#ifdef DISPLAY
//Display variables
int display_mode = -1;
int array_size;
REAL min_temp;
REAL max_temp;
REAL layer_min_temp;
REAL layer_max_temp;
float *color_field;
float transparency = 1.0f;
int current_slice = 0;

/**
 *  Parameters to control the camera angle so we can move where we're looking at
 *  the simulation from with the mouse. Kept for debugging.
 */
/*
int     ox                  = 0;
int     oy                  = 0;
int     buttonState         = 0; 
float   camera_trans[]      = {0, -0.2, -10};
float   camera_rot[]        = {0, 0, 0};
float   camera_trans_lag[]  = {0, -0.2, -10};
float   camera_rot_lag[]    = {0, 0, 0};
const float inertia         = 0.1f;
*/

/*
 * Sets max and min temperature values for simulation
 */
void array_minmax() {
	min_temp=temp[0][0][0];
	max_temp=temp[0][0][0];
    for(int k=0; k<num_slices; k++) {
		for (int i=0; i<num_rows; i++) {
			for(int j=0; j<num_cols; j++) {
				if(temp[i][j][k]<min_temp)
					min_temp = temp[i][j][k];
				if(temp[i][j][k]>max_temp)
					max_temp = temp[i][j][k];
			}
		}
	}
}

/**
 * Updates max_temp with the maximum temp of
 * the current slice.
 */
void array_max() {
    max_temp=temp[0][0][0];
    for(int i = 0; i < num_rows; i++) {
        for(int j = 0; j < num_cols; j++) {
            if(temp[i][j][current_slice] > max_temp) {
                max_temp = temp[i][j][current_slice];
            }
        }
    }
}

/* 
 * 3D to 1D indexing
 */
static int POSITION(int x, int y) {
	return (x*num_cols)+y;
}


/*
 * Colormap algorithm reproduces Matlab's RGB "Jet" plate
 * Concept based on: http://paulbourke.net/texture_colour/colourspace/ (11/21/12)
 */
void jet_color_set(int x, int y, int z) {
	REAL current_temp = temp[x][y][z];
	REAL delta_temp = max_temp - min_temp;
	
	if(current_temp < min_temp)
		current_temp = min_temp;
	if(current_temp > max_temp)
		current_temp = max_temp;

	if(current_temp < (min_temp + 0.25 * delta_temp)) {
		color_field[POSITION(x,y) * 3] = (GLfloat)0.0;	
		color_field[POSITION(x,y) * 3 + 1] = (GLfloat)(4*(current_temp - min_temp)/delta_temp);
		color_field[POSITION(x,y) * 3 + 2] = (GLfloat)1.0;
	}
	else if(current_temp < (min_temp + 0.5 * delta_temp)) {
		color_field[POSITION(x,y) * 3] = (GLfloat)0.0;	
		color_field[POSITION(x,y) * 3 + 1] = (GLfloat)1.0;	
		color_field[POSITION(x,y) * 3 + 2] = (GLfloat)(1.0 + 4 * (min_temp + 0.25 * delta_temp - current_temp) / delta_temp);
	}
	else if(current_temp < (min_temp + 0.75 * delta_temp)) {
		color_field[POSITION(x,y) * 3] = (GLfloat)(4 * (current_temp - min_temp - 0.5 * delta_temp) / delta_temp);	
		color_field[POSITION(x,y) * 3 + 1] = (GLfloat)1.0;
		color_field[POSITION(x,y) * 3 + 2] = (GLfloat)0.0;
	}
	else {
		color_field[POSITION(x,y) * 3] = (GLfloat)1.0;	
		color_field[POSITION(x,y) * 3 + 1] = (GLfloat)(1.0 + 4 * (min_temp + 0.75 * delta_temp - current_temp) / delta_temp);
		color_field[POSITION(x,y) * 3 + 2] = (GLfloat)0.0;
	}
}


/*
 * Draw temp surface via 1x1 faces. Dimensions constant for simplicity.
 * Originally drew cubes from Robert Bergmans voxel display code.
 */
void draw_cube(int x, int y, int z) {
	if(z == current_slice) {
		transparency = 1.0f;
	}
	else {
		transparency = 0.3f;
	}
	glBegin(GL_TRIANGLES);
		//front
		glColor4f(	color_field[POSITION(x,y) * 3],
					color_field[POSITION(x,y) * 3 + 1],
					color_field[POSITION(x,y) * 3 + 2],
					transparency);
		glVertex3f(0.0f,0.0f,1.0f);//5
		glVertex3f(1.0f,0.0f,1.0f);//6
		glVertex3f(0.0f,-1.0f,1.0f);//8
		
		glColor4f(	color_field[POSITION(x,y) * 3],
					color_field[POSITION(x,y) * 3 + 1],
					color_field[POSITION(x,y) * 3 + 2],
					transparency);
		glVertex3f(0.0f,-1.0f,1.0f);//8
		glVertex3f(1.0f,0.0f,1.0f);//6
		glVertex3f(1.0f,-1.0f,1.0f);//7

		//top
		glColor4f(	color_field[POSITION(x,y) * 3],
					color_field[POSITION(x,y) * 3 + 1],
					color_field[POSITION(x,y) * 3 + 2],
					transparency);
		glVertex3f(0.0f,0.0f,0.0f);//1
		glVertex3f(1.0f,0.0f,0.0f);//2
		glVertex3f(0.0f,0.0f,1.0f);//5
		
		glColor4f(	color_field[POSITION(x,y) * 3],
					color_field[POSITION(x,y) * 3 + 1],
					color_field[POSITION(x,y) * 3 + 2],
					transparency);
		glVertex3f(0.0f,0.0f,1.0f);//5
		glVertex3f(1.0f,0.0f,0.0f);//2
		glVertex3f(1.0f,0.0f,1.0f);//6
	
	
	//QUADS code left in case we need it later	
	/*
	glBegin(GL_QUADS);
		//front
		glColor4f(	color_field[POSITION(x,y) * 3],
					color_field[POSITION(x,y) * 3 + 1],
					color_field[POSITION(x,y) * 3 + 2],
					transparency);
		glVertex3f(0.0f,0.0f,1.0f);//5
		glVertex3f(1.0f,0.0f,1.0f);//6
		glVertex3f(1.0f,-1.0f,1.0f);//7
		glVertex3f(0.0f,-1.0f,1.0f);//8

		//top
		glColor4f(	color_field[POSITION(x,y) * 3],
					color_field[POSITION(x,y) * 3 + 1],
					color_field[POSITION(x,y) * 3 + 2],
					transparency);
		glVertex3f(0.0f,0.0f,0.0f);//1
		glVertex3f(1.0f,0.0f,0.0f);//2
		glVertex3f(1.0f,0.0f,1.0f);//6
		glVertex3f(0.0f,0.0f,1.0f);//5

		*//*//left
		glColor4f(	color_field[POSITION(x,y) * 3],
					color_field[POSITION(x,y) * 3 + 1],
					color_field[POSITION(x,y) * 3 + 2],
					transparency);
		glVertex3f(0.0f,0.0f,0.0f);//1
		glVertex3f(0.0f,0.0f,1.0f);//5
		glVertex3f(0.0f,-1.0f,1.0f);//8
		glVertex3f(0.0f,-1.0f,0.0f);//4

		//right
		glColor4f(	color_field[POSITION(x,y) * 3],
					color_field[POSITION(x,y) * 3 + 1],
					color_field[POSITION(x,y) * 3 + 2],
					transparency);
		glVertex3f(1.0f,0.0f,0.0f);//2
		glVertex3f(1.0f,0.0f,1.0f);//6
		glVertex3f(1.0f,-1.0f,1.0f);//7
		glVertex3f(1.0f,-1.0f,0.0f);//3

		//bottom
		glColor4f(	color_field[POSITION(x,y) * 3],
					color_field[POSITION(x,y) * 3 + 1],
					color_field[POSITION(x,y) * 3 + 2],
					transparency);
		glVertex3f(0.0f,-1.0f,0.0f);//4
		glVertex3f(1.0f,-1.0f,0.0f);//3
		glVertex3f(1.0f,-1.0f,1.0f);//7
		glVertex3f(0.0f,-1.0f,1.0f);//8

		//back
		glColor4f(	color_field[POSITION(x,y) * 3],
					color_field[POSITION(x,y) * 3 + 1],
					color_field[POSITION(x,y) * 3 + 2],
					transparency);
		glVertex3f(0.0f,0.0f,0.0f);//1
		glVertex3f(1.0f,0.0f,0.0f);//2
		glVertex3f(1.0f,-1.0f,0.0f);//3
		glVertex3f(0.0f,-1.0f,0.0f);//4
		*/

	glEnd();
}


/*
 *Draw all HUD/Overlay graphics. Quads for temp scale are hardcoded to Jet color map. Any new color map will require changes.
 */
void displayOverlay(){	
	int windowWidth = glutGet(GLUT_WINDOW_WIDTH);
	int windowHeight = glutGet(GLUT_WINDOW_HEIGHT);
    
	glMatrixMode( GL_PROJECTION );
    glPushMatrix();
        glLoadIdentity();
        glOrtho(0.0f,windowWidth,windowHeight,0.0f,0.0f,1.0f);

        glMatrixMode( GL_MODELVIEW );
        glPushMatrix();
            glLoadIdentity();
			glBegin( GL_QUADS );

				glColor3f( 0.0f, 0.0f, 1.0f );
                glVertex2f( (GLfloat)(windowWidth/2-75), 20.0f );
                glVertex2f( (GLfloat)(windowWidth/2-45), 20.0f );
                glVertex2f( (GLfloat)(windowWidth/2-45), 50.0f );
                glVertex2f( (GLfloat)(windowWidth/2-75), 50.0f );

				glColor3f( 0.0f, 1.0f, 1.0f );
                glVertex2f( (GLfloat)(windowWidth/2-45), 20.0f );
                glVertex2f( (GLfloat)(windowWidth/2-15), 20.0f );
                glVertex2f( (GLfloat)(windowWidth/2-15), 50.0f );
                glVertex2f( (GLfloat)(windowWidth/2-45), 50.0f );

				glColor3f( 0.0f, 1.0f, 0.0f );
                glVertex2f( (GLfloat)(windowWidth/2-15), 20.0f );
                glVertex2f( (GLfloat)(windowWidth/2+15), 20.0f );
                glVertex2f( (GLfloat)(windowWidth/2+15), 50.0f );
                glVertex2f( (GLfloat)(windowWidth/2-15), 50.0f );

				glColor3f( 1.0f, 1.0f, 0.0f );
                glVertex2f( (GLfloat)(windowWidth/2+15), 20.0f );
                glVertex2f( (GLfloat)(windowWidth/2+45), 20.0f );
                glVertex2f( (GLfloat)(windowWidth/2+45), 50.0f );
                glVertex2f( (GLfloat)(windowWidth/2+15), 50.0f );

				glColor3f( 1.0f, 0.0f, 0.0f );
                glVertex2f( (GLfloat)(windowWidth/2+45), 20.0f );
                glVertex2f( (GLfloat)(windowWidth/2+75), 20.0f );
                glVertex2f( (GLfloat)(windowWidth/2+75), 50.0f );
                glVertex2f( (GLfloat)(windowWidth/2+45), 50.0f );

            glEnd();
        glPopMatrix();
		glPushMatrix();
			glLoadIdentity();
			ostringstream str1;

			str1 << "Min Temp <                                          > Max Temp";
			glColor3f(1.0f, 1.0f, 1.0f); 
			glRasterPos2f((GLfloat)(windowWidth/2-150),35.0f);
			glutBitmapString(GLUT_BITMAP_HELVETICA_12, (const unsigned char *)str1.str().c_str());
			str1.str("");
			str1.clear();
			
			str1 << setw(4) << min_temp;
			glRasterPos2f((GLfloat)(windowWidth/2-100), 65.0f);
			glutBitmapString(GLUT_BITMAP_HELVETICA_12, (const unsigned char *)str1.str().c_str());
			str1.str("");
			str1.clear();

			str1 << setw(4) << (max_temp+min_temp)/2;
			glRasterPos2f((GLfloat)(windowWidth/2-20), 65.0f);
			glutBitmapString(GLUT_BITMAP_HELVETICA_12, (const unsigned char *)str1.str().c_str());
			str1.str("");
			str1.clear();


			str1 << setw(4) << max_temp;
			glRasterPos2f((GLfloat)(windowWidth/2+60), 65.0f);
			glutBitmapString(GLUT_BITMAP_HELVETICA_12, (const unsigned char *)str1.str().c_str());
			str1.str("");
			str1.clear();


		glPopMatrix();
        glPushMatrix();
            glLoadIdentity();
            ostringstream str;

			str << "Time Interval:" << endl;
			if(using_convection) {

				str << "Conv. Time Interval:" << endl;
				str << "Num Conv. Loops:" << endl;
			}

			str << endl << "Loop Total:" << endl;
			str << "Sim Time:" << endl;
			str << "Cum. Sim Time:" << endl << endl;
		
			str << "Model Dimensions:" << endl;
			str << "Current Slice:" << endl;
			str << "CHF:" << endl;
			
			glColor3f(1.0f, 1.0f, 1.0f); 
			glRasterPos2f(10.0f,(GLfloat)(windowHeight*3.0/4.0));
			glutBitmapString(GLUT_BITMAP_HELVETICA_12, (const unsigned char *)str.str().c_str());
		

			str.str("");
			str.clear();
		
			str << setprecision(4) << scientific;
			str << time_step << endl;
			if(using_convection) {
				str << time_inc << endl;
				str << num_conv_loops << endl;
			}

			str << endl << count << endl;
			str << sim_time << endl;
			str << initial_time + sim_time << endl << endl;
			
			str << num_rows << " X " << num_cols << " X " << num_slices << endl;
			str << current_slice+1 << endl;
			str << chf << endl;
			
			glColor3f(1.0f, 1.0f, 1.0f); 
			glRasterPos2f(150.0f,(GLfloat)(windowHeight*3.0/4.0));
			glutBitmapString(GLUT_BITMAP_HELVETICA_12, (const unsigned char *)str.str().c_str());
		            
        glPopMatrix();

    glMatrixMode( GL_PROJECTION );
    glPopMatrix();
}

/*
 * Helper function called from display3d. Broken out for readability.
 */
void display_helper() {
    array_max();
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_DEPTH_TEST);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	//  Handle the camera angle. Maintained for debugging
	/*
    for (int c = 0; c < 3; ++c)
    {
        camera_trans_lag[c] += (camera_trans[c] - camera_trans_lag[c]) * inertia;
        camera_rot_lag[c] += (camera_rot[c] - camera_rot_lag[c]) * inertia;
    }

    glTranslatef(camera_trans_lag[0], camera_trans_lag[1], camera_trans_lag[2]);
    glRotatef(camera_rot_lag[0], 1.0, 0.0, 0.0);
    glRotatef(camera_rot_lag[1], 0.0, 1.0, 0.0);
	*/

	//Draw the boundary lines
	glBegin(GL_LINES);
		glColor3f(1.0, 1.0, 1.0);

		glVertex3f(0.0f,0.0f,0.0f);
		glVertex3f((GLfloat)num_cols,0.0f,0.0f);
		glVertex3f((GLfloat)num_cols,0.0f,0.0f);
		glVertex3f((GLfloat)num_cols,(GLfloat)-num_rows,0.0f);
		glVertex3f((GLfloat)num_cols,(GLfloat)-num_rows,0.0f);
		glVertex3f(0.0f,(GLfloat)-num_rows,0.0f);
		glVertex3f(0.0f,(GLfloat)-num_rows,0.0f);
		glVertex3f(0.0f,0.0f,0.0f);

		glVertex3f(0.0f,0.0f,(GLfloat)num_slices);
		glVertex3f((GLfloat)num_cols,0.0f,(GLfloat)num_slices);
		glVertex3f((GLfloat)num_cols,0.0f,(GLfloat)num_slices);
		glVertex3f((GLfloat)num_cols,(GLfloat)-num_rows,(GLfloat)num_slices);
		glVertex3f((GLfloat)num_cols,(GLfloat)-num_rows,(GLfloat)num_slices);
		glVertex3f(0.0f,(GLfloat)-num_rows,(GLfloat)num_slices);
		glVertex3f(0.0f,(GLfloat)-num_rows,(GLfloat)num_slices);
		glVertex3f(0.0f,0.0f,(GLfloat)num_slices);

		glVertex3f(0.0f,0.0f,0.0f);
		glVertex3f(0.0f,0.0f,(GLfloat)num_slices);
		glVertex3f((GLfloat)num_cols,0.0f,0.0f);
		glVertex3f((GLfloat)num_cols,0.0f,(GLfloat)num_slices);
		glVertex3f((GLfloat)num_cols,(GLfloat)-num_rows,0.0f);
		glVertex3f((GLfloat)num_cols,(GLfloat)-num_rows,(GLfloat)num_slices);
		glVertex3f(0.0f,(GLfloat)-num_rows,0.0f);
		glVertex3f(0.0f,(GLfloat)-num_rows,(GLfloat)num_slices);
	glEnd();
	
	for(int i=0; i<num_rows; i++) {
		for (int j=0; j<num_cols; j++) {
			glPushMatrix();
			glTranslatef((GLfloat)j,(GLfloat)-i,(GLfloat)current_slice);
			jet_color_set(i,j,current_slice);
			draw_cube(i,j,current_slice);
			glPopMatrix();
		}
	}

	displayOverlay();
	glutSwapBuffers();
	glutPostRedisplay();
}

/*
 * Main simulation call during display. Makes required simulation calls e.g.- for convection, etc. 
 * Makes call to display helper function
 */
void display3D() {
	//Displays status information for the current loop
	if(count%num_loops == 0) {
        if(use_tolerance == 0) {
		    cout << setw(15) << count << setw(20) << fixed << setprecision(5) << sim_time << setw(20) << initial_time + sim_time << endl;
        }
        else {
            cout << setw(15) << count << setw(20) << fixed << setprecision(5) << sim_time << setw(20) << initial_time + sim_time << setw(20) << max_temp_diff << endl;
        }

		//Saves the current state of the simulation if the save_state flag is set
		if(save_state) {
			save_model_state();
		}
		display_helper();
	}
	if(sim_time <= run_time) {
		//Performs convection updates if the current simulation is using convection
		if(using_convection) {
			convection();
		}

		//Performs conduction calculations
		conduction();

		//Increments the simulation time and loop count
		sim_time += time_step;
		count++;
        
        if(use_tolerance == 1) {
            max_temp_diff = find_max_temp_diff();
            if(max_temp_diff < tolerance) {
                cout << "Maximum temperature change below the tolerance, stoping the simulation" << endl;
                //Saves the final result of the simulation
                if(save_state == 1 || save_result == 1) {
                    save_model_state();
                }
                save_surfer();
                cout << endl << "Simulation Complete" << endl;

                delete[] color_field;
                deallocate_memory();
                glutLeaveMainLoop();
            }
        }

        //Updates the moving source
        if(using_moving_source == 1) {
		    update_moving_sources();
        }
	}
	else {
		//Saves the final result of the simulation
		if(save_state == 1 || save_result == 1) {
			save_model_state();
		}
		save_surfer();
		cout << endl << "Simulation Complete" << endl;

		delete[] color_field;
        deallocate_memory();
		glutLeaveMainLoop();
	}
	glutPostRedisplay();
}


/*
 * This captures information when the mouse buttons are pressed.
 * Maintained for debugging.
 */
/*
void mouse_button(int button, int state, int x, int y) {
    int mods;

    if (state == GLUT_DOWN)
        buttonState |= 1<<button;
    else if (state == GLUT_UP)
        buttonState = 0;

    mods = glutGetModifiers();
    if (mods & GLUT_ACTIVE_SHIFT) 
    {
        buttonState = 2;
    } 
    else if (mods & GLUT_ACTIVE_CTRL) 
    {
        buttonState = 3;
    }

    ox = x; oy = y;

    glutPostRedisplay();
}
*/

/*
 * This captures mouse motion information.
 * Maintained for debugging
 */
/*
void mouse_move(int x, int y) {
    float dx = (float)(x - ox);
    float dy = (float)(y - oy);

    if (buttonState == 3) 
    {
        // left+middle = zoom
        camera_trans[2] += (dy / 100.0f) * 0.5f * fabs(camera_trans[2]);
    } 
    else if (buttonState & 2) 
    {
        // middle = translate
        camera_trans[0] += dx / 10.0f;
        camera_trans[1] -= dy / 10.0f;
    }
    else if (buttonState & 1) 
    {
        // left = rotate
        camera_rot[0] += dy / 5.0f;
        camera_rot[1] += dx / 5.0f;
    }

    ox = x; oy = y;
    glutPostRedisplay();
}
*/

/*
 *Standard keyboard character control
 */
void keyboard(unsigned char key, int x, int y) {
	switch(key) {
	case '-':
		if(current_slice > 0) {
			current_slice--;
			//camera_trans[1]-=0.5;//For camera mouse control, above
			//camera_trans[2]+=1.0;
			glMatrixMode(GL_PROJECTION);
			glLoadIdentity();
			gluPerspective(60, 1.77777f, 1.0, 20000.0);
			//gluLookAt(num_cols/2.0,num_rows*0.1,num_rows+current_slice,num_cols/2.0,-num_rows/3.0,current_slice,0.0,1.0,0.0);
            if(num_rows > num_cols) {
		        gluLookAt(num_cols/2.0,num_rows*0.1,num_rows+current_slice,num_cols/2.0,-num_rows/3.0,current_slice,0.0,1.0,0.0);
            }
            else {
                gluLookAt(num_cols/2.0,num_rows*0.1,num_cols+current_slice,num_cols/2.0,-num_rows/3.0,current_slice,0.0,1.0,0.0);
            }
		}
		break;
	case '+':
		if(current_slice < num_slices-1) {
			current_slice++;
			//camera_trans[1]+=0.5;//For camera mouse control, above
			//camera_trans[2]-=1.0;
			glMatrixMode(GL_PROJECTION);
			glLoadIdentity();
			gluPerspective(60, 1.77777f, 1.0, 20000.0);
			//gluLookAt(num_cols/2.0,num_rows*0.1,num_rows+current_slice,num_cols/2.0,-num_rows/3.0,current_slice,0.0,1.0,0.0);
            if(num_rows > num_cols) {
		        gluLookAt(num_cols/2.0,num_rows*0.1,num_rows+current_slice,num_cols/2.0,-num_rows/3.0,current_slice,0.0,1.0,0.0);
            }
            else {
                gluLookAt(num_cols/2.0,num_rows*0.1,num_cols+current_slice,num_cols/2.0,-num_rows/3.0,current_slice,0.0,1.0,0.0);
            }
		}
		break;
	case 'x':
		exit(0);
	default:
		break;
	}
	display_helper();
}

/*
 * Special keyboard control for arrows
 */
void keyboardSpecial(int key, int x, int y) {
	switch(key) {
	case GLUT_KEY_UP:
		if(current_slice > 0) {			
			current_slice--;
			//camera_trans[1]-=0.5;
			//camera_trans[2]+=1.0;
			glMatrixMode(GL_PROJECTION);
			glLoadIdentity();
			gluPerspective(60, 1.77777f, 1.0, 20000.0);
			//gluLookAt(num_cols/2.0,num_rows*0.1,num_rows+current_slice,num_cols/2.0,-num_rows/3.0,current_slice,0.0,1.0,0.0);
            if(num_rows > num_cols) {
		        gluLookAt(num_cols/2.0,num_rows*0.1,num_rows+current_slice,num_cols/2.0,-num_rows/3.0,current_slice,0.0,1.0,0.0);
            }
            else {
                gluLookAt(num_cols/2.0,num_rows*0.1,num_cols+current_slice,num_cols/2.0,-num_rows/3.0,current_slice,0.0,1.0,0.0);
            }
		}
		break;

	case GLUT_KEY_DOWN:
		if(current_slice < num_slices-1) {
			current_slice++;
			//camera_trans[1]+=0.5;
			//camera_trans[2]-=1.0;
			glMatrixMode(GL_PROJECTION);
			glLoadIdentity();
			gluPerspective(60, 1.77777f, 1.0, 20000.0);
			//gluLookAt(num_cols/2.0,num_rows*0.1,num_rows+current_slice,num_cols/2.0,-num_rows/3.0,current_slice,0.0,1.0,0.0);
            if(num_rows > num_cols) {
		        gluLookAt(num_cols/2.0,num_rows*0.1,num_rows+current_slice,num_cols/2.0,-num_rows/3.0,current_slice,0.0,1.0,0.0);
            }
            else {
                gluLookAt(num_cols/2.0,num_rows*0.1,num_cols+current_slice,num_cols/2.0,-num_rows/3.0,current_slice,0.0,1.0,0.0);
            }
		}
		break;
	default:
		break;
	}
	display_helper();
}
#endif


/**
 * Asynchronously Sends and receives edge data based on the processes rank
 */
void update_cells_mpi() {
	if(comm_sz > 1) {
		for(int j = 0; j < num_cols; j++) {
			if(my_rank == 0) {
				MPI_Isend(next_temp[chunk_size[my_rank]-1][j],num_slices,MPI_REAL_VALUE,my_rank+1,0,MPI_COMM_WORLD,&request1[j]);
				MPI_Irecv(next_temp[chunk_size[my_rank]][j],num_slices,MPI_REAL_VALUE,my_rank+1,0,MPI_COMM_WORLD,&request2[j]);
			}
			else if(my_rank == comm_sz-1) {
                MPI_Isend(next_temp[1][j],num_slices,MPI_REAL_VALUE,my_rank-1,0,MPI_COMM_WORLD,&request1[j]);
				MPI_Irecv(next_temp[0][j],num_slices,MPI_REAL_VALUE,my_rank-1,0,MPI_COMM_WORLD,&request2[j]);
			}
			else {
				MPI_Isend(next_temp[chunk_size[my_rank]][j],num_slices,MPI_REAL_VALUE,my_rank+1,0,MPI_COMM_WORLD,&request1[j]);
                MPI_Isend(next_temp[1][j],num_slices,MPI_REAL_VALUE,my_rank-1,0,MPI_COMM_WORLD,&request2[j]);
				MPI_Irecv(next_temp[0][j],num_slices,MPI_REAL_VALUE,my_rank-1,0,MPI_COMM_WORLD,&request3[j]);
				MPI_Irecv(next_temp[chunk_size[my_rank]+1][j],num_slices,MPI_REAL_VALUE,my_rank+1,0,MPI_COMM_WORLD,&request4[j]);
			}
		}
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
 * Swaps the temp arrays
 */
void swap_temp_array() {
    REAL ***tmp;
    
    tmp = temp;
    temp = next_temp;
    next_temp = tmp;
}

/*
 * Loads the input file into program memory and allocates
 * necessary memory to store the input variables
 */
void load_file() {
    ifstream source_file;       //Input file stream
    string temp_str;
    stringstream stringBuffer;
    REAL min_dim;
    REAL tmp_real;
	char tmp[500];
    
    //Ask for the input file names and displays an error message
    //if the file does not exist
	if(my_rank == 0) {
        cout << "Input File Name: " << source_filename << endl;
        source_file.open(source_filename.c_str(),ios::in);
        if(!source_file.is_open()) {
            cout << "File Not found!" << endl;
            MPI_Finalize();
            exit(1);
        }
		strcpy(tmp,source_filename.c_str());
		MPI_Bcast(tmp,500,MPI_CHAR,0,MPI_COMM_WORLD);
        MPI_Bcast(&save_result,1,MPI_INT,0,MPI_COMM_WORLD);

#ifdef TIME
#ifdef _WIN32
        GetSystemTime(&start);
#else
        gettimeofday(&start, NULL);
#endif
#endif
        /*
		do {
			cout << "Input File Name: ";
			//cin >> source_filename;
			source_file.open(source_filename.c_str(),ios::in);
			if(!source_file.is_open()) {
				cout << "File Not found!" << endl;
			}
		} while(!source_file.is_open());
		strcpy(tmp,source_filename.c_str());
		MPI_Bcast(tmp,500,MPI_CHAR,0,MPI_COMM_WORLD);
		
		//Asks the user if the state of the model should be saved every screen update
		/*cout << endl << "To save the state of the model every screen update enter 1, otherwise 0: ";
		while(!(cin >> save_state) || save_state < 0 || save_state > 1) {
			clear_cin();
			cout << "Incorrect input, to save the state of the model enter 1, else 0: ";
		}

		if(save_state == 0) {
			//Asks the user if the final result of the model should be saved
			cout << endl << "To save the final result of the model enter 1, otherwise 0: ";
			while(!(cin >> save_result) || save_result < 0 || save_result > 1) {
				clear_cin();
				cout << "Incorrect input, to save the final result of the model enter 1, else 0: ";
			}
		}
		
		//Ask for the state filename if the user specified that the file should be saved
		if(save_state == 1 || save_result == 1) {
			cout << "Output filename: ";
			cin >> output_filename;
		}
		*/
        if(save_result == 1) {
            cout << "Output filename: " << output_filename << endl;
        }
		//Asks for the DSAA surfer grid filenmae
		cout << "Surfer filename: " << output_su_filename << endl;;
		/*while(!(cin >> output_su_filename) || output_su_filename.length() < 5) {
			clear_cin();
			cout << "Please enter a filename at least 5 characters in length: ";
		}*/

		//Loads the input file
		cout << endl << endl << "Loading Input File";
	}
	else {
		MPI_Bcast(tmp,500,MPI_CHAR,0,MPI_COMM_WORLD);
        MPI_Bcast(&save_result,1,MPI_INT,0,MPI_COMM_WORLD);
		source_filename.assign(tmp);
		source_file.open(source_filename.c_str(),ios::in);
		if(!source_file.is_open()) {
			cout << "Process : " << my_rank << " File Not found!" << endl;
            exit(1);
		}
	}
    //Retrieves the simulation parameters from the input file
    source_file >> num_rows >> num_cols >> num_slices >> using_convection;
    source_file >> chf >> initial_time;
    getline(source_file,title);
    getline(source_file,title);
    
    //Allocates memory for MPI requests
    if(my_rank == 0 || my_rank == comm_sz -1) {
        request1 = new MPI_Request[num_cols];
        request2 = new MPI_Request[num_cols];
    }
    else {
        request1 = new MPI_Request[num_cols];
        request2 = new MPI_Request[num_cols];
        request3 = new MPI_Request[num_cols];
        request4 = new MPI_Request[num_cols];
    }
    
    
	if(my_rank == 0) {
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
	}
    chf *= 0.001;
    
    if(my_rank == 0) {
		//Calculates the number of characters for the surfer file index
		stringBuffer.str("");
        stringBuffer.clear();
        stringBuffer << num_slices;
        su_num_width = stringBuffer.str().length();
    }
	
    //Allocates memory for cell dimensions and distances
    dim_x = new REAL[num_cols];
    dim_y = new REAL[num_rows];
    dim_z = new REAL[num_slices];
    dist_x = new REAL[num_cols];
    dist_y = new REAL[num_rows];
    dist_z = new REAL[num_slices];
    
    //Reads in the Y (column) dimensions and finds the minimum column distance
    for(int i = 0; i < num_cols; i++) {
        source_file >> dim_x[i];
        if(i == 0) {
            min_col_dim = dim_x[0];
            dist_x[0] = dim_x[0]/2.0;
        }
        else {
            if(dim_x[i] < min_col_dim) {
                min_col_dim = dim_x[i];
            }
            dist_x[i] = dist_x[i-1] + dim_x[i-1]/2.0 + dim_x[i]/2.0;
        }
    }
    max_dist_x = dist_x[num_cols-1] + dim_x[num_cols-1]/2.0;

    //Reads in the X (row) dimensions and finds the minimum row distance
    for(int i = 0; i < num_rows; i++) {
        source_file >> dim_y[i];
        if(i == 0) {
            min_row_dim = dim_y[i];
            dist_y[0] = dim_y[0]/2.0;
        }
        else {
            if(dim_y[i] < min_row_dim) {
                min_row_dim = dim_y[i];
            }
            dist_y[i] = dist_y[i-1] + dim_y[i-1]/2.0 + dim_y[i]/2.0;
        }
    }
    max_dist_y = dist_y[num_rows-1] + dim_y[num_rows-1]/2.0;

    //Reads in the Z (slice depth) dimension and finds the minimum row distance
    for (int i = 0; i < num_slices; i++) {
        source_file >> dim_z[i];
        if (i == 0) {
            min_slice_dim = dim_z[i];
            dist_z[0] = dim_z[0]/2.0;
        }
        else {
            if (dim_z[i] < min_slice_dim) {
                min_slice_dim = dim_z[i];
            }
            dist_z[i] = dist_z[i-1] + dim_z[i-1]/2.0 + dim_z[i]/2.0;
        }
    }
    max_dist_z = dist_z[num_slices-1] + dim_z[num_slices-1]/2.0;

    //Reads in the conduction heat production values
    source_file >> num_hp;
    heat_production_values = new REAL[num_hp];
    for(int i = 0; i < num_hp; i++) {
        source_file >> heat_production_values[i];
        heat_production_values[i] /= 1E6;
    }
	if(my_rank == 0) {
		cout << "Read "<< num_hp << " heat production values" << endl;
    }
    //Reads in the thermal conduction difference values
    //Finds the minimum and maximum thermal conductivity differences and
    //performs some scaling of the conduction associated variables
    source_file >> num_tcd;
    thermal_conduct_diff = new REAL[num_tcd];
	if(my_rank == 0) {
		cout << "Converted " << num_tcd << " Thermal Conductivities to Diff. in m^2/y" << endl;
	}
    for(int i = 0; i < num_tcd; i++) {
        source_file >> thermal_conduct_diff[i];
        thermal_conduct_diff[i] *= 14.33;
        if(i == 0) {
            max_thermal_conduct_diff = thermal_conduct_diff[0];
            min_thermal_conduct_diff = thermal_conduct_diff[0];
        }
        else {
            if(thermal_conduct_diff[i] > max_thermal_conduct_diff) {
                max_thermal_conduct_diff = thermal_conduct_diff[i];
            }
            if(thermal_conduct_diff[i] < min_thermal_conduct_diff) {
                min_thermal_conduct_diff = thermal_conduct_diff[i];
            }
        }
		if(my_rank == 0) {
			cout << "  " << thermal_conduct_diff[i];
		}
    }

    //Calculates the maximum time step of the simulation
    min_dim = min_row_dim;
    if(min_col_dim < min_dim) {
        min_dim = min_col_dim;
    }
    if(num_slices > 1 && min_slice_dim < min_dim) {
        min_dim = min_slice_dim;
    }
    time_step = min_dim*min_dim/(5*max_thermal_conduct_diff);
    if(my_rank == 0) {
        if(temp_time_step < time_step) {
            time_step = temp_time_step;
        }
    }
    MPI_Bcast(&time_step,1,MPI_REAL_VALUE,0,MPI_COMM_WORLD);
    
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
		if(my_rank == 0) {
			cout << endl << "Read " << num_vel << " Velocities in m/yr" << endl;
        }
		
        //Finds the maximum convection velocity
        max_vel = vel[0];
		if(my_rank == 0) {
			cout << " " << vel[0];
		}
        for(int i = 1 ; i < num_vel; i++) {
            if(vel[i] > max_vel) {
                max_vel = vel[i];
            }
			if(my_rank == 0) {
				cout << " " << vel[i];
			}
        }
		if(my_rank == 0) {
			cout << endl;
		}
        
        /*
        T1 = max_thermal_conduct_diff;
        T2 = min_col_dim;
        T3 = min_row_dim;
        */
        //Finds the convection time increment
        tic = min_dim/max_vel;
        
        //Calculates the number of convection loops to perform per time step
        num_conv_loops = (int)(time_step/(10*tic));
        if(num_conv_loops > 5) {
            num_conv_loops = 5;
        }
        else if(num_conv_loops <= 0) {
            num_conv_loops = 1;
        }
    }
    
    //Calculates the chunk size and row displacements
	double count = 0.0; 
    double prev_count;
    chunk_size = new int[comm_sz];
	displacements = new int[comm_sz];
	int sum = 0;
    for(int i = 0; i < comm_sz; i++) {
		displacements[i] = (int)count;
        prev_count = count;
        count += (double)num_rows / (double)comm_sz;
        chunk_size[i] = (int)count - (int)prev_count;
		sum += chunk_size[i];
    }
	if(sum != num_rows) {
		chunk_size[comm_sz-1] += num_rows-sum;
		sum += num_rows-sum;
	}
    if(sum != num_rows) {
        cout << "Size mismatch" << endl;
        exit(1);
    }
	
    //Allocates memory for the conduction variables based on the previously read in simulation
    //parameters
    
    temp = new REAL**[chunk_size[my_rank]+2];
    next_temp = new REAL**[chunk_size[my_rank]+2];
    cond_codes = new int**[chunk_size[my_rank]+2];
    cond_hp_index = new int**[chunk_size[my_rank]+2];
    cond_tc_index = new int**[chunk_size[my_rank]+2];
    use_cond = new int**[chunk_size[my_rank]+2];
    for(int i = 0; i < chunk_size[my_rank]+2; i++) {
        temp[i] = new REAL*[num_cols];
        next_temp[i] = new REAL*[num_cols];
        cond_codes[i] = new int*[num_cols];
        cond_hp_index[i] = new int*[num_cols];
        cond_tc_index[i] = new int*[num_cols];
        use_cond[i] = new int*[num_cols];
        for (int j = 0; j < num_cols; j++) {
            temp[i][j] = new REAL[num_slices];
            next_temp[i][j] = new REAL[num_slices];
            cond_codes[i][j] = new int[num_slices];
            cond_hp_index[i][j] = new int[num_slices];
            cond_tc_index[i][j] = new int[num_slices];
            use_cond[i][j] = new int[num_slices];
        }
    }
    
	
    //Reads in the starting temperatures of the simulation from the input file
    for (int k = 0; k < num_slices; k++) {
        for(int i = 0; i < num_rows; i++) {
            for(int j = 0; j < num_cols; j++) {
				source_file >> tmp_real;
				if(my_rank == 0) {
					if(i >= 0 && i <= chunk_size[my_rank]) {
						temp[i][j][k] = tmp_real;
					}
				}
				else if(my_rank == comm_sz-1) {
					if(i >= displacements[my_rank]-1 && i <= displacements[my_rank] + chunk_size[my_rank]) {
						temp[i-displacements[my_rank]+1][j][k] = tmp_real;
					}
				}
				else {
					if(i >= displacements[my_rank]-1 && i <= displacements[my_rank] + chunk_size[my_rank]) {
						temp[i-displacements[my_rank]+1][j][k] = tmp_real;
					}
				}
            }
        }
    }
	if(my_rank == 0) {
		cout << "Read " << num_rows << " X " << num_cols << " X " << num_slices << " temps" << endl;
	}
	
	
    //Reads in the conduction codes for each cell of the simulation and parses
    //the array indexs from the codes
    //Unlike, the Fortran version of the program, the conduction direction codes
    //are ignored since the simulation accounts for them internally
    for (int k = 0; k < num_slices; k++) {
        for(int i = 0; i < num_rows; i++) {
            for(int j = 0; j < num_cols; j++) {
                source_file >> temp_str;
                if(temp_str.size() != 2*INDEX_WIDTH+1) {
                    if(my_rank == 0) {
                        cerr << "Incorrect Conduction code found at: " << i << "," << j << "," << k << endl;
                        cerr << "Code has character length of " << temp_str.size() << " instead of " << 2*INDEX_WIDTH+1 << endl;
                    }
                    exit(1);
                }
				if(my_rank == 0) {
					if(i >= 0 && i <= chunk_size[my_rank]) {
						cond_codes[i][j][k] = atoi(temp_str.c_str());
						cond_tc_index[i][j][k] = atoi(temp_str.substr(0*INDEX_WIDTH,INDEX_WIDTH).c_str())-1;
                        if(cond_tc_index[i][j][k] < 0 || cond_tc_index[i][j][k] >= num_tcd) {
                            cerr << "Incorrect Thermal Conductivity Diff. Index found at: " << i << "," << j << "," << k << endl;
                            exit(1);
                        }
						cond_hp_index[i][j][k] = atoi(temp_str.substr(1*INDEX_WIDTH,INDEX_WIDTH).c_str())-1;
                        if(cond_hp_index[i][j][k] < 0 || cond_hp_index[i][j][k] >= num_hp) {
                            cerr << "Incorrect Heat Production Value Index found at: " << i << "," << j << "," << k << endl;
                            exit(1);
                        }
                        use_cond[i][j][k] = atoi(temp_str.substr(2*INDEX_WIDTH,1).c_str());
                        if(use_cond[i][j][k] < 0 || use_cond[i][j][k] > 1) {
                            cerr << "Incorrect Conduction Enable flag found at: " << i << "," << j << "," << k << endl;
                            exit(1);
                        }
					}
				}
				else {
					if(i >= displacements[my_rank]-1 && i <= displacements[my_rank] + chunk_size[my_rank]) {
						cond_codes[i-displacements[my_rank]+1][j][k] = atoi(temp_str.c_str());
						cond_tc_index[i-displacements[my_rank]+1][j][k] = atoi(temp_str.substr(0*INDEX_WIDTH,INDEX_WIDTH).c_str())-1;
						if(cond_tc_index[i-displacements[my_rank]+1][j][k] < 0 || cond_tc_index[i-displacements[my_rank]+1][j][k] >= num_tcd) {
                            cerr << "Incorrect Thermal Conductivity Diff. Index found at: " << i << "," << j << "," << k << endl;
                            exit(1);
                        }
                        cond_hp_index[i-displacements[my_rank]+1][j][k] = atoi(temp_str.substr(1*INDEX_WIDTH,INDEX_WIDTH).c_str())-1;
						if(cond_hp_index[i-displacements[my_rank]+1][j][k] < 0 || cond_hp_index[i-displacements[my_rank]+1][j][k] >= num_hp) {
                            cerr << "Incorrect Heat Production Value Index found at: " << i << "," << j << "," << k << endl;
                            exit(1);
                        }
                        use_cond[i-displacements[my_rank]+1][j][k] = atoi(temp_str.substr(2*INDEX_WIDTH,1).c_str());
                        if(use_cond[i-displacements[my_rank]+1][j][k] < 0 || use_cond[i-displacements[my_rank]+1][j][k] > 1) {
                            cerr << "Incorrect Conduction Enable flag found at: " << i << "," << j << "," << k << endl;
                            exit(1);
                        }
					}
				}
            }
        }
    }
	if(my_rank == 0) {
		cout << "Read " << num_rows << " X " << num_cols << " X " << num_slices << " conduction codes" << endl;
	}
    
    //If convection is used for the user specified input file, memory is allocated for its
    //variables and they are read in from the input file
    if(using_convection) {     
        //Allocates memory for the convection variables based on the previously read in simulation
        //parameters
        conv_codes = new int**[chunk_size[my_rank]+2];
        conv_min_temp_index = new int**[chunk_size[my_rank]+2];
        conv_direction = new int**[chunk_size[my_rank]+2];
        conv_vel_index = new int**[chunk_size[my_rank]+2];
        conv_fluid_index = new int**[chunk_size[my_rank]+2];
        conv_rock_index = new int**[chunk_size[my_rank]+2];
        for(int i = 0; i < chunk_size[my_rank]+2; i++) {
            conv_codes[i] = new int*[num_cols];
            conv_min_temp_index[i] = new int*[num_cols];
            conv_direction[i] = new int*[num_cols];
            conv_vel_index[i] = new int*[num_cols];
            conv_fluid_index[i] = new int*[num_cols];
            conv_rock_index[i] = new int*[num_cols];
            for (int j = 0; j < num_cols; j++) {
                conv_codes[i][j] = new int[num_slices];
                conv_min_temp_index[i][j] = new int[num_slices];
                conv_direction[i][j] = new int[num_slices];
                conv_vel_index[i][j] = new int[num_slices];
                conv_fluid_index[i][j] = new int[num_slices];
                conv_rock_index[i][j] = new int[num_slices];
            }
        }
        
        //Reads in the convection codes for each cell of the simulation and parses the array
        //indexs from the ocdes
        for (int k = 0; k < num_slices; k++) {
            for(int i = 0; i < num_rows; i++) {
                for(int j = 0; j < num_cols; j++) {
                    source_file >> temp_str;
                    if(temp_str.size() != 5*INDEX_WIDTH) {
                        if(my_rank == 0) {
                            cerr << "Incorrect Convection code found at: " << i << "," << j << "," << k << endl;
                            cerr << "Code has character length of " << temp_str.size() << " instead of " << 5*INDEX_WIDTH << endl;
                        }
                        exit(1);
                    }
					if(my_rank == 0) {
						if(i >= 0 && i <= chunk_size[my_rank]) {
							conv_codes[i][j][k] = atoi(temp_str.c_str());
							conv_min_temp_index[i][j][k] = atoi(temp_str.substr(0*INDEX_WIDTH,INDEX_WIDTH).c_str())-1;
							if(conv_min_temp_index[i][j][k] < 0 || conv_min_temp_index[i][j][k] >= num_mtc) {
                                cerr << "Incorrect Minimum Convection Temp Index found at: " << i << "," << j << "," << k << endl;
                                exit(1);
                            }
                            conv_vel_index[i][j][k] = atoi(temp_str.substr(1*INDEX_WIDTH,INDEX_WIDTH).c_str())-1;
							if(conv_vel_index[i][j][k] < 0 || conv_vel_index[i][j][k] > num_vel) {
                                cerr << "Incorrect Convection Velocity Index found at: " << i << "," << j << "," << k << endl;
                                exit(1);
                            }
                            conv_fluid_index[i][j][k] = atoi(temp_str.substr(2*INDEX_WIDTH,INDEX_WIDTH).c_str())-1;
							if(conv_fluid_index[i][j][k] < 0 || conv_fluid_index[i][j][k] > num_hcf) {
                                cerr << "Incorrect Fluid Heat Capacity Index found at: " << i << "," << j << "," << k << endl;
                                exit(1);
                            }
                            conv_rock_index[i][j][k] = atoi(temp_str.substr(3*INDEX_WIDTH,INDEX_WIDTH).c_str())-1;
							if(conv_rock_index[i][j][k] < 0 || conv_rock_index[i][j][k] > num_vel) {
                                cerr << "Incorrect Rock Heat Capacity Index found at: " << i << "," << j << "," << k << endl;
                                exit(1);
                            }
                            conv_direction[i][j][k] = atoi(temp_str.substr(4*INDEX_WIDTH,2).c_str());
                            if(conv_direction[i][j][k] < 1 || conv_vel_index[i][j][k] > 27) {
                                cerr << "Incorrect Convection Direction found at: " << i << "," << j << "," << k << endl;
                                exit(1);
                            }
                        }
					}
					else {
						if(i >= displacements[my_rank]-1 && i <= displacements[my_rank] + chunk_size[my_rank]) {
							conv_codes[i-displacements[my_rank]+1][j][k] = atoi(temp_str.c_str());
							conv_min_temp_index[i-displacements[my_rank]+1][j][k] = atoi(temp_str.substr(0*INDEX_WIDTH,INDEX_WIDTH).c_str())-1;
							if(conv_min_temp_index[i-displacements[my_rank]+1][j][k] < 0 || conv_min_temp_index[i-displacements[my_rank]+1][j][k] >= num_mtc) {
                                cerr << "Incorrect Minimum Convection Temp Index found at: " << i << "," << j << "," << k << endl;
                                exit(1);
                            }
                            conv_vel_index[i-displacements[my_rank]+1][j][k] = atoi(temp_str.substr(1*INDEX_WIDTH,INDEX_WIDTH).c_str())-1;
							if(conv_vel_index[i-displacements[my_rank]+1][j][k] < 0 || conv_vel_index[i-displacements[my_rank]+1][j][k] > num_vel) {
                                cerr << "Incorrect Convection Velocity Index found at: " << i << "," << j << "," << k << endl;
                                exit(1);
                            }
                            conv_fluid_index[i-displacements[my_rank]+1][j][k] = atoi(temp_str.substr(2*INDEX_WIDTH,INDEX_WIDTH).c_str())-1;
							if(conv_fluid_index[i-displacements[my_rank]+1][j][k] < 0 || conv_fluid_index[i-displacements[my_rank]+1][j][k] > num_hcf) {
                                cerr << "Incorrect Fluid Heat Capacity Index found at: " << i << "," << j << "," << k << endl;
                                exit(1);
                            }
                            conv_rock_index[i-displacements[my_rank]+1][j][k] = atoi(temp_str.substr(3*INDEX_WIDTH,INDEX_WIDTH).c_str())-1;
							if(conv_rock_index[i-displacements[my_rank]+1][j][k] < 0 || conv_rock_index[i-displacements[my_rank]+1][j][k] > num_vel) {
                                cerr << "Incorrect Rock Heat Capacity Index found at: " << i << "," << j << "," << k << endl;
                                exit(1);
                            }
                            conv_direction[i-displacements[my_rank]+1][j][k] = atoi(temp_str.substr(4*INDEX_WIDTH,2).c_str());
                            if(conv_direction[i-displacements[my_rank]+1][j][k] < 1 || conv_direction[i-displacements[my_rank]+1][j][k] > 27) {
                                cerr << "Incorrect Convection Direction found at: " << i << "," << j << "," << k << endl;
                                exit(1);
                            }
                        }
					}
                }
            }
        }
		if(my_rank == 0) {
			cout << "Read " << num_rows << " X " << num_cols << " X " << num_slices << " convection codes" << endl;
		}
    }
    
    //Closes the input file
    source_file.close();

	if(my_rank == 0) {
		cout << endl << "Done Loading Input File" << endl;
	}
}

/*
 * Saves the current state of the simulation, using the same format
 * as the input file
 */
void save_model_state() {
    ofstream output_file;    //Output file stream
    REAL tmp_array[chunk_size[my_rank]*num_cols];
	REAL tmp_array1[num_rows*num_cols];
    int code_array[chunk_size[my_rank]*num_cols];
	int code_array1[num_rows*num_cols];
    //Opens the output file for writing
    if(my_rank == 0) {
        output_file.open(output_filename.c_str(),ios::out);
        if(!output_file.is_open()) {
            cerr << "Failed to write state to file" << endl;
            exit(1);
        }
        else {
            //Prints the simulation parameters to the output file
            output_file << setw(20) << num_rows << " " << setw(20) << num_cols << " " << setw(20) << num_slices << setw(20) << using_convection << endl;
            output_file << setw(20) << fixed << setprecision(OUT_PRECISION) << chf*1000.0 << " " << setw(20) << initial_time + sim_time << endl;
            output_file << title << endl;
            
            
            output_file << setfill(' ');
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
            
            // Prints the slice (Z) dimensions of the simulation to the oputput file

            for (int i = 0; i < num_slices; i++) {
                output_file << " " << dim_z[i];
            }
            output_file << endl;

            //Prints the heat production values of the simulation to the output file
            output_file << " " << num_hp;
            for(int i = 0; i < num_hp; i++) {
                output_file << " " << scientific << heat_production_values[i]*1E6;
            }
            output_file << endl;
            
            //Prints the thermal conductivity difference values to the output file
            output_file << " " << num_tcd;
            for(int i = 0; i < num_tcd; i++) {
                output_file << " " << thermal_conduct_diff[i]/14.33;
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
            for (int k = 0; k < num_slices; k++) {
                for(int i = 0; i < chunk_size[my_rank]; i++) {
                    for(int j = 0; j < num_cols; j++) {
                        tmp_array1[i*num_cols + j] = temp[i][j][k];
                    }
                }
                for(int i = 1; i < comm_sz; i++) {
                    MPI_Recv(&tmp_array1[displacements[i]*num_cols],num_cols*chunk_size[i],MPI_REAL_VALUE,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                }
                
                for(int i = 0; i < num_rows; i++) {
                    for(int j = 0; j < num_cols; j++) {
                        output_file << " " << setw(OUT_PRECISION+5) << tmp_array1[i*num_cols + j];
                    }
                    output_file << endl;
                }
                output_file << endl;
            }
            
            //Prints the conduction codes of the simulation to the output file
            output_file << setfill('0');
            for (int k = 0; k < num_slices; k++) {
                for(int i = 0; i < chunk_size[my_rank]; i++) {
                    for(int j = 0; j < num_cols; j++) {
                        code_array1[i*num_cols + j] = cond_codes[i][j][k];
                    }
                }
                for(int i = 1; i < comm_sz; i++) {
                    MPI_Recv(&code_array1[displacements[i]*num_cols],num_cols*chunk_size[i],MPI_INT,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                }
                
                for(int i = 0; i < num_rows; i++) {
                    for(int j = 0; j < num_cols; j++) {
                        output_file << " " << setw(2*INDEX_WIDTH+1) << code_array1[i*num_cols + j];
                    }
                    output_file << endl;
                }
                output_file << endl;
            }
            
            //Prints the convection codes to the output file if convection is being used
            if(using_convection) {
                for (int k = 0; k < num_slices; k++) {
                    for(int i = 0; i < chunk_size[my_rank]; i++) {
                        for(int j = 0; j < num_cols; j++) {
                            code_array1[i*num_cols + j] = conv_codes[i][j][k];
                        }
                    }
                    for(int i = 1; i < comm_sz; i++) {
                        MPI_Recv(&code_array1[displacements[i]*num_cols],num_cols*chunk_size[i],MPI_INT,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                    }
                    for(int i = 0; i < num_rows; i++) {
                        for(int j = 0; j < num_cols; j++) {
                            output_file << " " << setw(4*INDEX_WIDTH+2) << code_array1[i*num_cols + j];
                        }
                        output_file << endl;
                    }
                    output_file << endl;
                }
            }
            
            //Closes the output file
            output_file.close();
        }
    }
    else {
        for (int k = 0; k < num_slices; k++) {
            for(int i = 1; i <= chunk_size[my_rank]; i++) {
                for(int j = 0; j < num_cols; j++) {
                    tmp_array[(i-1)*num_cols + j] = temp[i][j][k];
                }
            }
            MPI_Send(tmp_array,num_cols*chunk_size[my_rank],MPI_REAL_VALUE,0,0,MPI_COMM_WORLD);
        }
        for (int k = 0; k < num_slices; k++) {
            for(int i = 1; i <= chunk_size[my_rank]; i++) {
                for(int j = 0; j < num_cols; j++) {
                    code_array[(i-1)*num_cols + j] = cond_codes[i][j][k];
                }
            }
            MPI_Send(code_array,num_cols*chunk_size[my_rank],MPI_INT,0,0,MPI_COMM_WORLD);
        }
        for (int k = 0; k < num_slices; k++) {
            for(int i = 1; i <= chunk_size[my_rank]; i++) {
                for(int j = 0; j < num_cols; j++) {
                    code_array[(i-1)*num_cols + j] = conv_codes[i][j][k];
                }
            }
            MPI_Send(code_array,num_cols*chunk_size[my_rank],MPI_INT,0,0,MPI_COMM_WORLD);
        }
    }
}

/*
 * Saves the current temperatures of the simulation to a DSAA surfer grid file
 */
void save_surfer() {
    ofstream output_file;            //Output file stream
    ostringstream oss;
    string filename, extension;
	REAL tmp_array[chunk_size[my_rank]*num_cols];
	REAL tmp_array1[num_rows*num_cols];
	if(my_rank == 0) {
		filename = output_su_filename.substr(0,output_su_filename.length()-4);
		extension = output_su_filename.substr(output_su_filename.length()-4,4);
    }
    for(int k = 0; k < num_slices; k++) {
		if(my_rank == 0) {
			oss.str("");
			oss.clear();
			oss << filename << setfill('0') << setw(su_num_width) << k << extension;
			for(int i = 0; i < chunk_size[my_rank]; i++) {
				for(int j = 0; j < num_cols; j++) {
					tmp_array1[i*num_cols + j] = temp[i][j][k];
				}
			}
			for(int i = 1; i < comm_sz; i++) {
				MPI_Recv(&tmp_array1[displacements[i]*num_cols],num_cols*chunk_size[i],MPI_REAL_VALUE,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			}
			//Opens the output file for writting
			output_file.open(oss.str().c_str(),ios::out);
			if(!output_file.is_open()) {
				cerr << "Failed to write surfer file" << endl;
				exit(1);
			}
			else {
				REAL min_temp, max_temp, temp_range;    //Minimum and maximum temperatures.
				REAL xmax,ymin;                         //Maximum x and minimum y distances
				
				//Finds the minimum and maximum temps in the temperature array
				min_temp = max_temp = tmp_array1[0];
				for(int i = 0; i < num_rows*num_cols; i++) {
					if(tmp_array1[i] > max_temp) {
						max_temp = tmp_array1[i];
					}
					if(tmp_array1[i] < min_temp) {
						min_temp = tmp_array1[i];
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
						output_file << " " << setw(OUT_PRECISION+5) << tmp_array1[i*num_cols + j];
					}
					output_file << endl;
				}
				
				//Closes the output file
				output_file.close();
			}
		}
		else {
			for(int i = 1; i <= chunk_size[my_rank]; i++) {
				for(int j = 0; j < num_cols; j++) {
					tmp_array[(i-1)*num_cols + j] = temp[i][j][k];
				}
			}
			MPI_Send(tmp_array,num_cols*chunk_size[my_rank],MPI_REAL_VALUE,0,0,MPI_COMM_WORLD);
		}
    }
}

/*
 * Calculates and returns the heat flow per year between two cells in the X direction
 * based on the provided indexes
 */
REAL cond_add_x(int row1, int col1, int slice1, int row2, int col2, int slice2) {
    REAL temp_diff;    //Temperature difference between the two cells
    REAL ad;           //
    
    temp_diff = temp[row1][col2][slice1] - temp[row1][col1][slice1];
    ad = dim_x[col2]/thermal_conduct_diff[cond_tc_index[row1][col2][slice1]] + dim_x[col1]/thermal_conduct_diff[cond_tc_index[row1][col1][slice1]];
    
    return 2*temp_diff/(ad*dim_x[col1]);
}

/*
 * Calculates and returns the heat flow per year between two cells in the Y direction
 * based on the provided indexes
 */
REAL cond_add_y(int row1, int col1, int slice1, int row2, int col2, int slice2) {
    REAL temp_diff;    //Temperature difference between the two cells
    REAL ad;           //
    
    temp_diff = temp[row2][col1][slice1] - temp[row1][col1][slice1];
    ad = dim_y[row2]/thermal_conduct_diff[cond_tc_index[row2][col1][slice1]] + dim_y[row1]/thermal_conduct_diff[cond_tc_index[row1][col1][slice1]];
    
    return 2*temp_diff/(ad*dim_y[row1]);
}

/*
 * Calculates and returns the heat flow per year between two cells in the Z direction
 * based on the provided indexes
 */
REAL cond_add_z(int row1, int col1, int slice1, int row2, int col2, int slice2) {
    if(num_slices == 1) {
        return 0.0;
    }
    REAL temp_diff;
    REAL ad;
    temp_diff = temp[row1][col1][slice2] - temp[row1][col1][slice1];
    ad = dim_z[slice2]/thermal_conduct_diff[cond_tc_index[row1][col1][slice2]] + dim_z[slice1]/thermal_conduct_diff[cond_tc_index[row1][col1][slice1]];
    return 2*temp_diff/(ad*dim_z[slice1]);
}

/*
 * Calculates the in-plane heat flow due to conduction in a given slice k.
 * If slices == 1
 *   2d simulation, return 0 for 3rd dimension heat transfer
 * else
 *   Calculate and return heat flow per year between two cells in the Z direction
 */
REAL in_plane_cond(int i, int j, int k) {
    REAL heat_flow_x;
    REAL heat_flow_y;

    /* k is fixed */
    if(my_rank == 0 && i == 0 && j == 0) { //Top left corner of slice
        heat_flow_x = cond_add_x(i,j,k,i,j+1,k);
        heat_flow_y = cond_add_y(i,j,k,i+1,j,k);
        next_temp[i][j][k] = 0.0;
    }
    else if(my_rank == 0 && i == 0 && j == num_cols-1) { //Top right corner of slice
        heat_flow_x = cond_add_x(i,j,k,i,j-1,k);
        heat_flow_y = cond_add_y(i,j,k,i+1,j,k);
        next_temp[i][j][k] = 0.0;
    }
    else if(my_rank == 0 && i == 0) { //Top of slice
        heat_flow_x = cond_add_x(i,j,k,i,j+1,k) + cond_add_x(i,j,k,i,j-1,k);
        heat_flow_y = cond_add_y(i,j,k,i+1,j,k);
        next_temp[i][j][k] = 0.0;
    }
    else if(((comm_sz == 1 && i == num_rows-1) || (my_rank == comm_sz-1 && i == chunk_size[my_rank])) && j == 0) { //Bottom left corner of slice
        heat_flow_x = cond_add_x(i,j,k,i,j+1,k);
        heat_flow_y = cond_add_y(i,j,k,i-1,j,k);
        next_temp[i][j][k] = DHF/dim_y[i];    //Constant heat flow at the bottom of the model
    }
    else if(((comm_sz == 1 && i == num_rows-1) || (my_rank == comm_sz-1 && i == chunk_size[my_rank])) && j == num_cols-1) { //Bottom right corner of slice
        heat_flow_x = cond_add_x(i,j,k,i,j-1,k);
        heat_flow_y = cond_add_y(i,j,k,i-1,j,k);
        next_temp[i][j][k] = DHF/dim_y[i];    //Constant heat flow at the bottom of the model
    }
    else if(((comm_sz == 1 && i == num_rows-1) || (my_rank == comm_sz-1 && i == chunk_size[my_rank]))) { //Bottom
        heat_flow_x = cond_add_x(i,j,k,i,j+1,k) + cond_add_x(i,j,k,i,j-1,k);
        heat_flow_y = cond_add_y(i,j,k,i-1,j,k);
        next_temp[i][j][k] = DHF/dim_y[i];    //Constant heat flow at the bottom of the model
    }
    else if(j == 0) { //Left side of slice
        heat_flow_x = cond_add_x(i,j,k,i,j+1,k);
        heat_flow_y = cond_add_y(i,j,k,i-1,j,k) + cond_add_y(i,j,k,i+1,j,k);
        next_temp[i][j][k] = 0.0;
    }
    else if(j == num_cols-1) { //Right side of slice
        heat_flow_x = cond_add_x(i,j,k,i,j-1,k);
        heat_flow_y = cond_add_y(i,j,k,i-1,j,k) + cond_add_y(i,j,k,i+1,j,k);
        next_temp[i][j][k] = 0.0;
    }
    else { //Middle of slice
        heat_flow_x = cond_add_x(i,j,k,i,j-1,k) + cond_add_x(i,j,k,i,j+1,k);
        heat_flow_y = cond_add_y(i,j,k,i-1,j,k) + cond_add_y(i,j,k,i+1,j,k);
        next_temp[i][j][k] = 0.0;
    }
    return (heat_flow_x + heat_flow_y);
}

void conduction_helper(int i, int j, int k) {
    if(use_cond[i][j][k] == 1) {
        REAL heatflow_in_plane;    //Heat flow occuring inside of plane
        REAL heatflow_cross_plane; //Heat flow into and out of plane/slice
        //Calculates heat flow in the X, Y, and Z direction into the current
        //cell based on its location within the model
        if (k == 0) { // First slice
            heatflow_in_plane = in_plane_cond(i,j,k);         // heat transfer inside of plane
            heatflow_cross_plane = cond_add_z(i,j,k,i,j,k+1); // slice-to-slice heat transfer.  first slice, so only from next slice transfers heat.
        }
        else if (k == num_slices - 1) {   // Last slice
            heatflow_in_plane = in_plane_cond(i,j,k);         // heat transfer inside of plane
            heatflow_cross_plane = cond_add_z(i,j,k,i,j,k-1); // slice-to-slice heat transfer.  last slice, so only previous slice transfers heat.
        }
        else {  // Middle
            heatflow_in_plane = in_plane_cond(i,j,k);                                     // you get the idea
            heatflow_cross_plane = cond_add_z(i,j,k,i,j,k+1) + cond_add_z(i,j,k,i,j,k-1); // slice-to-slice heat transfer. Middle, so both next and previous.
        }
        //Heat flow from the adjacent cells
        next_temp[i][j][k] += temp[i][j][k] + time_step*(heatflow_in_plane + heatflow_cross_plane);
        //Heat flow due to radioactive heat production
        next_temp[i][j][k] += heat_production_values[cond_hp_index[i][j][k]]*time_step/DTC;
    }
    else {
        next_temp[i][j][k] = temp[i][j][k];
    }
}

/*
 * Updates the temperature array using 3D conduction with finite
 * difference heat flow.
 */
void conduction() {
	for (int k = 0; k < num_slices; k++) {
        for(int j = 0; j < num_cols; j++) {
            conduction_helper(st_index,j,k);
            conduction_helper(en_index-1,j,k);
        }
    }
    update_cells_mpi();
    for (int k = 0; k < num_slices; k++) {
        for(int i = st_index+1; i < en_index-1; i++) {
            for(int j = 0; j < num_cols; j++) {
                conduction_helper(i,j,k);
            }
        }
    }
    for(int i = 0; i < num_cols; i++) {
        if(my_rank == 0 || my_rank == comm_sz - 1) {
            MPI_Wait(&request1[i],MPI_STATUS_IGNORE);
            MPI_Wait(&request2[i],MPI_STATUS_IGNORE);
        }
        else {
            MPI_Wait(&request1[i],MPI_STATUS_IGNORE);
            MPI_Wait(&request2[i],MPI_STATUS_IGNORE);
            MPI_Wait(&request3[i],MPI_STATUS_IGNORE);
            MPI_Wait(&request4[i],MPI_STATUS_IGNORE);
        }
    }
    swap_temp_array();    //Swaps the current and next temperature arrays
}

/*
 * Performs convection between two specified cells
 */
void perform_convection(int row1, int col1, int slice1, int row2, int col2, int slice2) {
    REAL avg_x_dim;    //distance between two temperature cells in the x direction
    REAL avg_y_dim;    //distance between two temperature cells in the y direction
    REAL avg_z_dim;    //distance between two temperature cells in the z direction
    REAL amt;          //
    REAL dist;         //Distance between the two cells
    REAL ratio;        //Ratio of amt to distance
    
    //Checks if the specified cell is within the bounds of the simulation and if it has a high enough
    //temperature to perform convection
	
    if((row2 >= 0) && ((my_rank != comm_sz-1) || ((comm_sz == 1) && (row2 < chunk_size[my_rank])) || ((comm_sz > 1) && (row2 <= chunk_size[my_rank]))) && (col2 >= 0) && (col2 < num_cols) && (slice2 >= 0) && (slice2 < num_slices) && (temp[row2][col2][slice2] - min_temp_conv[conv_min_temp_index[row1][col1][slice1]] >= 0)) {
        avg_x_dim = dist_x[col1] - dist_x[col2];
        avg_y_dim = dist_y[row1] - dist_y[row2];
        avg_z_dim = dist_z[slice1] - dist_z[slice2];

        amt = (vel[conv_vel_index[row1][col1][slice1]]*heat_capac_fluid[conv_fluid_index[row1][col1][slice1]]/heat_capac_rock[conv_rock_index[row1][col1][slice1]])*time_inc;
        dist = sqrt(avg_x_dim*avg_x_dim + avg_y_dim*avg_y_dim + avg_z_dim*avg_z_dim);
        ratio = amt/dist;
        if(ratio > 1) {
            ratio = 0.999999;
        }
        next_temp[row1][col1][slice1] = temp[row1][col1][slice1] + ratio *(temp[row2][col2][slice2]-temp[row1][col1][slice1]);
    }
    else {
        next_temp[row1][col1][slice1] = temp[row1][col1][slice1];
    }
}

void convection_helper(int i, int j, int k) {
    //Checks if convection can occur for the specified cell
    if((conv_codes[i][j][k] <= 0) || (i == 0) || (conv_direction[i][j][k] == 5) || (conv_direction[i][j][k] < 1) || (conv_direction[i][j][k] > 27)) {
        next_temp[i][j][k] = temp[i][j][k];
    }
    else {
        //Performs convection based on the convection direction code
        switch(conv_direction[i][j][k]) {
             /**
              * IN-PLANE convection -- 1 through 9. These codes are for convection taking place in the current, "k-th" plane
              *     1   2   3
              *     4   5   6
              *     7   8   9
              */
            case 1:
                perform_convection(i,j,k,i-1,j-1,k);
                break;
            case 2:                                                           
                perform_convection(i,j,k,i-1,j,k);                      
                break;                                                        
            case 3:                                                  
                perform_convection(i,j,k,i-1,j+1,k); 
                break;
            case 4:
                perform_convection(i,j,k,i,j-1,k);
                break;
            case 6:
                perform_convection(i,j,k,i,j+1,k);
                break;
            case 7:
                perform_convection(i,j,k,i+1,j-1,k);
                break;
            case 8:
                perform_convection(i,j,k,i+1,j,k);
                break;
            case 9:
                perform_convection(i,j,k,i+1,j+1,k);
                break;
            /**
             * CROSS-PLANE convection (previous "k-1th" plane) -- 10 through 18
             *      10  11  12
             *      13  14  15
             *      16  17  18      
             */
            case 10:                                   
                perform_convection(i,j,k,i-1,j-1,k-1);
                break;
            case 11:
                perform_convection(i,j,k,i-1,j,k-1);                      
                break;
            case 12:
                perform_convection(i,j,k,i-1,j+1,k-1); 
                break;
            case 13:
                perform_convection(i,j,k,i,j-1,k-1);
                break;
            case 14:
                perform_convection(i,j,k,i,j,k-1);
                break;
            case 15:
                perform_convection(i,j,k,i,j+1,k-1);
                break;
            case 16:
                perform_convection(i,j,k,i+1,j-1,k-1);
                break;
            case 17:
                perform_convection(i,j,k,i+1,j,k-1);
                break;
            case 18:
                perform_convection(i,j,k,i+1,j+1,k-1);
                break;
            /**
             * CROSS-PLANE convection ("k+1th" plane) -- 19 through 27 
             *      19  20  21
             *      22  23  24
             *      25  26  27
             */
            case 19:                                   
                perform_convection(i,j,k,i-1,j-1,k-1);
                break;
            case 20:
                perform_convection(i,j,k,i-1,j,k-1);                      
                break;
            case 21:
                perform_convection(i,j,k,i-1,j+1,k-1); 
                break;
            case 22:
                perform_convection(i,j,k,i,j-1,k-1);
                break;
            case 23:
                perform_convection(i,j,k,i,j,k-1);
                break;
            case 24:
                perform_convection(i,j,k,i,j+1,k-1);
                break;
            case 25:
                perform_convection(i,j,k,i+1,j-1,k-1);
                break;
            case 26:
                perform_convection(i,j,k,i+1,j,k-1);
                break;
            case 27:
                perform_convection(i,j,k,i+1,j+1,k-1);
                break;
        }
    }
}

/*
 * Updates the temperature array using convection
 */
void convection() { 
    //Performs the calculated number of convection updates per time step
    for(int n = 0; n < num_conv_loops; n++) {
        for(int k = 0; k < num_slices; k++) {
            for(int j = 0; j < num_cols; j++) {
                convection_helper(st_index,j,k);
                convection_helper(en_index-1,j,k);
            }
        }
        update_cells_mpi();
        for (int k = 0; k < num_slices; k++) {
            for(int i = st_index+1; i < en_index-1; i++) {
                for(int j = 0; j < num_cols; j++) {
                    convection_helper(i,j,k);
                }
            }
        }
        for(int i = 0; i < num_cols; i++) {
            if(my_rank == 0 || my_rank == comm_sz - 1) {
                MPI_Wait(&request1[i],MPI_STATUS_IGNORE);
                MPI_Wait(&request2[i],MPI_STATUS_IGNORE);
            }
            else {
                MPI_Wait(&request1[i],MPI_STATUS_IGNORE);
                MPI_Wait(&request2[i],MPI_STATUS_IGNORE);
                MPI_Wait(&request3[i],MPI_STATUS_IGNORE);
                MPI_Wait(&request4[i],MPI_STATUS_IGNORE);
            }
        }
        swap_temp_array();    //Swaps the current and next temperature arrays
    }
}

/*
 * Finds and returns the maximum temperature difference between
 * the current and next temperature arrays.
 */
REAL find_max_temp_diff() {
    REAL max_diff = fabs(next_temp[0][0][0] - temp[0][0][0]);
    REAL diff = 0.0;
    for(int k = 0; k < num_slices; k++) {
        for(int i = 0; i < num_rows; i++) {
            for(int j = 0; j < num_cols; j++) {
                diff = fabs(next_temp[i][j][k] - temp[i][j][k]);
                if(diff > max_diff) {
                    max_diff = diff;
                }
            }
        }
    }
    return max_diff;
}

/*
 * Finds the index of a given x, y, and z value in meters and
 * stores them in the index array
 */
void find_loc_index(REAL x_loc, REAL y_loc, REAL z_loc, int *index){
	if(x_loc < 0) {
		index[0] = -1;
	}
	else if(x_loc > max_dist_x) {
		index[0] = num_cols;
	}
	else {
		for(index[0] = 0; index[0] < num_cols; index[0]++) {
			if(x_loc <= dist_x[index[0]]+dim_x[index[0]]/2.0) {
				break;
			}
		}
	}

	if(y_loc < 0) {
		index[1] = -1;
	}
	else if(y_loc > max_dist_y) {
		index[1] = num_rows;
	}
	else {
		for(index[1] = 0; index[1] < num_rows; index[1]++) {
			if(y_loc <= dist_y[index[1]]+dim_y[index[1]]/2.0) {
				break;
			}
		}
	}

	if(z_loc < 0) {
		index[2] = -1;
	}
	else if(z_loc > max_dist_z) {
		index[2] = num_slices;
	}
	else {
		for(index[2] = 0; index[2] < num_slices; index[2]++) {
			if(z_loc <= dist_x[index[2]]+dim_x[index[2]]/2.0) {
				break;
			}
		}
	}
	
}


/*
 * Finds the indexes of two corners of the moving source
 * if either falls within the model. The valid parts of the
 * moving source are updated with the moving sources temperature
 */
void update_mvsrc(int index) {
    if(mvsrc_valid[index] == 1) {
        int loc_index[3], loc_offset_index[3];
        find_loc_index(mvsrc_x[index],mvsrc_y[index],mvsrc_z[index],loc_index);
        find_loc_index(mvsrc_x[index]+mvsrc_offset_x[index],mvsrc_y[index]+mvsrc_offset_y[index],mvsrc_z[index]+mvsrc_offset_z[index],loc_offset_index);

        if((loc_index[0] >= 0 && loc_index[0] < num_cols && loc_index[1] >= 0 && loc_index[1] < num_rows && loc_index[2] >= 0 && loc_index[2] < num_slices) || (loc_offset_index[0] >= 0 && loc_offset_index[0] < num_cols && loc_offset_index[1] >= 0 && loc_offset_index[1] < num_rows && loc_offset_index[2] >= 0 && loc_offset_index[2] < num_slices)) {
            for(int k = loc_index[2]; k <=  loc_offset_index[2]; k++) {
                for(int i = loc_index[1]; i <=  loc_offset_index[1]; i++) {
                    for(int j = loc_index[0]; j <=  loc_offset_index[0]; j++) {
                        if(i >= 0 && i < num_rows && j >= 0 && j < num_cols && k >= 0 && k < num_slices) {
                            temp[i][j][k] = mvsrc_temp[index];
                        }
                    }
                }
            }
        }
        else {
            mvsrc_valid[index] = 0;
        }
    }
}

/*
 * Updates the moving sources velocity and position vectors
 * then updates the temperatures in the current temp array
 */
void update_moving_sources() {
    for(int i = 0; i < num_mvsrc; i++) {
        mvsrc_vel_x[i] += mvsrc_accel_x[i]*time_step;
        mvsrc_vel_y[i] += mvsrc_accel_y[i]*time_step;
        mvsrc_vel_z[i] += mvsrc_accel_z[i]*time_step;

        mvsrc_x[i] += mvsrc_vel_x[i]*time_step;
        mvsrc_y[i] += mvsrc_vel_y[i]*time_step;
        mvsrc_z[i] += mvsrc_vel_z[i]*time_step;

        update_mvsrc(i);
    }
}


void usage() {
    cerr << "Arguments: " << endl;
    cerr << "   --input_filename       :       The filename of the input file" << endl;
    cerr << "   --surfer_filename      :       The filename of the output surfer file" << endl;
    cerr << "   --output_filename      :       The filename of the output state file" << endl;
    cerr << "   --save_result          :       Flag indicating if the final state of the simulation should be saved" << endl;
    cerr << "   --time_step            :       The new timestep of the model" << endl;
    cerr << "   --run_time             :       The run time of the simulation" << endl;
    cerr << "   --num_loops            :       The number of loops between screen updates" << endl;
}

/*
 * Performs a finite heat flow simulation using
 * conduction and convection.
 */
int main(int argc, char **argv) {

#ifdef DISPLAY
    cout << "\t\t Finite Difference Heat Flow Simulation" << endl;
	//Asks the user if they wish to visualize results
    cout << endl << "Press 1 to run visualization, otherwise 0: ";
    while(!(cin >> display_mode) || display_mode < 0 || display_mode > 1) {
		clear_cin();
        cout << "Incorrect input, to save the state of the model enter 1, else 0: ";
    }
#else

#endif
    int input_val;    //Temporary int value
    REAL temp_val = -1;    //Temporary REAL value
    
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--input_filename") == 0) {
            source_filename = argv[++i];
        } else if (strcmp(argv[i], "--surfer_filename") == 0) {
            output_su_filename = argv[++i];
        } else if (strcmp(argv[i], "--output_filename") == 0) {
            output_filename = argv[++i];
        } else if (strcmp(argv[i], "--save_result") == 0) {
            save_result = 1;
        } else if (strcmp(argv[i], "--time_step") == 0) {
            temp_time_step = atof(argv[++i]);
        } else if (strcmp(argv[i], "--run_time") == 0) {
            run_time = atof(argv[++i]);
        } else if (strcmp(argv[i], "--num_loops") == 0) {
            num_loops = atoi(argv[++i]);
        } else {
            cerr << "Unknown argument " << argv[i] << endl;
            usage();
            exit(0);
        }
    }
    
    if(source_filename.size() < 5) {
        cerr << "incorrect input filename specified" << endl;
        exit(0);
    }
    
    if(output_su_filename.size() == 0) {
        cerr << "incorrect surfer filename specified" << endl;
        exit(0);
    }
    
    if(save_result == 1 && output_filename.size() == 0) {
        cerr << "Output filename not specified, exiting program" << endl;
        exit(0);
    }
    
    if(temp_time_step <= 0) {
        cerr << "Time step not specified" << endl;
        exit(0);
    }
    
    if(run_time <= 0) {
        cerr << "Run time not specified" << endl;
        exit(0);
    }
    
    if(num_loops <= 0) {
        cerr << "Number of loops not specified" << endl;
        exit(0);
    }
    
    MPI_Init(NULL, NULL);						//Initializes MPI
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);	//Retrieves the number of processes running
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);	//Retrieves the rank of the current process	
    
	if(my_rank == 0) {
		cout << "\t\t Finite Difference Heat Flow Simulation" << endl;
    }
	
    //Loads the input file for the simulation
    load_file();
#ifdef TIME
    if(my_rank == 0) {
#ifdef _WIN32
        GetSystemTime(&end);
        cout << "Load and allocate time : " << (end.wHour-start.wHour)*3600.0 + (end.wMinute-start.wMinute)*60.0 + (end.wSecond-start.wSecond) + (end.wMilliseconds-start.wMilliseconds)/1000.0 << endl;
#else
        gettimeofday(&end, NULL);
        cout << "Load and allocate time : " << (end.tv_sec  - start.tv_sec) + (end.tv_usec - start.tv_usec)/1e6 << endl;
#endif
    }
#endif
    /*
     * Allows the user to change multiple rectangular blocks of temperatures
     * within the model
     */
	 /*
    cout << endl << endl << "To Change the Temp. on a Block, Enter 1, Else 0: ";
    while(!(cin >> input_val) || input_val < 0 || input_val > 1) {
		clear_cin();
        cout << "Incorrect Input, Enter 1 to Change, Else 0: ";
    }
    
    //Warning, the row column pairs need to be space seperated not comma seperated  
    if(input_val == 1) {
        int num_block, row1, row2, col1, col2, slice1, slice2;
        REAL new_temp;
        cout << "Enter the Number of Blocks to Change: ";
        while(!(cin >> num_block) || num_block < 0) {
            clear_cin();
			cout << "Enter a number greater than or equal to 0: ";
		}
        for(int i = 0; i < num_block; i++) {
            cout << endl << "Block " << i << endl;
            cout << "Enter the Coordinates of the Upper Left Corner <row> <column> <slice>: ";
			while(!(cin >> row1 >> col1 >> slice1) || row1 < 0 || col1 < 0 || slice1 < 0) {
                clear_cin();
				cout << "Incorrect input, enter three positive numbers with spaces: ";
			}
            cout << "Enter the Coordinates of the Lower Right Corner <row> <column> <slice>: ";
            while(!(cin >> row2 >> col2 >> slice2) || row2 < row1 || col2 < col1 || slice2 < slice1) {
                clear_cin();
				cout << "Incorrect input, enter three positive numbers with spaces: ";
			}
            cout << endl << "Current Block Temps" << endl;
            cout << setw(10) << "row" << " " << setw(10) << "col" << " " << setw(10) << "slice" << " " << setw(OUT_PRECISION+5) << "temp" << endl;
            cout << setw(10) << row1 << " " << setw(10) << col1 << " " << setw(10) << slice1 << setw(OUT_PRECISION+5) << fixed << setprecision(OUT_PRECISION) << temp[row1][col1][slice1] << endl;
            cout << setw(10) << row2 << " " << setw(10) << col2 << " " << setw(10) << slice2 << setw(OUT_PRECISION+5) << temp[row2][col2][slice2] << endl;
            cout << "Enter a New Temperature For the Block: ";
            while(!(cin >> new_temp)) {
                clear_cin();
                cout << "Incorrect input, enter a new temperature: ";
            }
            for(int i = row1; i < row2; i++) {
                for(int j = col1; j < col2; j++) {
                    for(int k = slice1; k < slice2; k++) {
                        if(i >= 0 && i < num_rows && j >= 0 && j < num_cols && k >= 0 && k < num_slices) {
                            temp[i][j][k] = new_temp;
                        }
                    }
                }
            }
        }
    }
    
    /*
     * Allows the user to start one or more moving sources.
     */
	 /*
    cout << endl << endl << "To Start One or More Moving Sources Enter 1, Else Enter 0: ";
    while(!(cin >> using_moving_source) || using_moving_source < 0 || using_moving_source > 1) {
        clear_cin();
        cout << "Incorrect Input, Enter 1 to Change, Else 0: ";
    }
    if(using_moving_source == 1) {
        REAL mag, angle1, angle2;
        cout << "Enter the number of moving sources: ";
        while(!(cin >> num_mvsrc) || num_mvsrc <= 0) {
            clear_cin();
            cout << "Incorrect input, enter a number greater than 0: ";
        }
        mvsrc_x = new REAL[num_mvsrc];
        mvsrc_y = new REAL[num_mvsrc];
        mvsrc_z = new REAL[num_mvsrc];
        mvsrc_offset_x = new REAL[num_mvsrc];
        mvsrc_offset_y = new REAL[num_mvsrc];
        mvsrc_offset_z = new REAL[num_mvsrc];
        mvsrc_vel_x = new REAL[num_mvsrc];
        mvsrc_vel_y = new REAL[num_mvsrc];
        mvsrc_vel_z = new REAL[num_mvsrc];
        mvsrc_accel_x = new REAL[num_mvsrc];
        mvsrc_accel_y = new REAL[num_mvsrc];
        mvsrc_accel_z = new REAL[num_mvsrc];
        mvsrc_temp = new REAL[num_mvsrc];
        mvsrc_valid = new int[num_mvsrc];
        for(int i = 0; i < num_mvsrc; i++) {
            cout << endl << "Moving source " << i << endl;
            cout << "Valid coordinates are x=0-"<<max_dist_x<<" y=0-"<<max_dist_y<<" z=0-"<<max_dist_z<<":" << endl;
            cout << "Enter the coordinates in meters for the corner closest to the origin, <x> <y> <z>: ";
            while(!(cin >> mvsrc_x[i] >> mvsrc_y[i] >> mvsrc_z[i]) || mvsrc_x[i] < 0 || mvsrc_x[i] > max_dist_x || mvsrc_y[i] < 0 || mvsrc_y[i] > max_dist_y || mvsrc_z[i] < 0 || mvsrc_z[i] > max_dist_z) {
                clear_cin();
                cout << "Incorrect input, enter a valid coordinate between x=0-"<<max_dist_x<<" y=0-"<<max_dist_y<<" z=0-"<<max_dist_z<<":";
            }
            cout << "Valid sizes are x=0-"<<max_dist_x-mvsrc_x[i]<<" y=0-"<<max_dist_y-mvsrc_y[i]<<" z=0-"<<max_dist_z-mvsrc_z[i]<<":"<<endl;
            cout << "Enter the size of the moving source in meters, <x size> <y size> <z size>: ";
            while(!(cin >> mvsrc_offset_x[i] >> mvsrc_offset_y[i] >> mvsrc_offset_z[i]) || mvsrc_offset_x[i] <= 0 || mvsrc_offset_x[i] > max_dist_x-mvsrc_x[i] || mvsrc_offset_y[i] <= 0 || mvsrc_offset_y[i] > max_dist_y-mvsrc_y[i] || mvsrc_offset_z[i] <= 0 || mvsrc_offset_z[i] > max_dist_z-mvsrc_z[i]) {
                clear_cin();
                cout << "Incorrect input, enter a valid distance between x=0-"<<max_dist_x-mvsrc_x[i]<<" y=0-"<<max_dist_y-mvsrc_y[i]<<" z=0-"<<max_dist_z-mvsrc_z[i]<<":";
            }
            cout << "Enter the angle of the moving sources vector in degrees from positve x towards negative y (0-360): ";
            while(!(cin >> angle1) || angle1 < 0 || angle1 > 360) {
                clear_cin();
                cout << "Incorrect input, enter a valid angle: ";
            }
            cout << "Enter the angle of the moving sources vector in degrees from positve z (0-180): ";
            while(!(cin >> angle2) || angle2 < 0 || angle2 > 180) {
                clear_cin();
                cout << "Incorrect input, enter a valid angle: ";
            }
            cout << "Enter the magnitude of the velocity vector in m/year: ";
            while(!(cin >> mag) || mag < 0) {
                clear_cin();
                cout << "Incorrect input, enter a velocity greater than 0: ";
            }
            mvsrc_vel_x[i] = mag*sin(angle2/180.0*M_PI)*cos(angle1/180.0*M_PI);
            mvsrc_vel_y[i] = mag*sin(angle2/180.0*M_PI)*sin(angle1/180.0*M_PI);
            mvsrc_vel_z[i] = mag*cos(angle2/180.0*M_PI);

            cout << "Enter the magnitude of the acceleration vector in m/year^2: ";
            while(!(cin >> mag) || mag < 0) {
                clear_cin();
                cout << "Incorrect input, enter an acceleration greater than 0: ";
            }
            mvsrc_accel_x[i] = mag*sin(angle2/180.0*M_PI)*cos(angle1/180.0*M_PI);
            mvsrc_accel_y[i] = mag*sin(angle2/180.0*M_PI)*sin(angle1/180.0*M_PI);
            mvsrc_accel_z[i] = mag*cos(angle2/180.0*M_PI);
            
            cout << "Enter the temperature of the moving source: ";
            while(!(cin >> mag)) {
                clear_cin();
                cout << "Incorrect input, enter a valid temperature: ";
            }
            mvsrc_temp[i] = mag;

            mvsrc_valid[i] = 1;

            update_mvsrc(i);
        }
    }
    */
	
    //Allows the user to decrease the size of the time step
	if(my_rank == 0) {
		cout << endl << endl << "Each Iteration in Time Spans " << scientific << time_step << " Years" << endl;
		cout << "Enter a Shorter Iteration Time in Years if Desired (any larger number otherwise): " << temp_time_step << endl;
		/*while(!(cin >> temp_val) || temp_val <= 0) {
			clear_cin();
			cout << "Incorrect input, enter a number greater than 0: ";
		}
		if(temp_val < time_step) {
			time_step = temp_val;
		}*/
    }
    
    DHF = chf * QFAC * time_step;

    //Calculates the time increment per convection loop
    time_inc = time_step/num_conv_loops;
    min_row_dim = 100.0;
    for(int i = 0; i < num_rows; i++) {
        if(dim_y[i] < min_row_dim) {
            min_row_dim = dim_y[i];
        }
    }
    
    //Asks the user for the runtime of the simulation
	if(my_rank == 0) {
		thermal_time_constant = min_row_dim*min_row_dim/max_thermal_conduct_diff;
		cout << endl << endl << "The Thermal Time Constant for the Vertical Dimension is " << thermal_time_constant << " Years" << endl;
		cout << "Enter Time Duration for Calculation in Years: " << run_time << endl;
		/*
		while(!(cin >> run_time) || run_time <= 0) {
			clear_cin();
			cout << "Incorrect input, enter a number greater than 0: ";
		}*/
		MPI_Bcast(&run_time,1,MPI_REAL_VALUE,0,MPI_COMM_WORLD);
		
		//Asks the user for the number of loops to perform between screen updates
		cout << endl << endl << "Enter the Number of Loops Between Screen Updates: " << num_loops << endl;
		/*while(!(cin >> num_loops) || num_loops <= 0) {
			clear_cin();
			cout << "Incorrect input, enter a number greater than 0: ";
		}*/
		
		MPI_Bcast(&num_loops,1,MPI_INT,0,MPI_COMM_WORLD);
		/*cout << endl << endl << "To have the simulation stop once the temperature change meets a tolerance, Enter 1 otherwise 0: ";
		while(!(cin >> use_tolerance) || use_tolerance < 0 || use_tolerance > 1) {
			clear_cin();
			cout << "Incorrect Input, Enter 1 to use a tolerance, Else 0: ";
		}*/
		use_tolerance = 0;
		MPI_Bcast(&use_tolerance,1,MPI_INT,0,MPI_COMM_WORLD);
		if(use_tolerance == 1) {
			cout << endl << "Enter the tolerance: ";
			while(!(cin >> tolerance) || tolerance  <= 0) {
				clear_cin();
				cout << "Incorrect input, enter a number greater than 0: ";
			}
		}
		MPI_Bcast(&tolerance,1,MPI_REAL_VALUE,0,MPI_COMM_WORLD);
    }
	else {
		MPI_Bcast(&run_time,1,MPI_REAL_VALUE,0,MPI_COMM_WORLD);
		MPI_Bcast(&num_loops,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(&use_tolerance,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(&tolerance,1,MPI_REAL_VALUE,0,MPI_COMM_WORLD);
	}
    //Initializes the simulation time to 0.0
    sim_time = 0.0;
    
    //Calculates the starting and ending index for each mpi process
    if(my_rank == 0) {
		st_index = 0;
		en_index = chunk_size[my_rank];
	}
	else {
		st_index = 1;
		en_index = chunk_size[my_rank]+1;
	}
    
	/*
    //Waits for the user to hit enter before beginning the simulation
    cout << endl;
    cin.ignore(numeric_limits <streamsize> ::max(), '\n' );
    PressEnterToContinue();
	*/
    /*
     * The main loop of the simulation
     */
    count = 0;    //Number of loops performed
	if(my_rank == 0) {
		cout << endl << endl << num_loops << " loops between screen updates" << endl << endl;
		if(use_tolerance == 0) {
			cout << setw(15) << "num loops" << setw(20) << "run time (years)" << setw(20) << "sim time (years)" << endl;
		}
		else {
			cout << setw(15) << "num loops" << setw(20) << "run time (years)" << setw(20) << "sim time (years)" << setw(20) << "Max temp diff" << endl;
		}
	}
#ifdef DISPLAY
	if(display_mode == 1) {
		array_minmax();
		array_size = num_cols * num_rows;
		color_field = new float[array_size * 3];
		for (int i=0; i<array_size *3; i++) {
			color_field[i] = 0.0;
		}

		glutInit(&argc, argv);
		int windowWidth = glutGet(GLUT_SCREEN_WIDTH);
		int windowHeight = glutGet(GLUT_SCREEN_HEIGHT);

		glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
		glutInitWindowSize(windowWidth, windowHeight);
		glutInitWindowPosition(0, 0);
		
		glutCreateWindow("ARC Simulation");

		glViewport(0, 0, windowWidth,windowHeight);

		glEnable (GL_BLEND);
		glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		gluPerspective(60, 1.77777f, 1.0, 20000.0);

		glutDisplayFunc(display3D);
		
		//glutMouseFunc(mouse_button);//Mouse motion and camera trans settings maintained for debugging
		//glutMotionFunc(mouse_move);
		glutKeyboardFunc(keyboard);
		glutSpecialFunc(keyboardSpecial);

		/*camera_trans[0] = -num_cols/2.0;
		camera_trans[1] = num_rows/3.0;
		camera_trans[2] = -num_rows*1.75*tan(28.0/180.0*M_PI);
		camera_rot[0] = 28.0;
		camera_trans_lag[0] = -num_cols/2.0;
		camera_trans_lag[1] = num_rows/3.0;
		camera_trans_lag[2] = -num_rows*1.75*tan(28.0/180.0*M_PI);
		camera_rot_lag[0] = 28.0;
		*/

		//gluLookAt(num_cols/2.0,num_rows*0.1,num_rows,num_cols/2.0,-num_rows/3.0,0.0,0.0,1.0,0.0);
        if(num_rows > num_cols) {
            gluLookAt(num_cols/2.0,num_rows*0.1,num_rows,num_cols/2.0,-num_rows/3.0,0.0,0.0,1.0,0.0);
        }
        else {
            gluLookAt(num_cols/2.0,num_rows*0.1,num_cols,num_cols/2.0,-num_rows/3.0,0.0,0.0,1.0,0.0);
        }
		glutMainLoop();
		//PressEnterToContinue();
	}
	else {
#endif
        MPI_Barrier(MPI_COMM_WORLD);
#ifdef TIME
        if(my_rank == 0) {
#ifdef _WIN32
            GetSystemTime(&start);
#else
            gettimeofday(&start, NULL);
#endif
        }
#endif
        while(sim_time <= run_time) {
            if(my_rank == 0) {
                //Displays status information for the current loop
                if(count%num_loops == 0) {
                    if(use_tolerance == 0) {
                        cout << setw(15) << count << setw(20) << fixed << setprecision(5) << sim_time << setw(20) << initial_time + sim_time << endl;
                    }
                    else {
                        cout << setw(15) << count << setw(20) << fixed << setprecision(5) << sim_time << setw(20) << initial_time + sim_time << setw(20) << max_temp_diff << endl;
                    }
                    //Saves the current state of the simulation if the save_state flag is set
                    /*if(save_state) {
                        save_model_state();
                    }*/
                }
            }
            
            //Performs convection updates if the current simulation is using convection
            if(using_convection) {
                convection();
            }
            
            //Performs conduction calculations
            conduction();
            //Increments the simulation time and loop count
            sim_time += time_step;
            count++;
            
            /*
            if(use_tolerance == 1) {
                max_temp_diff = find_max_temp_diff();
                if(max_temp_diff < tolerance) {
                    cout << "Maximum temperature change below the tolerance, stoping the simulation" << endl;
                    break;
                }
            }

            
            //Updates the moving source
            if(using_moving_source == 1) {
                update_moving_sources();
            }
            */
        }
        
    
        //Saves the final result of the simulation
        if(save_state == 1 || save_result == 1) {
            save_model_state();
        }
        
        save_surfer();

        if(my_rank == 0) {
            //Waits for the user to hit enter before ending the simulation
            cout << endl << "Simulation Complete" << endl;
            //PressEnterToContinue();
        }

        deallocate_memory();
       
        MPI_Barrier(MPI_COMM_WORLD);

#ifdef TIME
        if(my_rank == 0) {
#ifdef _WIN32
            GetSystemTime(&end);
            cout << "run time : " << (end.wHour-start.wHour)*3600.0 + (end.wMinute-start.wMinute)*60.0 + (end.wSecond-start.wSecond) + (end.wMilliseconds-start.wMilliseconds)/1000.0 << endl;
#else
            gettimeofday(&end, NULL);
            cout << "run time : " << (end.tv_sec  - start.tv_sec) + (end.tv_usec - start.tv_usec)/1e6 << endl;
#endif
        }
#endif
        
        //Finalizes the MPI program.
        if(my_rank == 0) {
            MPI_Finalize();
        }

        
#ifdef DISPLAY
	}
#endif
}
