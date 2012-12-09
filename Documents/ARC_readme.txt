ARC Heat Transfer Readme


***I. Running ARC***

Options are similar to the original Fortran program and those that are different are self-explanatory. Controls used under visualization are the up/down arrow keys and '+' and '-' keys. Full mouse control code is available but commented out in the source. You'll need to do some manipulation to get this working again if you want to use it, for example for further visualization development.


***II. Compiling on a Windows System using Microsoft Visual Studio Express***

A. Installing a compiler/IDE

You will need to install and register the Visual Studio IDE (not free or VS Express (free). On Windows 7/8, Visual Studio 2012 is current. For XP VS 2008 is recommended.
http://www.microsoft.com/en-us/download/details.aspx?id=6506 - v2008 (vcsetup.exe)
http://www.microsoft.com/visualstudio/eng/products/visual-studio-express-for-windows-desktop - v2012

Aanother option is to use Eclipse C++ in combination with minGW. The compiler is a port of Linux's gcc compiler and could make porting between operating systems easer. I'll give some relevant links, but only cover VS Express for now, assuming that is the more popular option.
http://www.eclipse.org/downloads/moreinfo/c.php
http://help.eclipse.org/indigo/index.jsp?topic=/org.eclipse.cdt.doc.user/concepts/cdt_c_before_you_begin.htm
http://www.mingw.org/


B. FreeGlut

1. Download the freeglut packages for windows from http://www.transmissionzero.co.uk/software/freeglut-devel/. (look for latest MSVC zip file unless you will be compiling via MinGW on another IDE.)
2. Extract to a common directory available to all users. Usually something close to the root directory should be fine like c:\freeglut\.
3. Make a copy of freeglut32.dll to either c:\Windows\System32\ (32 bit operating system, more likely if using XP) or c:\Windows\SysWOW64 (64 bit operating system, more likely if using Windows 7). This will allow Visual Studio and Windows to "see" the freeglut link library.


C. VSE Project Setup

Visual Studio likes organizing things into Solutions and Projects. Solutions are big overarching, well, 'projects' and Projects are sub-groups that encompass associated code in, say, a suite of products. Each project has its own code and own dependencies and header files. When you begin, you'll need to create a new project and solution. For ARC, it will need to be a  Win32 console application. Name both as you wish. In the wizard that pops up, select 'Console application' and check 'empty project' since we're inputting our own code, not using pre-generated stuff from Microsoft. This opens the project. 

Copy the arc.cpp file into your (My Documents\Visual Studio\Projects\'ProjectName'\'ProjectName'\) folder. If you want to use any input files during debugging, put those in here as well. Now right-click 'Source Files' on the Solution Explorer to the left. Add existing item and choose the source file we just copied. Now you have access to the program. NOTE: When working with any other program that has a header file, you'll need to add those just as you did your source files.

We now need to reference headers and libraries in the VS solution so that the #include section will call the proper files. Right click the project title and select properties. Change 'Configuration' at the top to 'all' so that the changes will affect both debug and release. Under configuration properties, select C/C++. In the 'Additional Include Directories' at right add two lines directing the compiler the freeglut include directory (for example: c:\freeglut\include\). Also under configuration properties, select Linker. Under 'General', change 'Additional Library Directories' at right to include your freeglut library directory (for example: c:\freeglut\lib). Then under 'Advanced', change 'Entry Point' at right to 'mainCRTStartup' to initialize C runtime library for a console program. Also under configuration properties, select Linker. Under 'General', change 'Additional Library Directories' at right to include your freeglut library directory (for example: c:\freeglut\lib). (NOTE: Older versions of Visual studio may have these settings in slightly different places. Their location should still be predictable and findable though.)


D. Running the program

Options should be self explanatory and very similar to the Fortran version. Any new options should be self-evident such as the option to hold maximum temperature on the color scale constant or let it update. If you want to debug, use 'CTRL-F5'. This will create a temporary build from source and give you feedback in the IDE, including appropriate compile errors line by line in the log below. If you want to create a release build, change the dropdown at the top menu from 'Debug' to 'Release'. Select Build->Build solution (or press F7) This will create an executable you can use day-to-day in whatever location you want. It will be created in your (My Documents\Visual Studio\Projects\'ProjectName'\Release\) folder.

NOTE: If you want to make any command-line arguments (for example -DDISPLAY to enable the simulation display's #ifdef checks), they can be specified under the project properties under Configuration Properties->Linker->Command Arguments.

When you want to have a second project in the solution (for example the conversion and comparison programs), you'll need to handle them a little differently. Right-click the solution title on the left and select Add->New Project and follow the wizard as before. When you want to debug a different project, you'll have to right-click the project and select 'Use as Startup Project' or it won't build. Doing a release build will compile all projects in the solution.

If you have multiple source files in your project that are unrelated, you'll want to exclude all but the one you want compiled in the build. Do this by right-clicking, going into the file's properties you *don't* want to compile and changing 'excluded from build' to 'yes'. Do not right-click and select 'exclude from project'. This will simply remove it from the project (clear as mud, no?).




***III. Running Under Linux***

A. Compiling for Linux can be in some ways easier, in some ways a little more difficult. The first thing you'll need to do is make sure freeGlut and OpenGL are installed. In a terminal command-line do the following in a debian based distribution like Ubuntu or Mint:
	1) sudo apt-get update 
	2) sudo apt-get install build-essential (these commands update the package repositories and make sure you have everything you need to compile the code)
	3) sudo apt-get install freeglut3-dev (installs all the libraries and headers you need for freeGlut)
	
Compile commands will usually look something like: 'g++ sourcefile.cpp -o sourcefile' (or gcc for pure c code), -o creates the interim object file before going to the linker to create the final binary executable. You'll need to use other command line switches at certain times. freeGlut will require the -lglut switch. #ifdef DISPLAY lines in the ARC source will require you to add a -DDISPLAY to compile the display code. Also there are certain optimization flags that will help your program run faster once debugging is done. -o3 will include most common optimizations.

B. To run in linux or unix, usually the command './filename' will do the trick.

C. To edit the source, you will need a text editor. Distributions commonly come with one installed that will highlight your code easily, such as Gedit in Ubuntu. Working further with Travis and his students you will also commonly see Vim being used, a very powerful text editor with an admittedly steep learning curve.


Other useful links:
http://www.boost.org/users/download/ - We used Boost extensively in our class and it is fairly joined at the hip to the standard C++ libraries in terms of acceptance. It has some data structures and containers that might be useful, but most of the time you might find its random number generators and distribution functions the most useful. Download the latest boost libraries from http://www.boost.org/users/download/(.zip file is fine.). Extract as before to a common directory like c:\boost\.
http://paulbourke.net/texture_colour/colourspace/ - A lengthy discussion of scientific color progressions. I used the 'Jet' style progression that uses the 24 bit RGB extents. It also has .dat files of 256 values you may find useful for future work.
http://netbeans.org/features/index.html - Alternative NetBeans IDE. Will do C++ but requires Cygwin to compile (See below).
http://www.cygwin.com/ - Alternative to Mingw and MSVC compilers. Useful if you want a unix/linux like command line environment in windows along with several standard compilers and programs like gcc and vim. Used by the NetBeans IDE.
