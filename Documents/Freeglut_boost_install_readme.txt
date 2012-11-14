Instructions for installing and compiling on MS Visual Studio (I'm assuming that this will work similarly under the free VS Express version):

1. Download the freeglut packages for windows from http://www.transmissionzero.co.uk/software/freeglut-devel/. (look for latest MSVC zip file unless you will be compiling via MinGW on another IDE.)
2. Extract to a common directory available to all users. This is not normally compliant with Windows installation conventions, but c:\program files\ ought to be fine.
3. Likewise download the latest boost libraries from http://www.boost.org/users/download/. (.zip file is fine.)
4. Extract as before to a common directory like c:\program files\.
5. We now need to reference the headers and libraries called by the #include section in the program.
6. Right click the project title and select properties. Change 'Configuration' at the top to 'all' so that the changes will affect both debug and release.
7. Under configuration properties, select C/C++. In the 'Additional Include Directories' at right add two lines directing the compiler to both your boost directory (for example: c:\Program Files\boost\boost_version\) and the freeglut include directory (for example: c:\Program Files\Freeglut\include\).
8. Also under configuration properties, select Linker. Under 'General', change 'Additional Library Directories' at right to include your freeglut library directory (for example: c:\Program Files\Freeglut\lib). Then under 'Advanced', change 'Entry Point' at right to 'mainCRTStartup' to initialize C runtime library for a console program. 

Notes: 
- Again, we acknowledge that installation to the program files directory is not considered best practice. Copying to c:\ or some other directory should be equally valid.
- Also, if any command line parameters need to be defined, especially when running the debugger (F5), they can be specified under the project properties under Configuration Properties->Debugging->Command Arguments.