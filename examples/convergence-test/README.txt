If you've read the main README and learned some basics about SymPhas, you
might be interested in getting some example projects set up.

This project will run convergence testing using the Allen Cahn manufactured
minimal solution model.

To compile:
    First, make a new directory `build'. Enter the directory and use the 
    following cmake command:

> cmake -DFFTW3_DIR:PATH="~/fftw3/lib/cmake/fftw3" -DUSE_IO:BOOL="True" -DUSE_VTK:BOOL="False" -DMODEL_INCLUDE_HEADER_NAME:FILEPATH="modelacmms.h" -DMODEL_INCLUDE_HEADER_DIR:PATH="../models" -DSOLVER_INCLUDE_HEADER_DIR:PATH="../solvers" -DSOLVER_INCLUDE_HEADER_NAME:FILEPATH="solverinclude.h" -DUSE_LATEX_FORMAT:BOOL="False" -DPRINTABLE_EQUATIONS:BOOL="True" -DCMAKE_BUILD_TYPE=Release ..

    Then run `make'. Keep in mind that you shouldn't be running any other
    exmaples at the same time within the same project heirarchy, since the 
    cmake command will change some defines that are used in the project 
    files. 

In order to run this project:
    Pass the configuration file "SymPhas/examples/config/allen_cahn_mms.in".
    Assuming that you're in the build directory, this can simply be given as
    a relative directory: "../../config/allen_cahn_mms.in". It will run about
    10 tests for changing the time step and the spatial step. It will generate
    a csv with the relevant information.


