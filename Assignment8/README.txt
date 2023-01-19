

# Before using any code in the directories compile the 
# "modules" folder using the bash script "modules/compile.sh"
cd modules ; bash compile.sh ; cd ..


# After all the .mod and .o files have been generated in the "modules"
# folder, proceed with the compilation of source codes in the folder
# using the bash scripts "*/compile.sh"
# 
# Example: 
cd RealSpaceRG ; bash compile.sh


# Finally, the executable "RealSpaceRG.out" will be find in the folder "RealSpaceRG"


# Examples of execution for each program are contained in the 
# documentation at the beginning of every "*.f90"
# source files and in the "execute_example.txt" files



# In order to use the LAPACK and BLAS libraries in the code, the variable
# $LAPACK_PATH is set in each 'compile.sh' script. This path needs to be
# changed manually if the libraries are in a different folder.
# Same for the FFTW3 library.