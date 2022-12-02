

# Before using any code in the directories compile the 
# "modules" folder using the bash script "modules/compile.sh"
cd modules ; bash compile.sh ; cd ..


# After all the .mod and .o files have been generated in the "modules"
# folder, proceed with the compilation of source codes in every subfolder
# "EX2*" using the bash scripts "EX2*/compile.sh"
# 
# Example: 
cd EX2a ; bash compile.sh


# Finally, the executable "a.out" will be find in the folder
./a.out


# Examples of execution for each program are contained in the 
# documentation at the beginning of every "*_Zinesi_CODE.f90"
# source files and in the "execute_example.txt" files


