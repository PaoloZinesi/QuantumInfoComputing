# 
# Automatic compilation of Fortran codes contained in the present folder
# 

# remove previous files
rm -f -- a.out
rm -f -- *.o
rm -f -- *.dat

# compile single codes without linking them
gfortran -I ../modules -Wall -Wextra \
-c EX2c_Zinesi_CODE.f90

# then link and produce executable a.out
gfortran -Wall -Wextra \
EX2c_Zinesi_CODE.o ../modules/checkpoint_mod.o ../modules/DCmatrix_mod.o