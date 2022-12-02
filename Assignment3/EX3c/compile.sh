# 
# Automatic compilation of Fortran codes contained in the present folder
# 

# set path which contains lapack libraries (liblapack.a, libblas.a, librefblas.a, libtmglib.a)
LAPACK_PATH=/usr/local/lib

# remove previous files
rm -f -- a.out
rm -f -- *.o
rm -f -- *.dat
rm -f -r results
rm -f errors.log

# compile single codes without linking them
gfortran -I ../modules -Wall -Wextra -O3 \
-c EX3c_Zinesi_CODE.f90

# then link and produce executable a.out
gfortran -Wall -Wextra -O3 \
EX3c_Zinesi_CODE.o ../modules/checkpoint_mod.o ../modules/DCmatrix_mod.o \
-L$LAPACK_PATH -llapack -lblas -ltmglib