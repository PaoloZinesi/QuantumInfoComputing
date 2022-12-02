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
-c HarmOsc_1D.f90 \
-L$LAPACK_PATH -llapack -lblas -ltmglib

# then link and produce executable a.out
gfortran -Wall -Wextra -O3 \
HarmOsc_1D.o ../modules/checkpoint_mod.o \
-L$LAPACK_PATH -llapack -lblas -ltmglib