# 
# Automatic compilation of Fortran codes contained in the present folder
# 

# set paths which contain libraries (LAPACK, FFTW)
#LAPACK_PATH=/usr/local/lib
FFTW_PATH=/usr/local/lib

# remove previous files
rm -f -- *.o
rm -f -- *.dat
rm -f -r results
rm -f *.log

# compile single codes without linking them
gfortran -I ../modules -Wall -Wextra -O3 \
-Wno-unused-parameter \
-c fftw_test.f90 TD_HarmOsc_1D.f90 \
-L$FFTW_PATH -lfftw3 -lm

# link and produce executable 'fftw_test.out'
gfortran -Wall -Wextra -O3 \
-o fftw_test.out \
fftw_test.o ../modules/checkpoint_mod.o ../modules/Zwavefunc_mod.o \
-L$FFTW_PATH -lfftw3 -lm

# link and produce executable 'TD_HarmOsc_1D.out'
gfortran -Wall -Wextra -O3 \
-o TD_HarmOsc_1D.out \
TD_HarmOsc_1D.o ../modules/checkpoint_mod.o ../modules/Zwavefunc_mod.o \
-L$FFTW_PATH -lfftw3 -lm