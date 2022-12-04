# 
# Automatic compilation of Fortran codes contained in the present folder
# 

# remove previous files
rm -f -- *.o
rm -f -- *.mod

# compile checkpoint and FFTW modules first
gfortran -Wall -Wextra -O3 \
-c checkpoint_mod.f90 FFTW_mod.f90

# compile all the modules
gfortran -Wall -Wextra -O3 \
-ffree-line-length-196 \
-c *.f90