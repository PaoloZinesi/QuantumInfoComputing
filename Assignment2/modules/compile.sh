# 
# Automatic compilation of Fortran codes contained in the present folder
# 

# remove previous files
rm -f -- *.o
rm -f -- *.mod

# compile checkpoint module first
gfortran -Wall -Wextra \
-c checkpoint_mod.f90

# compile all the modules
gfortran -Wall -Wextra \
-c *.f90