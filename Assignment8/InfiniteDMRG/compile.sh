# 
# Automatic compilation of Fortran codes contained in the present folder
# 

# set paths which contain libraries (LAPACK, FFTW)
LAPACK_PATH=/usr/local/lib
#FFTW_PATH=/usr/local/lib

# remove previous files
rm -f -- *.o
rm -f -- *.out
rm -f -- *.dat
rm -f -r results
rm -f *.log

# compile single codes without linking them
gfortran -I ../modules -Wall -Wextra -O3 \
-c InfiniteDMRG.f90 \
-L$LAPACK_PATH -llapack -lblas -ltmglib

# then link and produce executable *.out
gfortran -Wall -Wextra -O3 \
-o InfiniteDMRG.out \
InfiniteDMRG.o ../modules/checkpoint_mod.o ../modules/ManyBodyUtils_mod.o \
-L$LAPACK_PATH -llapack -lblas -ltmglib

mkdir -p results
echo mmax,Nitermax,lambda,iters,tf-ti >> results/timescalings.csv