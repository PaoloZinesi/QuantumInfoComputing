# 
# Automatic compilation of Fortran codes contained in the present folder
# 

# remove previous files
rm -f -- a.out
rm -f errors.log
rm -f -- *.o
rm -f -r results

# compile single codes without linking them
gfortran -I ../modules -Wall -Wextra -O3 \
-c EX3a_Zinesi_CODE.f90

# then link and produce executable a.out
gfortran -Wall -Wextra -O3 \
EX3a_Zinesi_CODE.o ../modules/checkpoint_mod.o ../modules/MatMul_mod.o