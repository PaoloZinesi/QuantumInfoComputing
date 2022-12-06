! Program to test the FFTW3 routines and the Zwavefunc module.
! The tests are performed by defining a wavefunction in psi_x for which
! the Fourier Transform (FT) is known analytically. Then, psi_x is FT to
! psi_p and finally transformed back to psi_x, gaining a global multiplicative
! factor of NN as explained in the documentation of FFTW3.
! 
! Example of execution:
! ./fftw_test
! will create three files in the folder 'results', which can be plotted using 
! EX5_plots.ipynb notebook. In these output files the first row contains the grid
! on which the WF is defined, while the other rows contain the values of the WF
! (real and imaginary part are the two values that appear in each row).
!
! Author: Paolo Zinesi
!

! main program
PROGRAM FFTW_test
      USE checkpoint_mod
      USE, INTRINSIC :: iso_c_binding
      USE FFTW_mod
      USE Zwavefunc_mod
      IMPLICIT NONE
      LOGICAL :: debug = .TRUE.

      REAL*8, PARAMETER :: PI = 4.D0*ATAN(1.D0)
      INTEGER :: ii
      REAL*4 :: ti, tf

      type(C_PTR) :: plan_direct, plan_inverse
      INTEGER :: NN = 500
      REAL*8 :: LL = 1.D0
      REAL*8 :: sigmax


      ! WF declarations and allocations
      TYPE(Zwavefunc) :: psi_x, psi_p
      psi_x = Zwavefunc(length=NN, need_fftw_alloc=.TRUE.)
      psi_p = Zwavefunc(length=NN, need_fftw_alloc=.TRUE.)

      ! create plans before filling the arrays
      CALL CPU_TIME(ti)
      plan_direct = fftw_plan_dft_1d(NN, psi_x%elem_fftw,psi_p%elem_fftw, FFTW_FORWARD,FFTW_MEASURE)
      plan_inverse = fftw_plan_dft_1d(NN, psi_p%elem_fftw,psi_x%elem_fftw, FFTW_BACKWARD,FFTW_MEASURE)
      CALL CPU_TIME(tf)
      CALL checkpoint(debug, str="Time [s] necessary to produce both plans =", val=tf-ti)

      ! fill array grids
      CALL create_grid(psi_x, xmin=0.D0, xmax=LL)
      CALL create_grid(psi_p, xmin=0.D0, xmax=(2*PI*NN)/LL)

      ! shift the grid of p in order to be compatible with both FFTW and physical
      ! conventions on the reciprocal-space array psi_p
      psi_p%grid(FLOOR((NN - 1)/2.D0)+2:NN) = psi_p%grid(FLOOR((NN - 1)/2.D0)+2:NN) - (2*PI*NN)/LL

      ! fill array
      sigmax = LL/20.D0
      psi_x%elem_fftw = (/ (  EXP(-0.5D0*((psi_x%grid(ii) - (LL/2.4D0))/sigmax)**2.D0) * &
                              EXP(COMPLEX(0,(2*PI/LL)*10*psi_x%grid(ii))), ii=1,NN) /)


      CALL checkpoint(debug, str="")


      ! write wavefunction into file
      CALL SYSTEM('mkdir -p results')
      CALL writeWFFile(psi_x, unit=10, file="results/init.dat", format="ES24.17")
      CALL checkpoint(debug, str="Norm^2 of psi_x =", val=norm2_WF(psi_x))
      CALL checkpoint(debug, str="")


      ! direct FFT
      CALL CPU_TIME(ti)
      CALL fftw_execute_dft(plan_direct, psi_x%elem_fftw, psi_p%elem_fftw)
      CALL CPU_TIME(tf)
      CALL checkpoint(debug, str="Time [s] necessary to execute direct plan =", val=tf-ti)

      CALL writeWFFile(psi_p, unit=20, file="results/FFT_init.dat", format="ES24.17")
      CALL checkpoint(debug, str="Norm^2 of psi_p =", val=norm2_WF(psi_p))
      CALL checkpoint(debug, str="")


      ! inverse FFT
      CALL CPU_TIME(ti)
      CALL fftw_execute_dft(plan_inverse, psi_p%elem_fftw, psi_x%elem_fftw)
      CALL CPU_TIME(tf)
      CALL checkpoint(debug, str="Time [s] necessary to execute inverse plan =", val=tf-ti)

      CALL writeWFFile(psi_x, unit=30, file="results/FFT_inverse.dat", format="ES24.17")
      CALL checkpoint(debug, str="Norm^2 of psi_x (after FT and inverse FT) =", val=norm2_WF(psi_x))
      CALL checkpoint(debug, str="")


      ! destroy plans and free WF memory
      CALL fftw_destroy_plan(plan_direct)
      CALL fftw_destroy_plan(plan_inverse)
      CALL freeWF(psi_x)
      CALL freeWF(psi_p)

END PROGRAM FFTW_test