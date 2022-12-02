! Module to easily import the FFTW3 subroutines without generating conflicts.
! 
! Author: Paolo Zinesi
!

MODULE FFTW_mod
      USE, INTRINSIC :: iso_c_binding
      IMPLICIT NONE

      INCLUDE "fftw3.f03"

END MODULE FFTW_mod

      