! Program to test the limits of INTEGER and REAL in fortran.
! Errors/warnings in the compilation of this code are the main purposes of this program.
! 
! Author: Paolo Zinesi
!

PROGRAM test_type_limits
      IMPLICIT NONE

      ! declarations
      INTEGER*2 :: i2add1
      INTEGER*2 :: i2add2
      INTEGER*4 :: i4add1
      INTEGER*4 :: i4add2

      REAL*4 :: r4add1
      REAL*4 :: r4add2
      REAL*8 :: r8add1
      REAL*8 :: r8add2


      ! gfortran compiler finds an arithmetic overflow here
      ! we can force the compilation by using the flag ‘-fno-range-check’
      i2add1 = 2000000
      i2add2 = 1
      i4add1 = 2000000
      i4add2 = 1

      ! gfortran compiler finds a type conversion from 'REAL(8)' to 'REAL(4)'
      ! only when using the flag ‘-Wconversion’, but it compiles anyway
      r4add1 = 4.D32*DATAN(1.D0)
      r4add2 = 1.D21*SQRT(2.D0)
      r8add1 = 4.D32*DATAN(1.D0)
      r8add2 = 1.D21*SQRT(2.D0)

      PRINT *, "Result of the sum with INTEGER*2:", i2add1+i2add2
      PRINT *, "Result of the sum with INTEGER*4:", i4add1+i4add2            

      PRINT *, "Result of the sum with REAL*4:", r4add1+r4add2
      PRINT *, "Result of the sum with REAL*8:", r8add1+r8add2

      STOP
END PROGRAM test_type_limits