! My first program in Fortran.
! This simple program is used to check if the Fortran code can compile.
! 
! Author: Paolo Zinesi
!

PROGRAM hello_world
      IMPLICIT NONE

      INTEGER :: N = 2

      IF (N .EQ. 2) THEN
            PRINT *, "Hello, world!"
      ELSE
            PRINT *, "Something's wrong here..."
      ENDIF

      STOP
END PROGRAM hello_world