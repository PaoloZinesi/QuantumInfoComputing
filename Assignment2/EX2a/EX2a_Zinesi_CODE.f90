! Program to test a checkpoint subroutine.
! The checkpoint subroutine takes as input a boolean variable to switch on the
! debugging mode, an optional string to print and an optional variable to print.
! 
! Example of execution:
! ./a.out
! will print the output of all the checkpoints 
!
! Author: Paolo Zinesi
!


! main program
PROGRAM checkpoint_test
      USE checkpoint_mod
      IMPLICIT NONE
      LOGICAL :: debug = .TRUE.

      CALL checkpoint(debug, str="Start checkpoints debugging"//NEW_LINE(""))
      CALL checkpoint(debug, val=0.E0)
      CALL checkpoint(debug, val=1)
      CALL checkpoint(debug, val=2.E0)

      ! integers
      CALL checkpoint(debug, val=500_2, str="INTEGER*2: ") ! INTEGER*2
      CALL checkpoint(debug, val=500000000, str="INTEGER*4: ") ! INTEGER*4
      CALL checkpoint(debug, val=50000000000_8, str="INTEGER*8: ") ! INTEGER*8

      ! reals
      CALL checkpoint(debug, val=1.E10, str="REAL*4: ") ! REAL*4
      CALL checkpoint(debug, val=1.1D200, str="REAL*8: ") ! REAL*8

      ! complex numbers
      CALL checkpoint(debug, val=COMPLEX(1.E0,1.E0), str="COMPLEX: ") ! COMPLEX
      CALL checkpoint(debug, val=COMPLEX(1.D0,1.D0), str="DOUBLE COMPLEX: ") ! DOUBLE COMPLEX

      ! unsupported types
      ! The error raised is "There is no specific subroutine for the generic ‘checkpoint’"
      ! CALL checkpoint(debug, val=.TRUE., str="LOGIC") ! LOGIC
      ! CALL checkpoint(debug, val="string", str="CHAR*") ! CHAR*
      
      CALL checkpoint(debug, str=NEW_LINE("")//"Stop checkpoints debugging")
      STOP
END PROGRAM checkpoint_test
