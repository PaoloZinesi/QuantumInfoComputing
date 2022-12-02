! Program to test performances of different functions.
! This program uses the subroutines "my_MatMul_naive" and "my_MatMul_opt" and 
! test their performances compared to the intrinsic function "MATMUL"
! 
! The dimensions of the matrices are given in the command-line when calling this executable
! in the order SIZE(A,1), SIZE(A,2), SIZE(B,1), SIZE(B,2). They must be integer and greater than zero.
! 
! Example of execution:
! ./a.out 100 200 200 400
! will create the matrices A(100,200) and B(200,400) and compute their matrix
! multiplication
!
! Author: Paolo Zinesi
!

! main program
PROGRAM MatMul_test_performances
      USE MatMul_mod
      USE checkpoint_mod
      IMPLICIT NONE
      LOGICAL :: debug = .TRUE.

      ! variables to read command-line arguments
      CHARACTER(LEN=32) :: arg
      INTEGER :: arg_int, arg_idx

      ! define matrices to be allocated later
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: AM
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: BM
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: CM

      ! variables to store matrix dimensions in the order
      ! SIZE(AM,1), SIZE(AM,2), SIZE(BM,1), SIZE(BM,2)
      INTEGER, DIMENSION(4) :: ABsizes
      
      ! utility and loop variables
      REAL :: ti, tf
      INTEGER :: ii, jj
      CHARACTER(LEN=1024) :: str



      ! -----------------------------------------
      ! read command-line arguments
      ! -----------------------------------------

      IF (COMMAND_ARGUMENT_COUNT() .EQ. 4) THEN
            DO arg_idx = 1, 4

                  ! get argument number "arg_idx" and put it into "arg_int"
                  CALL GET_COMMAND_ARGUMENT(arg_idx, arg)
                  READ(arg, *) arg_int

                  ABsizes(arg_idx) = arg_int

            END DO
            
      ELSE
            PRINT *, "Exactly 4 command-line arguments are expected"
            STOP
      END IF
      CALL checkpoint(debug, str="rows of A = ", val=ABsizes(1))
      CALL checkpoint(debug, str="columns of A = ", val=ABsizes(2))
      CALL checkpoint(debug, str="rows of B = ", val=ABsizes(3))
      CALL checkpoint(debug, str="columns of B = ", val=ABsizes(4))
      CALL checkpoint(debug, str="")



      ! -----------------------------------------
      ! creation and filling of matrices
      ! -----------------------------------------
      ! important checks on matrix dimensions
      ! matching dimensions
      IF (ABsizes(2) .NE. ABsizes(3)) THEN
            PRINT *, "Matrix dimensions do not match: SIZE(A,2) != SIZE(B,1)"
            STOP
      END IF
      CALL checkpoint(debug, str="Matrix dimensions match")

      ! positive dimensions
      DO ii = 1, 4
            IF(ABsizes(ii) .LE. 0) THEN
                  PRINT *, "Non-positive number given as matrix dimension"
                  STOP
            END IF
      END DO
      CALL checkpoint(debug, str="Matrix dimensions are positive")


      ! allocation of matrices (after dimension checks)
      ALLOCATE(AM(ABsizes(1), ABsizes(2)))
      ALLOCATE(BM(ABsizes(3), ABsizes(4)))
      ALLOCATE(CM(ABsizes(1), ABsizes(4)))
      CALL checkpoint(debug, str="Allocation of matrices perfomed successfully")


      

      ! fill matrix AM with numbers in [0,9] range
      DO jj = 1, SIZE(AM,2)
            DO ii = 1, SIZE(AM,1)
                  AM(ii,jj) = MOD(ii+jj, 10)
            END DO
      END DO
      ! print AM
      IF (debug .EQV. .TRUE.) THEN
            CALL checkpoint(debug, str="Matrix A: ")
            CALL print_Mat(AM)
      END IF

      ! fill matrix BM with numbers in [0,14] range
      DO jj = 1, SIZE(BM,2)
            DO ii = 1, SIZE(BM,1)
                  BM(ii,jj) = MOD(ii+jj, 15)
            END DO
      END DO
      ! print BM
      IF (debug .EQV. .TRUE.) THEN
            CALL checkpoint(debug, str="Matrix B: ")
            CALL print_Mat(BM)
      END IF


      ! -----------------------------------------
      ! start of the function testing part
      ! -----------------------------------------

      ! naive function testing
      CALL CPU_TIME(ti)
            CALL my_MatMul_naive(AM,BM,CM)
      CALL CPU_TIME(tf)

      ! print performances string
      WRITE (str, "(I2,I2,I2,I2,A10,ES15.8E2)") &
            SIZE(AM,1), SIZE(AM,2), SIZE(BM,1), SIZE(BM,2), " ; naive; ", tf-ti
      PRINT *, TRIM(str)
      ! print "my_MatMul_naive" result
      IF (debug .EQV. .TRUE.) THEN
            CALL checkpoint(debug, str="Matrix product with 'my_MatMul_naive':")
            CALL print_Mat(CM)
      END IF



      ! optimized function testing
      CALL CPU_TIME(ti)
            CALL my_MatMul_opt(AM,BM,CM)
      CALL CPU_TIME(tf)

      ! print performances string
      WRITE (str, "(I2,I2,I2,I2,A8,ES15.8E2)") &
            SIZE(AM,1), SIZE(AM,2), SIZE(BM,1), SIZE(BM,2), " ; opt; ", tf-ti
      PRINT *, TRIM(str)
      ! print "my_MatMul_opt" result
      IF (debug .EQV. .TRUE.) THEN
            CALL checkpoint(debug, str="Matrix product with 'my_MatMul_opt':")
            CALL print_Mat(CM)
      END IF


      ! Fortran intrinsic function testing
      CALL CPU_TIME(ti)
            CM = MATMUL(AM,BM)
      CALL CPU_TIME(tf)

      ! print performances string
      WRITE (str, "(I2,I2,I2,I2,A12,ES15.8E2)") &
            SIZE(AM,1), SIZE(AM,2), SIZE(BM,1), SIZE(BM,2), " ; builtin; ", tf-ti
      PRINT *, TRIM(str)
      ! print INTRINSIC MATMUL result
      IF (debug .EQV. .TRUE.) THEN
            CALL checkpoint(debug, str="Matrix product with intrinsic MATMUL function:")
            CALL print_Mat(CM)
      END IF

      

      ! -----------------------------------------
      ! final deallocation of matrices
      ! -----------------------------------------
      DEALLOCATE(AM)
      DEALLOCATE(BM)
      DEALLOCATE(CM)
      CALL checkpoint(debug, str="Deallocation of matrices perfomed successfully")
      
      STOP
END PROGRAM MatMul_test_performances