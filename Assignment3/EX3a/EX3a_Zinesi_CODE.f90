! Program to test performances of different functions.
! This program uses the subroutines "my_MatMul_naive" and "my_MatMul_opt" and 
! test their performances compared to the intrinsic function "MATMUL"
! 
! The dimensions of the matrices are given in the command-line when calling this executable
! in the order SIZE(A,1), SIZE(A,2), SIZE(B,1), SIZE(B,2). They must be integer and greater than zero.
! Next, in the command line the method of multiplication should be specified. Accepted strings are: {'naive','opt','builtin'}
! 
! Example of execution:
! ./a.out 100 200 200 400 naive
! will create the matrices A(100,200) and B(200,400) and compute their matrix
! multiplication using the 'my_MatMul_naive' subroutine
!
! Author: Paolo Zinesi
!

! main program
PROGRAM MatMul_test_performances
      USE MatMul_mod
      USE checkpoint_mod
      IMPLICIT NONE
      LOGICAL :: debug = .FALSE.

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

      ! variable to store chosen multiplication method
      CHARACTER(LEN=32) :: matmul_method
      
      ! utility and loop variables
      REAL :: ti, tf
      INTEGER :: ii, jj
      CHARACTER(LEN=1024) :: str



      ! -----------------------------------------
      ! read command-line arguments
      ! -----------------------------------------

      IF (COMMAND_ARGUMENT_COUNT() .EQ. 5) THEN
            DO arg_idx = 1, 4

                  ! get argument number "arg_idx" and put it into "arg_int"
                  CALL GET_COMMAND_ARGUMENT(arg_idx, arg)
                  READ(arg, *) arg_int

                  ABsizes(arg_idx) = arg_int

            END DO

            CALL GET_COMMAND_ARGUMENT(5, matmul_method)
            
      ELSE
            ! write error details into an error log file 
            PRINT *, "An error occurred, details on 'errors.log'"
            OPEN(unit=10, file='errors.log', access='APPEND')
                  WRITE(10, "('Exactly 5 command-line arguments are expected')")
            CLOSE(unit=10)
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
            ! write error details into an error log file
            PRINT *, "An error occurred, details on 'errors.log'"
            OPEN(unit=10, file='errors.log', access='APPEND')
                  WRITE(10, "('Matrix dimensions do not match: SIZE(A,2) != SIZE(B,1)')")
            CLOSE(unit=10)
            STOP
      END IF
      CALL checkpoint(debug, str="Matrix dimensions match")

      ! positive dimensions
      DO ii = 1, 4
            IF(ABsizes(ii) .LE. 0) THEN
                  ! write error details into an error log file 
                  PRINT *, "An error occurred, details on 'errors.log'"
                  OPEN(unit=10, file='errors.log', access='APPEND')
                        WRITE(10, "('Non-positive number given as matrix dimension')")
                  CLOSE(unit=10)
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
                  AM(ii,jj) = MOD(ii+3*jj, 10)
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
                  BM(ii,jj) = MOD(ii+7*jj, 15)
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

      ! create directory where to store results
      CALL SYSTEM('mkdir -p results')


      IF (TRIM(matmul_method) .EQ. "naive") THEN
            ! naive function testing
            CALL CPU_TIME(ti)
                  CALL my_MatMul_naive(AM,BM,CM)
            CALL CPU_TIME(tf)
      

            ! print performances string
            WRITE (str, "(I5,'; ',I5,'; ',I5,'; ',I5,'; ',ES15.8E2)") &
                  SIZE(AM,1), SIZE(AM,2), SIZE(BM,1), SIZE(BM,2), tf-ti

            ! write data into a dedicated file 
            OPEN(unit=10, file='results/'//TRIM(matmul_method)//'_performances.dat', access='APPEND')
                  WRITE(10, '(A)') TRIM(str)
            CLOSE(unit=10)


            ! print "my_MatMul_naive" result
            IF (debug .EQV. .TRUE.) THEN
                  CALL checkpoint(debug, str="Matrix product with 'my_MatMul_naive':")
                  CALL print_Mat(CM)
            END IF

      ELSE IF (TRIM(matmul_method) .EQ. "opt") THEN
            ! optimized function testing
            CALL CPU_TIME(ti)
                  CALL my_MatMul_opt(AM,BM,CM)
            CALL CPU_TIME(tf)
      

            ! print performances string
            WRITE (str, "(I5,'; ',I5,'; ',I5,'; ',I5,'; ',ES15.8E2)") &
                  SIZE(AM,1), SIZE(AM,2), SIZE(BM,1), SIZE(BM,2), tf-ti

            ! write data into a dedicated file 
            OPEN(unit=20, file='results/'//TRIM(matmul_method)//'_performances.dat', access='APPEND')
                  WRITE(20, '(A)') TRIM(str)
            CLOSE(unit=20)


            ! print "my_MatMul_opt" result
            IF (debug .EQV. .TRUE.) THEN
                  CALL checkpoint(debug, str="Matrix product with 'my_MatMul_opt':")
                  CALL print_Mat(CM)
            END IF

      ELSE IF (TRIM(matmul_method) .EQ. "builtin") THEN
            ! Fortran intrinsic function testing
            CALL CPU_TIME(ti)
                  CM = MATMUL(AM,BM)
            CALL CPU_TIME(tf)
      

            ! print performances string
            WRITE (str, "(I5,'; ',I5,'; ',I5,'; ',I5,'; ',ES15.8E2)") &
                  SIZE(AM,1), SIZE(AM,2), SIZE(BM,1), SIZE(BM,2), tf-ti

            ! write data into a dedicated file 
            OPEN(unit=30, file='results/'//TRIM(matmul_method)//'_performances.dat', access='APPEND')
                  WRITE(30, '(A)') TRIM(str)
            CLOSE(unit=30)


            ! print INTRINSIC MATMUL result
            IF (debug .EQV. .TRUE.) THEN
                  CALL checkpoint(debug, str="Matrix product with intrinsic MATMUL function:")
                  CALL print_Mat(CM)
            END IF

      ELSE
            ! write error details into an error log file 
            PRINT *, "An error occurred, details on 'errors.log'"
            OPEN(unit=40, file='errors.log', access='APPEND')
                  WRITE(40, "('No function corresponding to `', A, '`')") TRIM(matmul_method)
            CLOSE(unit=40)
            STOP

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