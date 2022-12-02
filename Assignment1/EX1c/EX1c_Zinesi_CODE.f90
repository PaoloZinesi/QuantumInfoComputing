! Program to test performances of different functions.
! This program uses the subroutines "my_MatMul_naive" and "my_MatMul_opt" and 
! test their performances compared to the intrinsic function "MATMUL"
!     
! Author: Paolo Zinesi
!

! main program
PROGRAM MatMul_test_performances
      IMPLICIT NONE

      REAL :: ti, tf
      INTEGER :: ii, jj

      ! defines matrix dimensions at compile time
      ! TODO in a future version: allocate matrices at runtime to decide their size
      ! without hard-coding them
      INTEGER, DIMENSION(1024,1024) :: AM
      INTEGER, DIMENSION(1024,1024) :: BM
      INTEGER, DIMENSION(1024,1024) :: CM


      INTERFACE
            SUBROUTINE my_MatMul_naive(leftM, rightM, outM)
                  INTEGER, DIMENSION(:,:), INTENT(IN) :: leftM
                  INTEGER, DIMENSION(:,:), INTENT(IN) :: rightM
                  INTEGER, DIMENSION(:,:), INTENT(OUT) :: outM
            END SUBROUTINE my_MatMul_naive

            SUBROUTINE my_MatMul_opt(leftM, rightM, outM)
                  INTEGER, DIMENSION(:,:), INTENT(IN) :: leftM
                  INTEGER, DIMENSION(:,:), INTENT(IN) :: rightM
                  INTEGER, DIMENSION(:,:), INTENT(OUT) :: outM
            END SUBROUTINE my_MatMul_opt

            SUBROUTINE print_Mat(M)
                  INTEGER, DIMENSION(:,:), INTENT(IN) :: M
            END SUBROUTINE print_Mat
      END INTERFACE


      ! fill matrix AM with numbers in [0,9] range
      DO jj = 1, SIZE(AM,2)
            DO ii = 1, SIZE(AM,1)
                  AM(ii,jj) = MOD(ii+jj, 10)
            END DO
      END DO

      ! fill matrix BM with numbers in [0,14] range
      DO jj = 1, SIZE(BM,2)
            DO ii = 1, SIZE(BM,1)
                  BM(ii,jj) = MOD(ii+jj, 15)
            END DO
      END DO


      ! -----------------------------------------
      ! start of the function testing part
      ! -----------------------------------------

      ! naive function testing
      CALL CPU_TIME(ti)
            CALL my_MatMul_naive(AM,BM,CM)
      CALL CPU_TIME(tf)

      ! print debug string
      PRINT *, SIZE(AM,1), ";naive;", tf-ti



      ! optimized function testing
      CALL CPU_TIME(ti)
            CALL my_MatMul_opt(AM,BM,CM)
      CALL CPU_TIME(tf)

      ! print debug string
      PRINT *, SIZE(AM,1), ";opt;", tf-ti



      ! Fortran intrinsic function testing
      CALL CPU_TIME(ti)
            CM = MATMUL(AM,BM)
      CALL CPU_TIME(tf)

      ! print debug string
      PRINT *, SIZE(AM,1), ";builtin;", tf-ti

      
      STOP
END PROGRAM MatMul_test_performances



! This subroutine implements the "naive" matrix multiplication defined as
! C_{i,k} = sum_{j} [A_{i,j} * B{j,k}] with
! A --> leftM
! B --> rightM
! C --> outM
! The index of summation (j) runs along the columns of A and along the rows of B
! 
! inputs:
! - leftM: left matrix in the multiplication
! - rightM: right matrix in the multiplication
! - outM: free matrix where to store the multiplication results
!
! outputs: None
! 
! TODO:
! implement an assignment operator to return the result of the multiplication.
! This requires to define an interface ...
SUBROUTINE my_MatMul_naive(leftM, rightM, outM)
      IMPLICIT NONE
      INTEGER, DIMENSION(:,:), INTENT(IN) :: leftM
      INTEGER, DIMENSION(:,:), INTENT(IN) :: rightM
      INTEGER, DIMENSION(:,:), INTENT(OUT) :: outM


      INTEGER :: ii, kk, jj
      INTEGER :: left_rowsize, right_colsize, internal_size

      left_rowsize = SIZE(leftM,1)
      right_colsize = SIZE(rightM,2)
      internal_size = SIZE(leftM,2)

      ! check for correctness of the matrix dimensions
      ! if dimensions are not correct print an error message
      IF (  (SIZE(leftM,1) .EQ. SIZE(outM,1)) .AND. &
            (SIZE(leftM,2) .EQ. SIZE(rightM,1)) .AND. &
            (SIZE(rightM,2) .EQ. SIZE(outM,2))) THEN

            ! compute explicitly the matrix multiplication
            ! and store the result into the matrix "outM"
            DO ii = 1, left_rowsize
                  DO kk = 1, right_colsize

                        outM(ii,kk) = 0
                        DO jj = 1, internal_size
                              outM(ii,kk) = outM(ii,kk) + (leftM(ii,jj)*rightM(jj,kk))
                        END DO
                  END DO
            END DO

      ELSE
            PRINT *, "Dimensions are not correct"
      ENDIF
END SUBROUTINE my_MatMul_naive



! This subroutine implements the optimized version of my_MatMul_naive subroutine.
! Since Fortran stores matrices by column, a consistent speedup is obtained if the 
! index of summation runs along the rows of both matrices. This implies that this
! subroutine must take care of the matrix transposition before computing the actual
! multiplication. However, the large speedup obtained motivates this additional 
! computational cost.
! 
! The summation ordering to speed up the multiplication is
! C_{i,k} = sum_{j} [A^T_{j,i} * B{j,k}] with
! A --> leftM
! B --> rightM
! C --> outM
! where A_T is the transposed matrix of A
! 
! inputs:
! - leftM: left matrix in the multiplication
! - rightM: right matrix in the multiplication
! - outM: free matrix where to store the multiplication results
!
! outputs: None
! 
! TODO:
! same as "my_MatMul_naive" TODO
SUBROUTINE my_MatMul_opt(leftM, rightM, outM)
      IMPLICIT NONE
      INTEGER, DIMENSION(:,:), INTENT(IN) :: leftM
      INTEGER, DIMENSION(:,:), INTENT(IN) :: rightM
      INTEGER, DIMENSION(:,:), INTENT(OUT) :: outM

      ! transposed matrix used in the multiplication
      INTEGER, DIMENSION(SIZE(leftM,2),SIZE(leftM,1)) :: leftM_T

      INTEGER :: ii, kk, jj
      INTEGER :: left_rowsize, right_colsize, internal_size

      left_rowsize = SIZE(leftM,1)
      right_colsize = SIZE(rightM,2)
      internal_size = SIZE(leftM,2)


      ! check for correctness of the matrix dimensions
      ! if dimensions are not correct print an error message
      IF (  (SIZE(leftM,1) .EQ. SIZE(outM,1)) .AND. &
            (SIZE(leftM,2) .EQ. SIZE(rightM,1)) .AND. &
            (SIZE(rightM,2) .EQ. SIZE(outM,2))) THEN


            ! calculate the transposed matrix
            DO ii = 1, internal_size
                  DO jj = 1, left_rowsize
                        leftM_T(ii,jj) = leftM(jj,ii)
                  END DO
            END DO


            ! compute explicitly the matrix multiplication
            ! and store the result into the matrix "outM"
            DO kk = 1, right_colsize
                  DO ii = 1, left_rowsize

                        outM(ii,kk) = 0
                        DO jj = 1, internal_size
                              outM(ii,kk) = outM(ii,kk) + (leftM_T(jj,ii)*rightM(jj,kk))
                        END DO
                  END DO
            END DO

      ELSE
            PRINT *, "Dimensions are not correct"
      ENDIF
END SUBROUTINE my_MatMul_opt



! print matrix M_{i,j} on the standard output (*),
! column by column as stored in memory
! 
! inputs:
! - leftM: left matrix in the multiplication
! - rightM: right matrix in the multiplication
! - outM: free matrix where to store the multiplication results
!
! outputs: 
! - printed matrix on *
!
! TODO:
! rearrange the matrix print in a more readable form
SUBROUTINE print_Mat(M)
      IMPLICIT NONE
      INTEGER, DIMENSION(:,:), INTENT(IN) :: M

      INTEGER :: ii, jj

      ! cycle over all columns and print each of them
      DO jj = 1, SIZE(M,2)
            DO ii = 1, SIZE(M,1)
                  PRINT *, M(ii,jj), ' '
            END DO
            PRINT *, "End of column", jj
      END DO
      PRINT *, "End of matrix"
END SUBROUTINE print_Mat

! TODO:
! create a function that returns .TRUE. if two matrices are equal and .FALSE.
! otherwise, for this simple algorithm it is not necessary but it might be 
! useful in the future