! Module that contains various matrix utilities, in particular two
! matrix multiplication subroutines.
! 
! Author: Paolo Zinesi
!

MODULE MatMul_mod
      USE checkpoint_mod
      IMPLICIT NONE

      LOGICAL, PRIVATE :: debug = .FALSE.
      ! The interface in this module is NOT needed, because FORTRAN
      ! generates a dedicated one when importing this module in the main program

      CONTAINS


      ! This subroutine implements for INTEGER*4 matrices the "naive"  
      ! matrix multiplication, defined as
      ! C_{i,k} = sum_{j} [A_{i,j} * B{j,k}] with
      ! A --> leftM
      ! B --> rightM
      ! C --> outM
      ! The index of summation (j) runs along the columns of A and along the rows of B
      ! 
      ! inputs:
      ! - leftM [integer(:,:)]: left matrix in the multiplication
      ! - rightM [integer(:,:)]: right matrix in the multiplication
      ! - outM [integer(:,:)]: matrix where to store the multiplication results
      !
      ! outputs: None
      ! 
      ! TODO: None
      SUBROUTINE my_MatMul_naive(leftM, rightM, outM)
            IMPLICIT NONE
            INTEGER, DIMENSION(:,:), INTENT(IN) :: leftM
            INTEGER, DIMENSION(:,:), INTENT(IN) :: rightM
            INTEGER, DIMENSION(:,:), INTENT(OUT) :: outM
      
            ! utility and loop variables
            INTEGER :: ii, kk, jj
            INTEGER :: left_rowsize, right_colsize, internal_size
            CHARACTER(LEN=1024) :: str
      
            left_rowsize = SIZE(leftM,1)
            right_colsize = SIZE(rightM,2)
            internal_size = SIZE(leftM,2)
      
            CALL checkpoint(debug, str=NEW_LINE("")//"Entering 'my_MatMul_naive' subroutine")

            ! print dimensions debug string
            WRITE (str, "(A12,I2,I2,A13,I2,I2,A11,I2,I2)") &
                  "DIM(leftM) =", SIZE(leftM,1), SIZE(leftM,2), &
                  "; DIM(rightM) =", SIZE(rightM,1), SIZE(rightM,2), &
                  "; DIM(outM) =", SIZE(outM,1), SIZE(outM,2)
            CALL checkpoint(debug, str=TRIM(str))


            ! check for correctness of the matrix dimensions
            ! if dimensions are not correct print an error message
            IF (  (SIZE(leftM,1) .EQ. SIZE(outM,1)) .AND. &
                  (SIZE(leftM,2) .EQ. SIZE(rightM,1)) .AND. &
                  (SIZE(rightM,2) .EQ. SIZE(outM,2))) THEN
                  
                  CALL checkpoint(debug, str="'my_MatMul_naive': Matrix dimensions match")

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
                  PRINT *, "'my_MatMul_naive': Matrix dimensions do NOT match"

            ENDIF
            CALL checkpoint(debug, str="Exiting 'my_MatMul_naive' subroutine"//NEW_LINE(""))
      END SUBROUTINE my_MatMul_naive


      ! This subroutine implements for INTEGER*4 matrices the optimized version of 
      ! my_MatMul_naive subroutine.
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
      ! - leftM [integer(:,:)]: left matrix in the multiplication
      ! - rightM [integer(:,:)]: right matrix in the multiplication
      ! - outM [integer(:,:)]: free matrix where to store the multiplication results
      !
      ! outputs: None
      ! 
      ! TODO: None
      SUBROUTINE my_MatMul_opt(leftM, rightM, outM)
            IMPLICIT NONE
            INTEGER, DIMENSION(:,:), INTENT(IN) :: leftM
            INTEGER, DIMENSION(:,:), INTENT(IN) :: rightM
            INTEGER, DIMENSION(:,:), INTENT(OUT) :: outM
      
            ! transposed matrix used in the multiplication
            INTEGER, DIMENSION(:,:), ALLOCATABLE :: leftM_T
            
            ! utility and loop variables
            INTEGER :: ii, kk, jj
            INTEGER :: left_rowsize, right_colsize, internal_size
            CHARACTER(LEN=1024) :: str
      
            left_rowsize = SIZE(leftM,1)
            right_colsize = SIZE(rightM,2)
            internal_size = SIZE(leftM,2)
      
            CALL checkpoint(debug, str=NEW_LINE("")//"Entering 'my_MatMul_opt' subroutine")
            ALLOCATE(leftM_T(SIZE(leftM,2),SIZE(leftM,1)))

            ! print dimensions debug string
            WRITE (str, "(A12,I2,I2,A13,I2,I2,A11,I2,I2)") &
                  "DIM(leftM) =", SIZE(leftM,1), SIZE(leftM,2), &
                  "; DIM(rightM) =", SIZE(rightM,1), SIZE(rightM,2), &
                  "; DIM(outM) =", SIZE(outM,1), SIZE(outM,2)
            CALL checkpoint(debug, str=TRIM(str))


            ! check for correctness of the matrix dimensions
            ! if dimensions are not correct print an error message
            IF (  (SIZE(leftM,1) .EQ. SIZE(outM,1)) .AND. &
                  (SIZE(leftM,2) .EQ. SIZE(rightM,1)) .AND. &
                  (SIZE(rightM,2) .EQ. SIZE(outM,2))) THEN
                  
                  CALL checkpoint(debug, str="'my_MatMul_opt': Matrix dimensions match")
      
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
                  PRINT *, "'my_MatMul_opt': Matrix dimensions do NOT match"

            ENDIF

            DEALLOCATE(leftM_T)
            CALL checkpoint(debug, str="Exiting 'my_MatMul_opt' subroutine"//NEW_LINE(""))
      END SUBROUTINE my_MatMul_opt


      ! print matrix M_{i,j} on the standard output (*),
      ! row by row in a readable form
      ! 
      ! inputs:
      ! - M [integer(:,:)]: matrix to print
      !
      ! outputs: None
      !
      ! TODO: None
      !
      SUBROUTINE print_Mat(M)
            IMPLICIT NONE
            INTEGER, DIMENSION(:,:), INTENT(IN) :: M
      
            INTEGER :: ii
      
            ! cycle over all rows and print each of them
            PRINT *
            DO ii = 1, SIZE(M,1)
                  PRINT *, M(ii,:)
            END DO
            PRINT *
      END SUBROUTINE print_Mat



      ! TODO:
      ! create a function that returns .TRUE. if two matrices are equal and .FALSE.
      ! otherwise, for this simple algorithm it is not necessary but it might be 
      ! useful in the future



END MODULE MatMul_mod