! Module that contains various many body Hamiltonians utilities, such as subroutines for left and right tensor products
! with identities for many-body Hamiltonians. These subroutines can then be applied to simulations with very large Hamiltonians
! (such as RG, DMRG,...).
! 
! Author: Paolo Zinesi
!

MODULE ManyBodyUtils_mod
      USE checkpoint_mod
      IMPLICIT NONE

      LOGICAL, PRIVATE :: debug = .FALSE.
      ! The interface in this module is NOT needed, because FORTRAN
      ! generates a dedicated one when importing this module in the main program

      CONTAINS


      ! This subroutine fills an input matrix of dimension (2**k,2**k) with sub-matrices of the form
      ! D=delta*id_{2**(k-2)} on the 4 antidiagonal blocks. The typical form of the filled matrix AA is
      !     AA =  0 0 0 D 
      !           0 0 D 0
      !           0 D 0 0
      !           D 0 0 0
      !
      ! inputs:
      ! - AA [REAL*8, dimension(:,:)]: input REAL*8 (sub)matrix to fill as explained above
      ! - k [integer]: base-2 logarithm of matrix size, such that 2**k = SIZE(AA,1) = SIZE(AA,2)
      ! - delta [REAL*8]: values used to fill the antidiagonal matrix
      !                              
      ! outputs: None
      ! 
      ! TODO: None
      !
      SUBROUTINE fillAntidiag(AA, k, delta)
            ! matrix to fill
            REAL*8, DIMENSION(:,:), INTENT(INOUT) :: AA

            ! 2**k = size(AA,1)
            INTEGER, INTENT(IN) :: k
            REAL*8, OPTIONAL :: delta
            REAL*8 :: delta_actual

            ! loop variables
            INTEGER :: ii

            ! check if square matrix
            IF(SIZE(AA,1) .NE. SIZE(AA,2)) THEN
                  PRINT "('Not square matrix')"
                  STOP
            END IF

            ! check if 2**k = size(AA,1)
            IF(SIZE(AA,1) .NE. 2**k) THEN
                  PRINT "('Dimensions are not correct')"
                  PRINT "('size of AA = ', I0, ', ', I0, ' != ', I0)", SIZE(AA,1), SIZE(AA,2), 2**k
                  STOP
            END IF

            ! check if matrix is at least 4x4
            IF(k .LT. 2) THEN
                  PRINT "('The matrix size cannot be smaller than 4x4')"
                  STOP
            END IF

            ! set optional value for delta
            IF(PRESENT(delta)) THEN
                  delta_actual = delta
            ELSE
                  delta_actual = 1.D0
            END IF


            ! -------------------------
            ! filling of submatrices
            ! -------------------------

            ! leftermost identity sub-matrix
            DO ii = 1, 2**(k-2)
                  AA(3*(2**(k-2))+ii, ii) = AA(3*(2**(k-2))+ii, ii) + delta_actual
            END DO

            ! second leftermost identity sub-matrix
            DO ii = 1, 2**(k-2)
                  AA(2*(2**(k-2))+ii, 1*(2**(k-2))+ii) = AA(2*(2**(k-2))+ii, 1*(2**(k-2))+ii) + delta_actual
            END DO

            ! second rightermost identity sub-matrix
            DO ii = 1, 2**(k-2)
                  AA(1*(2**(k-2))+ii, 2*(2**(k-2))+ii) = AA(1*(2**(k-2))+ii, 2*(2**(k-2))+ii) + delta_actual
            END DO

            ! rightermost identity sub-matrix
            DO ii = 1, 2**(k-2)
                  AA(ii, 3*(2**(k-2))+ii) = AA(ii, 3*(2**(k-2))+ii) + delta_actual
            END DO

      END SUBROUTINE fillAntidiag



      ! This subroutine fills an input matrix of dimension (2**NN,2**NN) with a lambda*sigma_i^z Hamiltonian
      ! (1 <= i <= NN) of the subsystem 'i' over all the total NN subsystems.
      !
      ! inputs:
      ! - MM [REAL*8, dimension(:,:)]: input REAL*8 matrix to fill with lambda*sigma_i^z
      ! - ii [integer]: index of subsystem to compute the sigma_z Hamiltonian of (1 <= ii <= NN)
      ! - NN [integer]: base-2 logarithm of matrix size, such that 2**NN = SIZE(MM,1) = SIZE(MM,2)
      ! - lambda [REAL*8]: multiplicative constant in front of the Hamiltonian term
      !                              
      ! outputs: None
      ! 
      ! TODO: None
      !
      SUBROUTINE fillSigmaZ(MM, ii, NN, lambda)
            ! matrix to fill
            REAL*8, DIMENSION(:,:), INTENT(INOUT) :: MM

            ! parameters to fill Hamiltonian
            INTEGER, INTENT(IN) :: NN
            INTEGER, INTENT(IN) :: ii
            REAL*8, OPTIONAL :: lambda
            REAL*8 :: lambda_actual


            ! loop variables
            INTEGER :: jj

            ! check if square matrix
            IF(SIZE(MM,1) .NE. SIZE(MM,2)) THEN
                  PRINT "('Not square matrix')"
                  STOP
            END IF

            ! check if 2**NN = size(MM,1)
            IF(SIZE(MM,1) .NE. 2**NN) THEN
                  PRINT "('Dimensions are not correct')"
                  PRINT "('size of MM = ', I0, ', ', I0, ' != ', I0)", SIZE(MM,1), SIZE(MM,2), 2**NN
                  STOP
            END IF

            ! check if value of ii is valid
            IF((ii .LT. 1) .OR. (ii .GT. NN)) THEN
                  PRINT "('Index of subsystem for sigmaZ is not valid')"
                  STOP
            END IF

            ! set optional value for lambda
            IF(PRESENT(lambda)) THEN
                  lambda_actual = lambda
            ELSE
                  lambda_actual = 1.D0
            END IF

            
            ! -------------------------
            ! filling of sigma_z
            ! -------------------------

            DO jj = 1, 2**NN

                  ! exploit the alternating signs for the sigma_z Hamiltonians
                  ! to easily compute the sigma_z contributions
                  MM(jj,jj) = MM(jj,jj) + lambda_actual * ((-1)**((jj-1)/(2**(NN-ii))))
            END DO
            CALL checkpoint(debug, "SigmaZ contribution filled successfully")


      END SUBROUTINE fillSigmaZ



      ! This subroutine fills an input matrix of dimension (2**NN,2**NN) with a lambda*sigma_i^x*sigma_{i+1}^x
      ! Hamiltonian (1 <= i <= NN-1) of the interaction between subsystems 'i' and 'i+1' over all the total 
      ! NN subsystems.
      !
      ! inputs:
      ! - MM [REAL*8, dimension(:,:)]: input REAL*8 matrix to fill with lambda*sigma_i^x*sigma_{i+1}^x
      ! - ii [integer]: index of first subsystem to compute the sigma_x sigma_x interaction (1 <= ii <= NN-1)
      ! - NN [integer]: base-2 logarithm of matrix size, such that 2**NN = SIZE(MM,1) = SIZE(MM,2)
      ! - lambda [REAL*8]: multiplicative constant in front of the Hamiltonian term
      !                              
      ! outputs: None
      ! 
      ! TODO: None
      !
      SUBROUTINE fillSigmaXSigmaX(MM, ii, NN, lambda)
            ! matrix to fill
            REAL*8, DIMENSION(:,:), INTENT(INOUT) :: MM

            ! parameters to fill Hamiltonian
            INTEGER, INTENT(IN) :: NN
            INTEGER, INTENT(IN) :: ii
            REAL*8, OPTIONAL :: lambda
            REAL*8 :: lambda_actual

            ! loop variables
            INTEGER :: bb

            ! check if square matrix
            IF(SIZE(MM,1) .NE. SIZE(MM,2)) THEN
                  PRINT "('Not square matrix')"
                  STOP
            END IF

            ! check if 2**NN = size(MM,1)
            IF(SIZE(MM,1) .NE. 2**NN) THEN
                  PRINT "('Dimensions are not correct')"
                  PRINT "('size of MM = ', I0, ', ', I0, ' != ', I0)", SIZE(MM,1), SIZE(MM,2), 2**NN
                  STOP
            END IF

            ! check if value of ii is valid
            IF((ii .LT. 1) .OR. (ii .GT. NN-1)) THEN
                  PRINT "('Index of subsystem for sigmaX sigmaX is not valid')"
                  STOP
            END IF

            ! set optional value for lambda
            IF(PRESENT(lambda)) THEN
                  lambda_actual = lambda
            ELSE
                  lambda_actual = 1.D0
            END IF

            
            ! --------------------------------
            ! filling of sigma_x sigma_x
            ! --------------------------------

            ! fill an antidiaognal matrix for all the 2**(ii-1) diagonal blocks
            ! to obtain the sigma_x sigma_x n.n. interaction
            DO bb = 1, 2**(ii-1)
                  CALL fillAntidiag(MM((bb-1)*2**(NN+1-ii)+1:bb*2**(NN+1-ii), &
                                       (bb-1)*2**(NN+1-ii)+1:bb*2**(NN+1-ii)), &
                                    NN+1-ii, lambda_actual)
            END DO
            CALL checkpoint(debug, "SigmaX SigmaX contribution filled successfully")


      END SUBROUTINE fillSigmaXSigmaX



      ! This subroutine fills an input matrix of dimension (2**NN,2**NN) with a id_{Nleft} subMM id_{Nright} matrix.
      ! The size of the two identities and of subMM must match to give the size of MM. Nright is inferred from 
      ! NN, subNN and Nleft (Nright = NN - subNN - Nleft). The submatrix subMM can be of the same size of MM, 
      ! in that case Nleft must be = 0 (also Nright will be = 0) and the martix MM will be simply filled with subMM.
      ! This subroutine is optimized for small values of subNN, while it performs loosely (FOR NOW) on bigger tensor products.
      !
      ! inputs:
      ! - MM [REAL*8, dimension(:,:)]: input REAL*8 matrix to fill with tensor product of identities and subMM
      ! - subMM [REAL*8, dimension(:,:)]: input REAL*8 matrix to be used in the tensor product with identities
      ! - NN [integer]: base-2 logarithm of MM size, such that 2**NN = SIZE(MM,1) = SIZE(MM,2)
      ! - subNN [integer]: base-2 logarithm of subMM size, such that 2**subNN = SIZE(subMM,1) = SIZE(subMM,2)
      ! - Nleft [integer]: dimension of the left identity
      !                              
      ! outputs: None
      ! 
      ! TODO: Improve the tensor product computation for large subMM
      !
      SUBROUTINE tensorProductIdentity(MM, subMM, NN, subNN, Nleft)
            ! matrices
            REAL*8, DIMENSION(:,:), INTENT(INOUT) :: MM, subMM

            ! parameters to fill matrix
            INTEGER, INTENT(IN) :: NN, subNN, Nleft
            INTEGER :: Nright

            ! loop variables
            INTEGER :: leftB, rightB, rowsub, colsub
            INTEGER :: row_i, col_i

            ! check if square matrices
            IF(SIZE(MM,1) .NE. SIZE(MM,2)) THEN
                  PRINT "('Matrix to fill is not a square matrix')"
                  STOP
            END IF
            IF(SIZE(subMM,1) .NE. SIZE(subMM,2)) THEN
                  PRINT "('Matrix for tensor product is not a square matrix')"
                  STOP
            END IF

            ! check if 2**NN = size(MM,1)
            IF(SIZE(MM,1) .NE. 2**NN) THEN
                  PRINT "('Dimensions are not correct')"
                  PRINT "('size of MM = ', I0, ', ', I0, ' != ', I0)", SIZE(MM,1), SIZE(MM,2), 2**NN
                  STOP
            END IF
            ! check if 2**subNN = size(subMM,1)
            IF(SIZE(subMM,1) .NE. 2**subNN) THEN
                  PRINT "('Dimensions are not correct')"
                  PRINT "('size of subMM = ', I0, ', ', I0, ' != ', I0)", SIZE(subMM,1), SIZE(subMM,2), 2**subNN
                  STOP
            END IF

            ! check if value of dimensions is correct to perform the tensor product
            IF((NN .LT. 1) .OR. (subNN .LT. 1) .OR. (subNN .GT. NN)) THEN
                  PRINT "('Sizes of matrices are not valid')"
                  STOP
            END IF
            IF((Nleft .LT. 0) .OR. (NN .LT. subNN+Nleft)) THEN
                  PRINT "('Sizes of identity matrices for the tensor products are not valid')"
                  STOP
            END IF


            ! dimension of right identity matrix
            Nright =  NN - subNN - Nleft
            

            ! --------------------------------
            ! filling of tensor product
            ! --------------------------------

            ! do the filling for each element of the left identity matrix
            DO leftB = 0, (2**Nleft)-1

                  ! fix a row and column value for the subMM matrix
                  DO colsub = 1, 2**subNN
                  DO rowsub = 1, 2**subNN

                        ! do the filling for each element of the right identity matrix
                        DO rightB = 0, (2**Nright)-1

                              ! compute multi-index
                              row_i = leftB*2**(subNN+Nright) + (rowsub-1)*2**Nright + rightB + 1
                              col_i = leftB*2**(subNN+Nright) + (colsub-1)*2**Nright + rightB + 1

                              ! fill element
                              MM(row_i, col_i) = MM(row_i, col_i) + subMM(rowsub, colsub)

                        END DO
                  END DO
                  END DO
            END DO

            CALL checkpoint(debug, "Tensor products with identities filled successfully")


      END SUBROUTINE tensorProductIdentity



      ! This subroutine fills an input matrix of dimension (2**(sizeAA+sizeBB),2**(sizeAA+sizeBB)) with
      ! AA * BB tensor product of square matrices.
      ! The size of AABB must match the input sizes sizeAA and sizeBB of AA and BB, respectively.
      ! This subroutine is the most general subroutine to compute the tensor product of two input matrices
      ! AA and BB, and for this reason it has poor performances when AA, BB have many zeroes (like identities).
      ! For more specific tasks, other subroutines need to be used.
      !
      ! inputs:
      ! - AABB [REAL*8, dimension(:,:)]: input REAL*8 matrix to fill with tensor product of AA and BB
      ! - AA [REAL*8, dimension(:,:)]: input REAL*8 matrix to be used in the tensor product (on the left)
      ! - BB [REAL*8, dimension(:,:)]: input REAL*8 matrix to be used in the tensor product (on the right)
      ! - sizeAA [integer]: base-2 logarithm of AA size, such that 2**sizeAA = SIZE(AA,1) = SIZE(AA,2)
      ! - sizeBB [integer]: base-2 logarithm of BB size, such that 2**sizeBB = SIZE(BB,1) = SIZE(BB,2)
      !                              
      ! outputs: None
      ! 
      ! TODO: None
      !
      SUBROUTINE generalTensorProduct(AABB, AA, BB, sizeAA, sizeBB)
            ! matrices
            REAL*8, DIMENSION(:,:), INTENT(INOUT) :: AABB, AA, BB

            ! parameters to fill tensor product
            INTEGER, INTENT(IN) :: sizeAA, sizeBB

            ! loop variables
            INTEGER :: rowA, colA, rowAB, colAB

            ! check if square matrices
            IF(SIZE(AABB,1) .NE. SIZE(AABB,2)) THEN
                  PRINT "('Matrix to fill is not a square matrix')"
                  STOP
            END IF
            IF(SIZE(AA,1) .NE. SIZE(AA,2)) THEN
                  PRINT "('Left matrix for tensor product is not a square matrix')"
                  STOP
            END IF
            IF(SIZE(BB,1) .NE. SIZE(BB,2)) THEN
                  PRINT "('Right matrix for tensor product is not a square matrix')"
                  STOP
            END IF

            ! check if 2**(sizeAA+sizeBB) = size(AABB,1)
            IF(SIZE(AABB,1) .NE. 2**(sizeAA+sizeBB)) THEN
                  PRINT "('Dimensions are not correct')"
                  PRINT "('size of AABB = ', I0, ', ', I0, ' != ', I0)", SIZE(AABB,1), SIZE(AABB,2), 2**(sizeAA+sizeBB)
                  STOP
            END IF
            ! check if 2**sizeAA = size(AA,1)
            IF(SIZE(AA,1) .NE. 2**sizeAA) THEN
                  PRINT "('Dimensions are not correct')"
                  PRINT "('size of AA = ', I0, ', ', I0, ' != ', I0)", SIZE(AA,1), SIZE(AA,2), 2**sizeAA
                  STOP
            END IF
            ! check if 2**sizeAA = size(AA,1)
            IF(SIZE(BB,1) .NE. 2**sizeBB) THEN
                  PRINT "('Dimensions are not correct')"
                  PRINT "('size of BB = ', I0, ', ', I0, ' != ', I0)", SIZE(BB,1), SIZE(BB,2), 2**sizeBB
                  STOP
            END IF

            ! check if value of dimensions is correct to perform the tensor product
            IF((sizeAA .LT. 1) .OR. (sizeBB .LT. 1)) THEN
                  PRINT "('Sizes of matrices are not valid')"
                  STOP
            END IF
            

            ! --------------------------------
            ! filling of tensor product
            ! --------------------------------

            ! for each value of the left matrix AA compute the tensor product with the matrix BB
            DO colA = 1, 2**sizeAA
            DO rowA = 1, 2**sizeAA
                  ! position of first element of the submatrix corresponding to (rowA,colA)
                  rowAB = (rowA-1)*(2**sizeBB)+1
                  colAB = (colA-1)*(2**sizeBB)+1

                  ! fill entire submatrix in a single instruction
                  AABB(rowAB:rowAB+2**sizeBB-1, colAB:colAB+2**sizeBB-1) = &
                        AABB(rowAB:rowAB+2**sizeBB-1, colAB:colAB+2**sizeBB-1) + BB*AA(rowA,colA)
                        
            END DO
            END DO

            CALL checkpoint(debug, "Tensor products filled successfully")


      END SUBROUTINE generalTensorProduct
      


      ! This subroutine prints on a file a REAL*8 matrix M in readable form.
      ! Unit, filename, and output format are given as input for flexibility.
      ! Only units different from 0,5,6 can be used to write into file (otherwise 
      ! collisions with default writing units arise).
      !
      ! inputs:
      ! - M [REAL*8, dimension(:,:)]: input REAL*8 to write into a file
      ! - unit [integer]: identifier of the unit connected to the file
      ! - file [character(len=*)]: filename where matrix has to be written
      ! - format [character(len=*)]: string to specify how to format written numbers.
      !                              A suggestion is to use format="(ES24.17)"
      !                              
      ! outputs: None
      ! 
      ! TODO: None
      ! 
      SUBROUTINE writeMatFile(M, unit, file, format)
            IMPLICIT NONE

            REAL*8, DIMENSION(:,:), INTENT(IN) :: M
            INTEGER :: unit
            CHARACTER (LEN=*) :: file, format

            ! utility and loop variables
            INTEGER :: ii, jj

            ! do not use units dedicated to standard fortran I/O (PRINT * for ex.)
            IF((unit .EQ. 0) .OR. (unit .EQ. 5) .OR. (unit .EQ. 6)) THEN
                  PRINT "('Not allowed to use unit', I2)", unit
                  PRINT "(A, ' has not been written')", file

            ELSE
                  ! write data into the file 
                  OPEN(unit=unit, file=file)
                        DO ii = 1, SIZE(M,1)
                              DO jj = 1, SIZE(M,2)
                                    ! advance='no' to remain in the present line                        
                                    WRITE(unit, format, advance="no") M(ii,jj)
                              END DO

                              ! new line
                              WRITE(unit,*)
                        END DO
                  CLOSE(unit=unit)
            END IF
            CALL checkpoint(debug, "Written matrix on file "//file)

      END SUBROUTINE writeMatFile




END MODULE ManyBodyUtils_mod