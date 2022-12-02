! Program to study the normalized spacings distribution of eigenvalues of different
! random matrices.
! This program uses the following subroutines to study different matrices:
! - DCmatrix 'RandMat' to generate random Hermitian matrices with complex entries uniformly 
! distributed in {-1 < Re(z),Im(z) < 1} in the complex plane
! - 'eigenvaluesMat' subroutine to compute the eigenvalues of random Hermitian matrices
! - 'DLARNV' subroutine to generate the diagonal of a real diagonal matrix, which
! eigenvalues are equivalent to the  diagonal entries.
!
! The normalized spacings between eigenvalues are stored in an output file depending
! on the type of random matrix considered. The (linear) dimension and the type of 
! the random square matrix to be generated is given in the command-line.
! The linear size must be greater than 1 to compute at least one eigenvalue spacing,
! while accepted strings for random matrix type are: {'HS', 'DS'} for, respectively,
! Hermitian matrices with uniform real and immaginary parts in [-1,1] and
! real diagonal matrices with uniform entries in [-1,1].
! 
! Examples of execution:
! ./a.out 1000 HS
! will create a 1000x1000 random Hermitian matrix, compute its eigenvalues and append the 
! normalized spacings in 'results/norm_spac_HS.dat'.
!
! ./a.out 1000 DS
! will create a 1000x1000 random real diagonal matrix and append the normalized spacings 
! in 'results/norm_spac_DS.dat'.
!
! Author: Paolo Zinesi
!

! main program
PROGRAM eigen_test
      USE checkpoint_mod
      USE DCmatrix_mod
      IMPLICIT NONE
      LOGICAL :: debug = .FALSE.


      ! variables to read command-line arguments
      CHARACTER(LEN=32) :: arg
      INTEGER :: arg_idx

      ! define matrices and vectors to be allocated
      INTEGER :: size
      TYPE(DCmatrix) :: AM_H
      REAL*8, DIMENSION(:), ALLOCATABLE :: eigenvals, norm_spacings
      REAL*8 :: mean_spac
      CHARACTER(LEN=4) :: mattype

      
      ! utility and loop variables
      CHARACTER(LEN=1024) :: str
      INTEGER :: info

      ! seeds variables
      INTEGER, DIMENSION(4) :: iseed
      REAL*4, DIMENSION(4) :: rseed
      

      ! -----------------------------------------
      ! read command-line arguments
      ! -----------------------------------------

      IF (COMMAND_ARGUMENT_COUNT() .EQ. 2) THEN

            ! fill "size"
            arg_idx = 1
            CALL GET_COMMAND_ARGUMENT(arg_idx, arg)
            READ(arg, *) size

            ! fill "mattype"
            arg_idx = 2
            CALL GET_COMMAND_ARGUMENT(arg_idx, mattype)
            
      ELSE
            ! write error details into an error log file 
            PRINT *, "An error occurred, details on 'errors.log'"
            OPEN(unit=10, file='errors.log', access='APPEND')
                  WRITE(10, "('Exactly 2 command-line arguments are expected')")
            CLOSE(unit=10)
            STOP
      END IF


      ! check if at least N > 1
      IF(size .LE. 1) THEN
            ! write error details into an error log file 
            PRINT *, "An error occurred, details on 'errors.log'"
            OPEN(unit=10, file='errors.log', access='APPEND')
                  WRITE(10, "('N > 1 is requires to compute normalized spacing')")
            CLOSE(unit=10)
            STOP
      END IF
      CALL checkpoint(debug, str="Linear dimension of A = ", val=size)
      CALL checkpoint(debug, str="Chosen matrix type = "//TRIM(mattype)//NEW_LINE(""))
      

      ! -----------------------------------------
      ! allocations of matrices/vectors
      ! -----------------------------------------
      
      ! allocations of complex matrix
      AM_H = DCmatrix(rows=size, cols=size)

      ! allocation of eigenvalue and normalized eigenvalue spacing arrays
      ALLOCATE(norm_spacings(size-1))
      ALLOCATE(eigenvals(size))

      CALL checkpoint(debug, str="Allocations perfomed successfully")

      ! create directory where to store results
      CALL SYSTEM('mkdir -p results')


      ! -----------------------------------------
      ! eigendecomposition
      ! -----------------------------------------

      ! initialize random seed (iseed(4) is required to be odd by LAPACK routines)
      CALL RANDOM_NUMBER(rseed)
      iseed(1:3) = INT(rseed(1:3)*4096)
      iseed(4) = 2 * INT(rseed(4)*2048) + 1

      ! Hermitian matrix
      IF (TRIM(mattype) .EQ. 'HS') THEN

            ! filling of random DCmatrix
            ! random complex numbers are generated such that real and imaginary parts belong to [-1,1] 
            ! the diagonal of the hermitian matrix is real and uniformly distributed in [-1,1]
            CALL RandMat(AM_H, dist='S', iseed=iseed, sym='H')
            CALL checkpoint(debug, str="Random matrix successfully filled.")


            ! call eigenvalues subroutine to get sorted eigenvalues
            CALL eigenvaluesMat(AM_H, eigenvals)
            CALL checkpoint(debug, str="Eigenvalues successfully calculated.")


      ! Real diagonal matrix
      ELSE IF (TRIM(mattype) .EQ. 'DS') THEN

            ! filling of random diagonal using LAPACK routine DLARN
            ! random real numbers are generated uniformly in [-1,1]
            CALL DLARNV(1, iseed, size, eigenvals)
            CALL checkpoint(debug, str="Random real diagonal successfully filled.")

            CALL DLASRT('I', size, eigenvals, info)
            CALL checkpoint(debug, str="Eigenvalues successfully calculated. INFO value = ", val=info)

      ELSE
            ! write error details into an error log file 
            PRINT *, "An error occurred, details on 'errors.log'"
            OPEN(unit=10, file='errors.log', access='APPEND')
                  WRITE(10, "('Invalid matrix type. Accepted string are `HS` and `DS`.')")
            CLOSE(unit=10)
            STOP
      END IF


      
      ! -----------------------------------------
      ! computation of normalized spacing
      ! -----------------------------------------

      ! normalized spacing
      norm_spacings(:) = eigenvals(2:size) - eigenvals(1:size-1)
      mean_spac = SUM(norm_spacings) / (size-1)
      norm_spacings(:) = norm_spacings(:) / mean_spac
      CALL checkpoint(debug, str="Mean spacing = ", val=SUM(norm_spacings)/(size-1))

      ! print normalized spacing on file
      WRITE(str, "('results/norm_spac_', I0.1, '_', A, '.dat')") size, TRIM(mattype)
      OPEN(unit=10, file=TRIM(str), access='APPEND')
            WRITE(10, "(ES13.7)") norm_spacings
      CLOSE(unit=10)
      CALL checkpoint(debug, str="Printed normalized spacing on file"//NEW_LINE(""))



      ! -----------------------------------------
      ! final deallocation of matrices/arrays
      ! -----------------------------------------
      
      ! deallocation of matrices/arrays
      DEALLOCATE(AM_H%elem, eigenvals, norm_spacings)
      CALL checkpoint(debug, str="Deallocations perfomed successfully")
      
END PROGRAM eigen_test