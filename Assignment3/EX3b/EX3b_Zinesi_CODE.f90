! Program to test random matrix and eigenvalue subroutines on LAPACK.
! This program uses the DCmatrix 'RandMat' subroutine to generate random Hermitian matrices
! with complex entries uniformly distributed in {-1 < Re(z),Im(z) < 1} in the complex plane
! and 'eigenvaluesMat' subroutine to compute the eigenvalues of such matrix.
! The normalized spacings between eigenvalues are stored in an output file.
! 
! The (linear) dimension of the random square matrix to be generated is given in the command-line.
! It must be greater than 1, otherwise the eigenvalue spacings are not defined.
! 
! Example of execution:
! ./a.out 1000
! will create a 1000x1000 random Hermitian matrix, compute its eigenvalues and store the 
! normalized spacings in 'results/norm_spac.dat'.
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

      ! define matrices to be allocated later
      TYPE(DCmatrix) :: AM
      INTEGER :: AMsize

      ! eigenvalues arrays and normalized spacings
      REAL*8, DIMENSION(:), ALLOCATABLE :: eigenvals, norm_spacings
      REAL*8 :: mean_spac


      ! initialize seed for reproducibility
      INTEGER, DIMENSION(4) :: iseed = (/ 0, 0, 0, 0 /) 


      ! -----------------------------------------
      ! read command-line arguments
      ! -----------------------------------------

      IF (COMMAND_ARGUMENT_COUNT() .EQ. 1) THEN
            arg_idx = 1

            ! get argument number "arg_idx" and put it into "AMsize"
            CALL GET_COMMAND_ARGUMENT(arg_idx, arg)
            READ(arg, *) AMsize
            
      ELSE
            PRINT *, "Exactly 1 command-line argument is expected"
            STOP
      END IF


      ! check if at least N > 1
      IF(AMsize .LE. 1) THEN
            PRINT *, "N > 1 is requires to compute normalized spacing"
            STOP
      END IF
      CALL checkpoint(debug, str="Linear dimension of A = ", val=AMsize)
      CALL checkpoint(debug, str="")
      

      ! -----------------------------------------
      ! creation and filling of matrix AM
      ! -----------------------------------------
      
      ! allocation and creation of random DCmatrix using the dedicated subroutine
      AM = DCmatrix(rows=AMsize, cols=AMsize)
      CALL RandMat(M_rand=AM, dist='S', iseed=iseed, sym='H')


      ! create directory where to store results
      CALL SYSTEM('mkdir -p results')

      ! write the random matrix on file for debugging purposes
      IF (debug .EQV. .TRUE.) THEN
            CALL writeMatFile(AM, unit=10, file='results/AM.dat', format="(ES12.5,SP,ES13.5,' i   ')")
            PRINT '(A)', "Matrix printed on 'results/AM.dat'"
      END IF



      ! -----------------------------------------
      ! eigendecomposition of matrix AM
      ! -----------------------------------------

      ! allocation of eigenvalues array
      ALLOCATE(eigenvals(AMsize))
      
      ! computation of eigenvalues with the dedicated subroutine
      CALL eigenvaluesMat(AM, eigenvals)


      
      ! -----------------------------------------
      ! computation of normalized spacing
      ! -----------------------------------------

      ! allocation of normalized eigenvalue spacing
      ALLOCATE(norm_spacings(AMsize-1))
      
      ! normalized spacing
      norm_spacings(:) = eigenvals(2:AMsize) - eigenvals(1:AMsize-1)
      mean_spac = SUM(norm_spacings) / (AMsize-1)
      norm_spacings(:) = norm_spacings(:) / mean_spac
      CALL checkpoint(debug, str="Mean spacing = ", val=SUM(norm_spacings)/(AMsize-1))

      ! print normalized spacing on file
      OPEN(unit=10, file="results/norm_spac.dat", access='APPEND')
            WRITE(10, "(ES13.7)") norm_spacings
      CLOSE(unit=10)
      PRINT '(A)', "Normalized spacings printed on 'results/norm_spac.dat'"



      ! -----------------------------------------
      ! final deallocation of matrices/arrays
      ! -----------------------------------------
      DEALLOCATE(AM%elem, eigenvals, norm_spacings)
      CALL checkpoint(debug, str="Deallocation of matrices perfomed successfully")
      
END PROGRAM eigen_test