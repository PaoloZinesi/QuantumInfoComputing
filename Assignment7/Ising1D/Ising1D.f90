! Program to simulate the one-dimensional Ising model with Hamiltonian
! H = \lambda \sum_{i=1}^N \sigma_i^z + \sum_{i=1}^{N-1} \sigma_i^x \sigma_{i+1}^x
! with interaction strength \lambda and N subsystem with local dimension d=2.
! The computational cost increases exponentially with the number of subsystems N given the
! fact that an Hamiltonian of N subsystems is (2**N,2**N) dimensional.
! The program uses a dedicated subroutine to fill the antidiagonal elements to fill the sigma_x sigma_x
! interaction terms, but the sigma_z terms are computed with a single DO loop without subroutines.
! Results and timescaling are saved into external files to be processed later.
! 
! Example of execution:
! ./Ising1D 6 10 input/lambda.dat
! will create a (2**6,2**6) = (64,64) Hamiltonian and compute the 10 lowest eigenvalues using the
! value of lambda stored into 'input/lambda.dat'. Results are stored in the folder "results".
! For machines with 8 GB of RAM, the maximum allowed N is 14.
!
! Author: Paolo Zinesi
!

! main program
PROGRAM Ising1D
      USE checkpoint_mod
      IMPLICIT NONE
      LOGICAL :: debug = .FALSE.

      ! variables to store Hamiltonian and to compute eigenvalues
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: HN
      INTEGER :: NN, K_max
      REAL*8 :: lambda

      ! loop and utility variables
      INTEGER :: ii, jj, bb
      REAL :: ti, tf

      ! variables to read command-line arguments
      CHARACTER(LEN=100) :: arg, lambdafile

      ! variables to call DSYEVX subroutine
      REAL*8, EXTERNAL :: DLAMCH
      REAL*8, DIMENSION(:), ALLOCATABLE :: eigenvals, work
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: eigenvects
      INTEGER :: eig_num, lwork, info
      INTEGER, DIMENSION(:), ALLOCATABLE :: iwork, ifail

      ! variables to store timescalings
      REAL :: t_creation, t_diag



      INTERFACE
            ! subroutine to fill matrix with sigma_x sigma_x interaction
            SUBROUTINE fillAntidiag(AA, k)
                  REAL*8, DIMENSION(:,:), INTENT(INOUT) :: AA
                  INTEGER, INTENT(IN) :: k
            END SUBROUTINE fillAntidiag

            ! subroutine to print matrix on file
            SUBROUTINE writeMatFile(M, unit, file, format)
                  REAL*8, DIMENSION(:,:), INTENT(IN) :: M
                  INTEGER :: unit
                  CHARACTER (LEN=*) :: file, format
            END SUBROUTINE writeMatFile

      END INTERFACE


      ! open 'errors.log' and 'timescalings.log
      OPEN(unit=10, file='errors.log', access='APPEND')

      ! create directory where to store results
      CALL SYSTEM('mkdir -p results')
      CALL SYSTEM('mkdir -p results/spectra')



      ! -----------------------------------------
      ! read command-line arguments
      ! -----------------------------------------

      IF (COMMAND_ARGUMENT_COUNT() .EQ. 3) THEN

            ! NN: number of spins in the lattice
            CALL GET_COMMAND_ARGUMENT(1, arg)
            READ(arg, *) NN

            ! K_max: number of eigenvalues to compute
            CALL GET_COMMAND_ARGUMENT(2, arg)
            READ(arg, *) K_max

            ! lambdafile: filename where to read the value of interaction strength lambda
            ! gives an error if input file does not exist
            CALL GET_COMMAND_ARGUMENT(3, lambdafile)
            OPEN(unit=20, file=TRIM(lambdafile), status='OLD')
                  READ(20,*) lambda
            CLOSE(unit=20)
            
      ELSE
            ! write error details into an error log file 
            PRINT *, "An error occurred, details on 'errors.log'"
            WRITE(10, "('Exactly 3 command-line arguments are expected:')")
            WRITE(10, "('NN, K_max, lambdafile')")
            STOP
      END IF


      ! checks on input validity
      IF ((NN .LT. 2) .OR. (K_max .LE. 0) .OR. (K_max .GT. 2**NN)) THEN

            ! write error details into an error log file
            PRINT *, "An error occurred, details on 'errors.log'"
            WRITE(10, "('Variables must be `NN`>=2 and 0<`K_max`<=2**NN')")
            STOP
      END IF

      CALL checkpoint(debug, str="NN = ", val=NN)
      CALL checkpoint(debug, str="K_max = ", val=K_max)
      CALL checkpoint(debug, str="Lambda = ", val=lambda)
      CALL checkpoint(debug, str="")



      ! -----------------------------------------
      ! filling of matrix
      ! -----------------------------------------
      ALLOCATE(HN(2**NN,2**NN))
      CALL checkpoint(debug, "Allocation performed successfully")
      CALL CPU_TIME(ti)

      ! initialize all elements to zero
      HN = 0.D0

      ! fill sigma_z contributions
      DO jj = 1, 2**NN

            ! exploit the alternating signs for the sigma_z Hamiltonians
            ! to easily compute the diagonal elements
            HN(jj,jj) = lambda * SUM((/ ((-1)**((jj-1)/(2**(NN-ii))), ii=1,NN) /))
      END DO
      CALL checkpoint(debug, "HN(1,1) (should be equal to lambda*N) = ", val=HN(1,1))

      
      ! fill sigma_i^x sigma_{i+1}^x interactions
      DO ii = 1, NN-1

            ! repeat the same procedure for all the 2**(ii-1) diagonal blocks
            DO bb = 1, 2**(ii-1)
                  CALL fillAntidiag(HN((bb-1)*2**(NN+1-ii)+1:bb*2**(NN+1-ii), &
                                       (bb-1)*2**(NN+1-ii)+1:bb*2**(NN+1-ii)), &
                                    NN+1-ii)
            END DO
      END DO
      CALL CPU_TIME(tf)
      t_creation = tf - ti
      CALL checkpoint(debug, "HN(1,3*2**(N-2)+1) (should be 1) = ", val=HN(1,3*2**(NN-2)+1))
      CALL checkpoint(debug, "Matrix has been filled completely")


      ! printing matrix on file in debugging mode
      IF(debug .EQV. .TRUE.) THEN
            CALL writeMatFile(HN, unit=30, file="results/HN.dat", format="(ES13.5)")
      END IF



      ! -----------------------------------------
      ! diagonalization of matrix
      ! -----------------------------------------

      ! empirical lwork
      lwork = 35*(2**NN)

      ! allocation of arrays to call DSYEVX
      ALLOCATE(eigenvals(2**NN), eigenvects(1,1))
      ALLOCATE(work(lwork), iwork(5*(2**NN)), ifail(1))
      CALL CPU_TIME(ti)

      ! call diagonalization routine
      CALL DSYEVX('N', 'I', &             ! compute only eigenvalues with index in [IL,IU]
                  'L', 2**NN, HN, 2**NN,& ! store lower triangular matrix of 2**NN x 2**NN matrix HN
                  0.D0, 1.D0, &           ! VL, VU are not referenced since RANGE="I"
                  1, K_max, &             ! find eigenvalues from IL=1 to IU=K_max
                  2*DLAMCH('S'), &        ! absolute error tolerance for the eigenvalues
                  eig_num, eigenvals, &   ! eigenvalues and eigenvalues number
                  eigenvects, 1, &        ! eigenvectors and leading dimension, not referenced
                  work, lwork, iwork, &   ! work, lwork, iwork variables
                  ifail, info &           ! fail indices and info value
      )

      CALL CPU_TIME(tf)
      t_diag = tf - ti
      CALL checkpoint(debug, "INFO value of DSYEVX = ", val=info)
      CALL checkpoint(debug, "Optimal size of LWORK = ", val=work(1))
      CALL checkpoint(debug, "Used size of LWORK = ", val=lwork)



      ! -----------------------------------------
      ! write results
      ! -----------------------------------------

      ! write timescaling
      OPEN(unit=20, file="results/timescalings.csv", access='APPEND')
            WRITE(20, "(I0, ',', I0, ',', ES24.17, ',', ES13.7, ',', ES13.7)") &
                  NN, K_max, lambda, t_creation, t_diag
      CLOSE(unit=20)


      ! write spectrum
      WRITE(arg, "('results/spectra/spectrum_', I0, '_', I0, '.dat')") NN, K_max
      OPEN(unit=21, file=TRIM(arg), access='APPEND')

            ! print lambda and the spectrum
            WRITE(21, "(ES24.17)") lambda

            DO ii = 1, K_max
                  WRITE(21, "(ES24.17, ' ')", advance="no") eigenvals(ii)
            END DO
            WRITE(21, *)
            WRITE(21, *)

      CLOSE(unit=21)



      ! -----------------------------------------
      ! final deallocations
      ! -----------------------------------------

      DEALLOCATE(HN)
      DEALLOCATE(eigenvals, eigenvects)
      DEALLOCATE(work, iwork, ifail)
      CALL checkpoint(debug, str="Deallocations perfomed successfully")

      ! close 'errors.log'
      CLOSE(unit=10)


END PROGRAM Ising1D





! This subroutine fills an input matrix of dimension (2**k,2**k) with sub-matrices of the form
! I=id_{2**(k-2)} on the 4 antidiagonal blocks. The typical form of the filled matrix AA is
!     AA =  0 0 0 I 
!           0 0 I 0
!           0 I 0 0
!           I 0 0 0
!
! inputs:
! - AA [REAL*8, dimension(:,:)]: input REAL*8 (sub)matrix to fill as explained above
! - k [integer]: base-2 logarithm of matrix size, such that 2**k = SIZE(AA,1) = SIZE(AA,2)
!                              
! outputs: None
! 
! TODO: None
!
SUBROUTINE fillAntidiag(AA, k)
      ! matrix to fill
      REAL*8, DIMENSION(:,:), INTENT(INOUT) :: AA

      ! 2**k = size(AA,1)
      INTEGER, INTENT(IN) :: k

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


      ! -------------------------
      ! filling of submatrices
      ! -------------------------

      ! leftermost identity sub-matrix
      DO ii = 1, 2**(k-2)
            AA(3*(2**(k-2))+ii, ii) = AA(3*(2**(k-2))+ii, ii) + 1
      END DO

      ! second leftermost identity sub-matrix
      DO ii = 1, 2**(k-2)
            AA(2*(2**(k-2))+ii, 1*(2**(k-2))+ii) = AA(2*(2**(k-2))+ii, 1*(2**(k-2))+ii) + 1
      END DO

      ! second rightermost identity sub-matrix
      DO ii = 1, 2**(k-2)
            AA(1*(2**(k-2))+ii, 2*(2**(k-2))+ii) = AA(1*(2**(k-2))+ii, 2*(2**(k-2))+ii) + 1
      END DO

      ! rightermost identity sub-matrix
      DO ii = 1, 2**(k-2)
            AA(ii, 3*(2**(k-2))+ii) = AA(ii, 3*(2**(k-2))+ii) + 1
      END DO

END SUBROUTINE fillAntidiag


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

END SUBROUTINE writeMatFile