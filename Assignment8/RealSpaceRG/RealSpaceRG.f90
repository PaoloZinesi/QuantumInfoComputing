! Program to simulate the one-dimensional Ising model with Hamiltonian
! H = \lambda \sum_{i=1}^N \sigma_i^z + \sum_{i=1}^{N-1} \sigma_i^x \sigma_{i+1}^x
! with interaction strength \lambda using a coarse-graining procedure to find the ground state energy.
! The program uses the module...
! 
! Example of execution:
! ./RealSpaceRG 5 10 input/lambda.dat
! will ...
!
! Author: Paolo Zinesi
!

! main program
PROGRAM RealSpaceRG
      USE checkpoint_mod
      USE ManyBodyUtils_mod
      IMPLICIT NONE
      LOGICAL :: debug = .FALSE.

      ! variables to store Hamiltonian and to compute eigenvalues
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: HN, mat2N
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: AN, BN
      INTEGER :: NN, Nitermax
      REAL*8 :: lambda

      ! loop and utility variables
      INTEGER :: ii, iterRG
      REAL*8, DIMENSION(2,2) :: sigmaX
      REAL :: ti, tf

      ! variables to store timescalings
      REAL :: t_H2N_creation, t_H2N_diag, t_N_matmul

      ! variables to read command-line arguments
      CHARACTER(LEN=100) :: arg, lambdafile

      ! variables to call DSYEVX subroutine
      REAL*8, EXTERNAL :: DLAMCH
      REAL*8, DIMENSION(:), ALLOCATABLE :: eigenvals, work
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: eigenvects
      INTEGER :: eig_num, lwork, info
      INTEGER, DIMENSION(:), ALLOCATABLE :: iwork, ifail

      ! variables to call DGEMM subroutine
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: tmp_matmul

      
      ! open 'errors.log'
      OPEN(unit=10, file='errors.log', access='APPEND')

      ! create file where to store spectra
      CALL SYSTEM('mkdir -p results')
      CALL SYSTEM('mkdir -p results/spectra')


      ! -----------------------------------------
      ! read command-line arguments
      ! -----------------------------------------

      IF (COMMAND_ARGUMENT_COUNT() .EQ. 3) THEN

            ! NN: number of spins in the lattice before doubling the size 
            CALL GET_COMMAND_ARGUMENT(1, arg)
            READ(arg, *) NN

            ! Nitermax: number of iterations in the RG algorithm
            CALL GET_COMMAND_ARGUMENT(2, arg)
            READ(arg, *) Nitermax

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
            WRITE(10, "('NN, Nitermax, lambdafile')")
            STOP
      END IF


      ! checks on input validity
      IF ((NN .LT. 2) .OR. (Nitermax .LE. 0)) THEN

            ! write error details into an error log file
            PRINT *, "An error occurred, details on 'errors.log'"
            WRITE(10, "('Variables must be `NN`>=2 and `Nitermax`>0')")
            STOP
      END IF

      CALL checkpoint(debug, str="NN = ", val=NN)
      CALL checkpoint(debug, str="Nitermax = ", val=Nitermax)
      CALL checkpoint(debug, str="Lambda = ", val=lambda)
      CALL checkpoint(debug, str="")


      ! open file where to store spectrum and print lambda
      WRITE(arg, "('results/spectra/spectrum_', I0, '_', I0, '.dat')") NN, Nitermax
      OPEN(unit=20, file=TRIM(arg), access='APPEND')
      WRITE(20, "(ES24.17)") lambda


      ! -----------------------------------------
      ! allocations
      ! -----------------------------------------
      ! matrix allocations
      ALLOCATE(HN(2**NN,2**NN))
      ALLOCATE(mat2N(2**(2*NN),2**(2*NN)))
      ALLOCATE(AN(2**NN,2**NN),BN(2**NN,2**NN))

      ! empirical lwork
      lwork = 35*(2**(2*NN))
      ! allocation of arrays to call DSYEVX
      ALLOCATE(eigenvals(2**(2*NN)), eigenvects(2**(2*NN),2**NN))
      ALLOCATE(work(lwork), iwork(5*(2**(2*NN))), ifail(2**(2*NN)))

      ! allocation of temporary array to call DGEMM
      ALLOCATE(tmp_matmul(2**NN,2**(2*NN)))
      CALL checkpoint(debug, "Allocations performed successfully")



      ! -----------------------------------------
      ! 0-th iteration
      ! -----------------------------------------

      ! initialize elements
      HN = 0.D0
      AN = 0.D0
      BN = 0.D0
      sigmaX = 0.D0
      sigmaX(1,2) = 1.D0
      sigmaX(2,1) = 1.D0


      ! fill HN of the first iteration (iterRG=0)
      DO ii = 1, NN
            ! sigma_z contribution
            CALL fillSigmaZ(HN, ii, NN, lambda=lambda)
      END DO
      DO ii = 1, NN-1
            ! sigma_x sigma_x contribution
            CALL fillSigmaXSigmaX(HN, ii, NN, lambda=1.D0)
      END DO
      CALL checkpoint(debug, "Filled starting HN")


      ! fill AN of the first iteration (iterRG=0)
      CALL tensorProductIdentity(AN, sigmaX, 2**(NN-1))
      CALL checkpoint(debug, "Filled starting AN")


      ! fill BN of the first iteration (iterRG=0)
      CALL tensorProductIdentity(BN, sigmaX, 2**0)
      CALL checkpoint(debug, "Filled starting BN")
      CALL checkpoint(debug, str="")
      

      
      ! -----------------------------------------
      ! RG iterations
      ! -----------------------------------------

      DO iterRG = 1, Nitermax
            CALL checkpoint(debug, str="Iteration ", val=iterRG)
            
            ! --------------------------------- 
            ! creation of H2N matrix
            ! ---------------------------------
            mat2N = 0.D0
            CALL CPU_TIME(ti)

            ! fill HN * id_NN
            CALL tensorProductIdentity(mat2N, HN, 2**0)

            ! fill id_NN * HN
            CALL tensorProductIdentity(mat2N, HN, 2**NN)

            ! fill AN * BN
            CALL generalTensorProduct(mat2N, AN, BN)

            CALL CPU_TIME(tf)
            t_H2N_creation = tf-ti



            ! ---------------------------------
            ! diagonalization of H2N matrix
            ! ---------------------------------
            CALL CPU_TIME(ti)
            ! call diagonalization routine
            CALL DSYEVX('V', 'I', &                         ! compute eigenvalues and eigenvectors with index in [IL,IU]
                        'L', 2**(2*NN), mat2N, 2**(2*NN), & ! store lower triangular matrix of 2**NN x 2**NN matrix HN
                        0.D0, 1.D0, &                       ! VL, VU are not referenced since RANGE="I"
                        1, 2**NN, &                         ! find eigenvalues from IL=1 to IU=2**NN
                        2*DLAMCH('S'), &                    ! absolute error tolerance for the eigenvalues
                        eig_num, eigenvals, &               ! eigenvalues and eigenvalues number
                        eigenvects, 2**(2*NN), &            ! eigenvectors and their leading dimension
                        work, lwork, iwork, &               ! work, lwork, iwork variables
                        ifail, info &                       ! fail indices and info value
            )
            CALL CPU_TIME(tf)
            t_H2N_diag = tf-ti
            CALL checkpoint(debug, "INFO value of DSYEVX = ", val=info)
            CALL checkpoint(debug, "Optimal size of LWORK = ", val=work(1))
            CALL checkpoint(debug, "Used size of LWORK = ", val=lwork)
            

            ! ---------------------------
            ! matrices for next iteration
            ! ---------------------------
            ! fill HN (diagonal) matrix of the next step
            HN = 0.D0
            DO ii = 1, 2**NN
                  HN(ii,ii) = 0.5D0 * eigenvals(ii)
            END DO

            ! fill AN matrix of the next step
            mat2N = 0.D0
            CALL tensorProductIdentity(mat2N, AN, 2**NN)

            CALL CPU_TIME(ti)
            ! call matrix multiplication routines to compute AN matrix of the next step
            CALL DGEMM('T', 'N', &                          ! transpose A, B not transposed 
                       2**NN, 2**(2*NN), 2**(2*NN), 1.D0, & ! M=rows of A**T, N=columns of B, K=cols of A**T, alpha
                       eigenvects, 2**(2*NN), &             ! A matrix (eigenvectors) with its leading dimension
                       mat2N, 2**(2*NN), &                  ! B matrix (H2N) with its leading dimension
                       0.D0, tmp_matmul, 2**NN &            ! buffer matrix and its leading dimension
            )
            CALL DGEMM('N', 'N', &                                ! A and B not transposed 
                       2**NN, 2**NN, 2**(2*NN), SQRT(0.5D0), &    ! M=rows of A, N=columns of B, K=cols of A, alpha
                       tmp_matmul, 2**NN, &                       ! A matrix (tmp_matmul) with its leading dimension
                       eigenvects, 2**(2*NN), &                   ! B matrix (eigenvectors) with its leading dimension
                       0.D0, AN, 2**NN &                          ! new AN matrix and its leading dimension
            )
            CALL CPU_TIME(tf)
            t_N_matmul = 0.5*(tf-ti)

            ! fill BN matrix of the next step
            mat2N = 0.D0
            CALL tensorProductIdentity(mat2N, BN, 2**0)

            CALL CPU_TIME(ti)
            ! call matrix multiplication routines to compute BN matrix of the next step
            CALL DGEMM('T', 'N', &                          ! transpose A, B not transposed 
                       2**NN, 2**(2*NN), 2**(2*NN), 1.D0, & ! M=rows of A**T, N=columns of B, K=cols of A**T, alpha
                       eigenvects, 2**(2*NN), &             ! A matrix (eigenvectors) with its leading dimension
                       mat2N, 2**(2*NN), &                  ! B matrix (H2N) with its leading dimension
                       0.D0, tmp_matmul, 2**NN &            ! buffer matrix and its leading dimension
            )
            CALL DGEMM('N', 'N', &                                ! A and B not transposed 
                       2**NN, 2**NN, 2**(2*NN), SQRT(0.5D0), &    ! M=rows of A, N=columns of B, K=cols of A, alpha
                       tmp_matmul, 2**NN, &                       ! A matrix (tmp_matmul) with its leading dimension
                       eigenvects, 2**(2*NN), &                   ! B matrix (eigenvectors) with its leading dimension
                       0.D0, BN, 2**NN &                          ! new BN matrix and its leading dimension
            )
            CALL CPU_TIME(tf)
            t_N_matmul = t_N_matmul + 0.5*(tf-ti)
            CALL checkpoint(debug, str="Time needed [s] =  ", val=tf-ti)
            CALL checkpoint(debug, str="E0/N =  ", val=HN(1,1)/NN)
            CALL checkpoint(debug, str="")


            ! -----------------------------------------
            ! write results
            ! -----------------------------------------
            ! write timescaling
            OPEN(unit=21, file="results/timescalings.csv", access='APPEND')
                  !WRITE(21, "(I0, ',', I0, ',', ES24.17, ',', ES13.7, ',', ES13.7, ',', ES13.7)") &
                  WRITE(21, "(I0, ',', I0, ',', ES24.17, 3(',', ES13.7))") &
                        NN, iterRG, lambda, t_H2N_creation, t_H2N_diag, t_N_matmul
            CLOSE(unit=21)

            ! write ground state energy for the given iterRG
            WRITE(20, "(ES24.17, ' ')", advance="no") HN(1,1)/NN

      END DO
      WRITE(20, *)
      WRITE(20, *)


      ! -----------------------------------------
      ! final deallocations
      ! -----------------------------------------

      DEALLOCATE(HN, mat2N, AN, BN)
      DEALLOCATE(eigenvals, eigenvects)
      DEALLOCATE(work, iwork, ifail)
      DEALLOCATE(tmp_matmul)
      CALL checkpoint(debug, str="Deallocations perfomed successfully")

      ! close 'errors.log'
      CLOSE(unit=10)

      ! close spectrum file
      CLOSE(unit=20)


END PROGRAM RealSpaceRG








