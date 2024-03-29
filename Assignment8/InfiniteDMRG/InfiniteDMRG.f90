! Program to simulate the one-dimensional Ising model with Hamiltonian
! H = \lambda \sum_{i=1}^N \sigma_i^z + \sum_{i=1}^{N-1} \sigma_i^x \sigma_{i+1}^x
! with interaction strength \lambda using a coarse-graining procedure to find the ground state energy. 
! The program uses the module "ManyBodyUtils_mod.f90" to perform the many-body functions compactly and 
! the LAPACK routines DSYEVX, DGEMM to perform matrix diagonalization and matrix-matrix multiplication,
! respectively.
! The simulation is performed using the Infinite Density Matrix Renormalization Group procedure.
! Specifically, at each iteration step the system size is increased by two and the exact total Hamiltonian of 
! size ((m_old*d)**2,(m_old*d)**2) is transformed into an effective Hamiltonian of size ((m*d)**2,(m*d)**2) 
! by projecting it in the space of the largest 'm' populations of the reduced density matrices obtained from the 
! total Hamiltonian. At each step, the value of 'm' is updated so that m = min(m_max, 2*m_old).
! The starting system always contains 4 sites, for a total dimension of (2**4,2**4) for the starting total
! Hamiltonian.
! 
! Example of execution:
! ./InfiniteDRMG.out 10000 16 input/lambda.dat input/eps.dat
! will find the energy density of the ground state of the system, until the maximum number of iterations = 10000 
! is reached or the updates of the energy goes below the threshold contained in "input/eps.dat". The value of lambda is
! contained in the file "input/lambda.dat". The evolution of the energy density over the iterations is saved into the file
! "results/spectra/spectrum_10000_16.dat" and the time needed for the computations is appended to "results/timescalings.csv".
!
! Author: Paolo Zinesi
!

! main program
PROGRAM InfiniteDRMG
      USE checkpoint_mod
      USE ManyBodyUtils_mod
      IMPLICIT NONE
      LOGICAL :: debug = .FALSE.

      ! variables to store Hamiltonians and eigenvectors
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: H1, H2, H3, H4
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: H12buff, H34buff
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: H12, H23, H34
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: Htot
      INTEGER :: Nitermax, mmax
      INTEGER :: mm_old, mm
      INTEGER, PARAMETER :: dd = 2
      REAL*8 :: lambda, eps
      REAL*8 :: energydensity = 0.D0, energy_update = 1.D0
      INTEGER :: N_eff = 4

      ! variables to compute reduced density matrices
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: psi_GS
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: rho
      
      ! loop and utility variables
      INTEGER :: ii, iip, jj, jjp, iterRG
      INTEGER :: row_i, col_i
      REAL*8, DIMENSION(2,2) :: sigmaX
      REAL :: ti, tf

      ! variables to read command-line arguments
      CHARACTER(LEN=100) :: arg, lambdafile, epsfile

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

      ! create folder where to store spectra
      CALL SYSTEM('mkdir -p results')
      CALL SYSTEM('mkdir -p results/spectra')


      ! -----------------------------------------
      ! read command-line arguments
      ! -----------------------------------------

      IF (COMMAND_ARGUMENT_COUNT() .EQ. 4) THEN

            ! Nitermax: number of iterations in the RG algorithm
            CALL GET_COMMAND_ARGUMENT(1, arg)
            READ(arg, *) Nitermax

            ! mmax: maximum number of block states to keep at each iteration step
            CALL GET_COMMAND_ARGUMENT(2, arg)
            READ(arg, *) mmax

            ! lambdafile: filename where to read the value of interaction strength lambda
            ! gives an error if input file does not exist
            CALL GET_COMMAND_ARGUMENT(3, lambdafile)
            OPEN(unit=20, file=TRIM(lambdafile), status='OLD')
                  READ(20,*) lambda
            CLOSE(unit=20)

            ! epsfile: filename where to read the value of goal precision eps
            ! gives an error if input file does not exist
            CALL GET_COMMAND_ARGUMENT(4, epsfile)
            OPEN(unit=20, file=TRIM(epsfile), status='OLD')
                  READ(20,*) eps
            CLOSE(unit=20)
            
      ELSE
            ! write error details into an error log file 
            PRINT *, "An error occurred, details on 'errors.log'"
            WRITE(10, "('Exactly 4 command-line arguments are expected:')")
            WRITE(10, "('Nitermax, mmax, lambdafile, epsfile')")
            STOP
      END IF


      ! checks on input validity
      IF ((Nitermax .LE. 0) .OR. (mmax .LT. 2) .OR. (eps .LE. 0.D0)) THEN

            ! write error details into an error log file
            PRINT *, "An error occurred, details on 'errors.log'"
            WRITE(10, "('Variables must be `Nitermax`>0, `mmax`>=2, and `eps`>0')")
            STOP
      END IF

      CALL checkpoint(debug, str="Nitermax = ", val=Nitermax)
      CALL checkpoint(debug, str="mmax = ", val=mmax)
      CALL checkpoint(debug, str="Lambda = ", val=lambda)
      CALL checkpoint(debug, str="Eps = ", val=eps)
      CALL checkpoint(debug, str="")


      ! open file where to store spectrum and print lambda
      WRITE(arg, "('results/spectra/spectrum_', I0, '_', I0, '.dat')") Nitermax, mmax
      OPEN(unit=20, file=TRIM(arg), access='APPEND')
      WRITE(20, "(ES24.17)") lambda


      ! -----------------------------------------
      ! allocations (matrix are allocated with
      ! the maximum required storage)
      ! -----------------------------------------
      ! Hamiltonian allocations
      ALLOCATE(H1(mmax,mmax), H2(dd,dd), H3(dd,dd), H4(mmax,mmax))
      ALLOCATE(H12(mmax*dd,mmax*dd), H23(dd*dd,dd*dd), H34(mmax*dd,mmax*dd))
      ALLOCATE(Htot((mmax*dd)**2,(mmax*dd)**2))
      ALLOCATE(H12buff(mmax*dd,mmax*dd), H34buff(mmax*dd,mmax*dd))

      ! reduced density matrices allocations
      ALLOCATE(psi_GS((mmax*dd)**2,1))
      ALLOCATE(rho(mmax*dd,mmax*dd))

      ! empirical lwork
      lwork = 35*(mmax*dd)**2
      ! allocation of arrays to call DSYEVX
      ALLOCATE(eigenvals((mmax*dd)**2), eigenvects(mmax*dd,mmax))
      ALLOCATE(work(lwork), iwork(5*(mmax*dd)**2), ifail((mmax*dd)**2))

      ! allocation of temporary array to call DGEMM
      ALLOCATE(tmp_matmul(mmax,mmax*dd))
      CALL checkpoint(debug, "Allocations performed successfully")



      ! -----------------------------------------
      ! 0-th iteration
      ! -----------------------------------------
      CALL CPU_TIME(ti)

      ! initialize matrices
      mm = 2
      H1 = 0.D0
      H2 = 0.D0
      H3 = 0.D0
      H4 = 0.D0
      H12 = 0.D0
      H23 = 0.D0
      H34 = 0.D0
      sigmaX = 0.D0
      sigmaX(1,2) = 1.D0
      sigmaX(2,1) = 1.D0

      
      ! fill Hamiltonians of the first iteration (iterRG=0)
      CALL fillSigmaZ(H1(1:mm,1:mm), 1, 1, lambda)
      CALL fillSigmaZ(H2(1:dd,1:dd), 1, 1, lambda)
      CALL fillSigmaZ(H3(1:dd,1:dd), 1, 1, lambda)
      CALL fillSigmaZ(H4(1:mm,1:mm), 1, 1, lambda)
      CALL fillSigmaXSigmaX(H12(1:mm*dd,1:mm*dd), 1, 2)
      CALL fillSigmaXSigmaX(H23(1:dd*dd,1:dd*dd), 1, 2)
      CALL fillSigmaXSigmaX(H34(1:mm*dd,1:mm*dd), 1, 2)

      ! fill interaction buffer Hamiltonians
      H12buff = 0.D0
      H12buff(1:mm*dd,1:mm*dd) = H12buff(1:mm*dd,1:mm*dd) + H12(1:mm*dd,1:mm*dd)
      CALL tensorProductIdentity(H12buff(1:mm*dd,1:mm*dd), H1(1:mm,1:mm), 1)
      CALL tensorProductIdentity(H12buff(1:mm*dd,1:mm*dd), H2(1:dd,1:dd), mm)

      H34buff = 0.D0
      H34buff(1:mm*dd,1:mm*dd) = H34buff(1:mm*dd,1:mm*dd) + H34(1:mm*dd,1:mm*dd)
      CALL tensorProductIdentity(H34buff(1:mm*dd,1:mm*dd), H3(1:dd,1:dd), 1)
      CALL tensorProductIdentity(H34buff(1:mm*dd,1:mm*dd), H4(1:mm,1:mm), dd)
      CALL checkpoint(debug, "Filled single Hamiltonians (first iteration)")

      ! fill total Hamiltonian of the first iteration
      Htot = 0.D0
      CALL tensorProductIdentity(Htot(1:(mm*dd)**2,1:(mm*dd)**2), H1(1:mm,1:mm), 1)
      CALL tensorProductIdentity(Htot(1:(mm*dd)**2,1:(mm*dd)**2), H2(1:dd,1:dd), mm)
      CALL tensorProductIdentity(Htot(1:(mm*dd)**2,1:(mm*dd)**2), H12(1:mm*dd,1:mm*dd), 1)
      CALL tensorProductIdentity(Htot(1:(mm*dd)**2,1:(mm*dd)**2), H23(1:mm*dd,1:mm*dd), mm)
      CALL tensorProductIdentity(Htot(1:(mm*dd)**2,1:(mm*dd)**2), H3(1:dd,1:dd), mm*dd)
      CALL tensorProductIdentity(Htot(1:(mm*dd)**2,1:(mm*dd)**2), H4(1:mm,1:mm), mm*dd*dd)
      CALL tensorProductIdentity(Htot(1:(mm*dd)**2,1:(mm*dd)**2), H34(1:mm*dd,1:mm*dd), mm*dd)
      CALL checkpoint(debug, "Filled total Hamiltonian (first iteration)")
      CALL checkpoint(debug, str="")
      

      
      ! -----------------------------------------
      ! RG iterations
      ! -----------------------------------------
      iterRG = 1
      DO WHILE((iterRG .LE. Nitermax) .AND. (ABS(energy_update) .GT. eps))
            CALL checkpoint(debug, str="Iteration ", val=iterRG)

            ! the updated value of mm cannot be greater than mmax
            mm_old = mm
            mm = MIN(mmax, mm_old*dd)
            CALL checkpoint(debug, str="Prevoius value of m = ", val=mm_old)
            CALL checkpoint(debug, str="New value of m = ", val=mm)


            ! ---------------------------------
            ! diagonalization of Htot
            ! ---------------------------------
            ! call diagonalization routine
            CALL DSYEVX('V', 'I', &                                     ! compute eigenvalues and eigenvectors with index in [IL,IU]
                        'L', &                                          ! store lower triangular matrix
                        (mm_old*dd)**2, &                               ! order of Htot
                        Htot(1:(mm_old*dd)**2,1:(mm_old*dd)**2),&       ! subset of matrix Htot
                        (mm_old*dd)**2, &                               ! leading dimension of Htot
                        0.D0, 1.D0, &                                   ! VL, VU are not referenced since RANGE="I"
                        1, 1, &                                         ! find only first eigenvalue (IL=IU=1)
                        2*DLAMCH('S'), &                                ! absolute error tolerance for the eigenvalues
                        eig_num, eigenvals(1:1),&                       ! eigenvalues and eigenvalues number
                        psi_GS(1:(mm_old*dd)**2,1:1), (mm_old*dd)**2, & ! eigenvectors and their leading dimension
                        work, lwork, iwork, &                           ! work, lwork, iwork variables
                        ifail, info &                                   ! fail indices and info value
            )
            
            ! computation of energy density
            energy_update = energydensity - eigenvals(1)/N_eff
            energydensity = eigenvals(1)/N_eff
            IF(MOD(iterRG, MAX(1,Nitermax/1000)) .EQ. 1) THEN
                  ! write ground state energy for the given iterRG
                  WRITE(20, "(ES24.16, ' ')", advance="no") energydensity
            END IF
            
            !t_H2N_diag = tf-ti
            CALL checkpoint(debug, "INFO value of DSYEVX = ", val=info)
            CALL checkpoint(debug, "Optimal size of LWORK = ", val=work(1))
            CALL checkpoint(debug, "Used size of LWORK = ", val=lwork)
            CALL checkpoint(debug, "Diagonalization of Htot performed successfully")
            

            ! ---------------------------
            ! reduced density matrix
            ! ---------------------------
            ! left reduced density matrix (only upper triangular part)
            rho = 0.D0
            DO iip = 1, mm_old*dd
            DO ii = 1, iip
                  DO jj = 0, dd-1
                  DO jjp = 0, mm_old-1
                        ! compute multi-indices for partial trace
                        row_i = (ii-1) * dd*mm_old + jj * mm_old + jjp + 1
                        col_i = (iip-1) * dd*mm_old + jj * mm_old + jjp + 1

                        ! add value to reduced density matrix
                        rho(ii,iip) = rho(ii,iip) + psi_GS(row_i,1)*psi_GS(col_i,1)
                  END DO
                  END DO
            END DO
            END DO
            CALL checkpoint(debug, "Upper triangular part of rho computed successfully")

            ! ----------------------------------
            ! diagonalization of left subsystem
            ! ----------------------------------
            ! call diagonalization routine
            CALL DSYEVX('V', 'I', &                                     ! compute eigenvalues and eigenvectors with index in [IL,IU]
                        'U', &                                          ! store upper triangular matrix
                        mm_old*dd, &                                    ! order of rho
                        rho(1:mm_old*dd,1:mm_old*dd),&                  ! subset of matrix rho
                        mm_old*dd, &                                    ! leading dimension of rho
                        0.D0, 1.D0, &                                   ! VL, VU are not referenced since RANGE="I"
                        mm_old*dd-mm+1, mm_old*dd, &                    ! find the last mm eigenvalues (from IL=mm_old*dd-mm+1 to UL=mm_old*dd)
                        2*DLAMCH('S'), &                                ! absolute error tolerance for the eigenvalues
                        eig_num, eigenvals(1:mm), &                     ! eigenvalues and eigenvalues number
                        eigenvects(1:mm_old*dd,1:mm), mm_old*dd, &      ! eigenvectors and their leading dimension
                        work, lwork, iwork, &                           ! work, lwork, iwork variables
                        ifail, info &                                   ! fail indices and info value
            )
            CALL checkpoint(debug, "INFO value of DSYEVX = ", val=info)
            CALL checkpoint(debug, "Optimal size of LWORK = ", val=work(1))
            CALL checkpoint(debug, "Used size of LWORK = ", val=lwork)
            CALL checkpoint(debug, "Largest population of  rho = ", val=eigenvals(mm))
            CALL checkpoint(debug, "Diagonalization of rho performed successfully")


            !-----------------------------------------------------
            ! projection of Hamiltonians in the restricted 
            ! basis for left subsystem
            !-----------------------------------------------------
            ! call matrix multiplication routines to find H1 of the next iteration
            CALL DGEMM('T', 'N', &                                      ! transpose A, B not transposed 
                       mm, mm_old*dd, mm_old*dd, 1.D0, &                ! M=rows of A**T, N=columns of B, K=cols of A**T, alpha
                       eigenvects(1:mm_old*dd,1:mm), mm_old*dd, &       ! A matrix (eigenvectors) with its leading dimension
                       H12buff(1:mm_old*dd,1:mm_old*dd), mm_old*dd, &   ! B matrix (H12buff) with its leading dimension
                       0.D0, tmp_matmul(1:mm,1:mm_old*dd), mm &         ! buffer matrix and its leading dimension
            )
            CALL DGEMM('N', 'N', &                                      ! A and B not transposed 
                       mm, mm, mm_old*dd, 1.D0, &                       ! M=rows of A, N=columns of B, K=cols of A, alpha
                       tmp_matmul(1:mm,1:mm_old*dd), mm, &              ! A matrix (tmp_matmul) with its leading dimension
                       eigenvects(1:mm_old*dd,1:mm), mm_old*dd, &       ! B matrix (eigenvectors) with its leading dimension
                       0.D0, H1(1:mm,1:mm), mm &                        ! NEW H1 matrix and its leading dimension
            )


            H12 = 0.D0
            CALL tensorProductIdentity(H12(1:mm_old*dd,1:mm_old*dd), sigmaX, mm_old)
            ! call matrix multiplication routines to find H12 of the next iteration
            CALL DGEMM('T', 'N', &                                      ! transpose A, B not transposed 
                       mm, mm_old*dd, mm_old*dd, 1.D0, &                ! M=rows of A**T, N=columns of B, K=cols of A**T, alpha
                       eigenvects(1:mm_old*dd,1:mm), mm_old*dd, &       ! A matrix (eigenvectors) with its leading dimension
                       H12(1:mm_old*dd,1:mm_old*dd), mm_old*dd, &       ! B matrix (H12) with its leading dimension
                       0.D0, tmp_matmul(1:mm,1:mm_old*dd), mm &         ! buffer matrix and its leading dimension
            )
            CALL DGEMM('N', 'N', &                                      ! A and B not transposed 
                       mm, mm, mm_old*dd, 1.D0, &                       ! M=rows of A, N=columns of B, K=cols of A, alpha
                       tmp_matmul(1:mm,1:mm_old*dd), mm, &              ! A matrix (tmp_matmul) with its leading dimension
                       eigenvects(1:mm_old*dd,1:mm), mm_old*dd, &       ! B matrix (eigenvectors) with its leading dimension
                       0.D0, H12buff(1:mm,1:mm), mm &                   ! buffer matrix and its leading dimension
            )

            ! filling of NEW H12 matrix using sigmaX and the buffer H12buff
            H12 = 0.D0
            CALL generalTensorProduct(H12(1:mm*dd,1:mm*dd), H12buff(1:mm,1:mm), sigmaX)
            CALL checkpoint(debug, str="Successful calculation of next step matrices H1,H12")

            ! filling symmetric Hamiltonians
            H4 = H1
            H34 = 0.D0
            CALL generalTensorProduct(H34(1:mm*dd,1:mm*dd), sigmaX, H12buff(1:mm,1:mm))

            ! ----------------------------------------
            ! creation of Hamiltonians for next step
            ! ----------------------------------------
            ! up to know H1, H12, H4, H34 have been updated, while H2, H3, H23 are constant
            ! we now complete the algorithm by filling the block matrices H12buff, H23, Htot

            ! fill interaction buffer Hamiltonians
            H12buff = 0.D0
            H12buff(1:mm*dd,1:mm*dd) = H12buff(1:mm*dd,1:mm*dd) + H12(1:mm*dd,1:mm*dd)
            CALL tensorProductIdentity(H12buff(1:mm*dd,1:mm*dd), H1(1:mm,1:mm), 1)
            CALL tensorProductIdentity(H12buff(1:mm*dd,1:mm*dd), H2(1:dd,1:dd), mm)

            H34buff = 0.D0
            H34buff(1:mm*dd,1:mm*dd) = H34buff(1:mm*dd,1:mm*dd) + H34(1:mm*dd,1:mm*dd)
            CALL tensorProductIdentity(H34buff(1:mm*dd,1:mm*dd), H3(1:dd,1:dd), 1)
            CALL tensorProductIdentity(H34buff(1:mm*dd,1:mm*dd), H4(1:mm,1:mm), dd)
            CALL checkpoint(debug, "Filled single Hamiltonians, iteration = ", val=iterRG)

            ! fill total Hamiltonian
            Htot = 0.D0
            CALL tensorProductIdentity(Htot(1:(mm*dd)**2,1:(mm*dd)**2), H1(1:mm,1:mm), 1)
            CALL tensorProductIdentity(Htot(1:(mm*dd)**2,1:(mm*dd)**2), H2(1:dd,1:dd), mm)
            CALL tensorProductIdentity(Htot(1:(mm*dd)**2,1:(mm*dd)**2), H12(1:mm*dd,1:mm*dd), 1)
            CALL tensorProductIdentity(Htot(1:(mm*dd)**2,1:(mm*dd)**2), H23(1:dd*dd,1:dd*dd), mm)
            CALL tensorProductIdentity(Htot(1:(mm*dd)**2,1:(mm*dd)**2), H3(1:dd,1:dd), mm*dd)
            CALL tensorProductIdentity(Htot(1:(mm*dd)**2,1:(mm*dd)**2), H4(1:mm,1:mm), mm*dd*dd)
            CALL tensorProductIdentity(Htot(1:(mm*dd)**2,1:(mm*dd)**2), H34(1:mm*dd,1:mm*dd), mm*dd)
            CALL checkpoint(debug, "Filled total Hamiltonian, iteration = ", val=iterRG)
            CALL checkpoint(debug, str="")

            ! effective system size described by the current Hamiltonian
            N_eff = 4 + 2*iterRG
            iterRG = iterRG + 1

      END DO

      ! write ground state energy of the last given iterRG
      WRITE(20, "(ES24.16, ' ')", advance="no") energydensity
      WRITE(20, *)
      WRITE(20, *)
      CALL CPU_TIME(tf)

      ! -----------------------------------------
      ! write results
      ! -----------------------------------------
      ! write timescaling
      OPEN(unit=21, file="results/timescalings.csv", access='APPEND')
      WRITE(21, "(I0, ',', I0, ',', ES24.17, ',' I0, ',', ES13.7)") &
            mmax, Nitermax, lambda, iterRG, tf-ti
      CLOSE(unit=21)


      ! -----------------------------------------
      ! final deallocations
      ! -----------------------------------------
      ! Hamiltonians
      DEALLOCATE(H1, H2, H3, H4)
      DEALLOCATE(H12, H23, H34)
      DEALLOCATE(Htot)
      DEALLOCATE(H12buff, H34buff)

      ! reduced density matrices
      DEALLOCATE(psi_GS)
      DEALLOCATE(rho)

      ! arrays to call DSYEVX
      DEALLOCATE(eigenvals, eigenvects)
      DEALLOCATE(work, iwork, ifail)

      ! temporary array to call DGEMM
      DEALLOCATE(tmp_matmul)

      CALL checkpoint(debug, "Deallocations performed successfully")

      ! close 'errors.log'
      CLOSE(unit=10)

      ! close spectrum file
      CLOSE(unit=20)


END PROGRAM InfiniteDRMG








