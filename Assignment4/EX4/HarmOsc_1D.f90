! Program to generate selected eigenvalues and eigenfunctions of the 1-dimensional harmonic oscillator.
! This program uses the subroutines DSTEVX to compute selected eigenvalues and eigenvectors of a
! discretized real, symmetric, tridiagonal matrix (i.e., the Hamiltonian). A further improvement would
! be the restriction of the system for x>=0, to remove half of the computational cost. However, I found
! no system of linear equations that can be casted into a symmetric tridiagonal form, which is particularly
! easy to compute compared to the full nonsymmetric problem.
! 
! The number of function evaluation points (N), the number of eigenvalues (k), the file containing the
! angular frequency (omega), the file where to store results are given as command-line arguments.
! 
! Example of execution:
! ./a.out 10000 10 input/omega.dat results/out.dat
! will evaluate N=10000 function points of the 1D harmonic oscillator, find k=10 eigenvalues,
! take the value of omega from 'input/omega.dat' file and write the results on 'results/out.dat'.
! The maximum length L of the system is computed as 4*sqrt(k_max/omega) and it adapts automatically
! on the system typical size to increase flexibility.
!
! Author: Paolo Zinesi
!

! main program
PROGRAM HarmOsc_1D
      USE checkpoint_mod
      IMPLICIT NONE
      LOGICAL :: debug = .FALSE.

      ! variables to read command-line arguments
      CHARACTER(LEN=100) :: arg, omega_file, outfile

      ! input information
      INTEGER :: Ntot, k_max
      REAL*8 :: omega

      ! space variables
      REAL*8 :: Ltot, deltax
      REAL*8, DIMENSION(:), ALLOCATABLE :: xgrid

      ! define matrix elements
      REAL*8, DIMENSION(:), ALLOCATABLE :: diag, upper_diag

      ! variables to call DSTEVX subroutine
      REAL*8, EXTERNAL :: DLAMCH
      REAL*8, DIMENSION(:), ALLOCATABLE :: eigenvals, work
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: eigenvects
      INTEGER :: eig_num, info
      INTEGER, DIMENSION(:), ALLOCATABLE :: iwork, ifail

      ! normalization
      REAL*8 :: norm2

      ! utility and loop variables
      INTEGER :: ii, jj

      ! open 'errors.log'
      OPEN(unit=10, file='errors.log', access='APPEND')


      ! -----------------------------------------
      ! read command-line arguments
      ! -----------------------------------------

      IF (COMMAND_ARGUMENT_COUNT() .EQ. 4) THEN
            
            ! N: number of discretization points
            CALL GET_COMMAND_ARGUMENT(1, arg)
            READ(arg, *) Ntot

            ! k: number of eigenvalues
            CALL GET_COMMAND_ARGUMENT(2, arg)
            READ(arg, *) k_max

            ! omega: filename where to read the value of omega
            ! gives an error if input file does not exist
            CALL GET_COMMAND_ARGUMENT(3, omega_file)
            OPEN(unit=20, file=TRIM(omega_file), status='OLD')
                  READ(20,*) omega
            CLOSE(unit=20)

            ! outifle: filename where to write the eigenvalues and eigenvectors
            CALL GET_COMMAND_ARGUMENT(4, outfile)
            
      ELSE
            ! write error details into an error log file 
            PRINT *, "An error occurred, details on 'errors.log'"
            WRITE(10, "('Exactly 4 command-line arguments are expected:')")
            WRITE(10, "('N, k, omega_file, out_file')")
            STOP
      END IF

      ! checks on positivity
      IF ((Ntot .LE. 1) .OR. (k_max .LE. 0) .OR. (omega .LE. 0)) THEN
            ! write error details into an error log file
            PRINT *, "An error occurred, details on 'errors.log'"
            WRITE(10, "('Variables `N`, `k`, `omega` must be strictly positive')")
            STOP
      END IF

      ! checks on eigenvalue number
      IF (k_max .GT. Ntot) THEN
            ! write error details into an error log file
            PRINT *, "An error occurred, details on 'errors.log'"
            WRITE(10, "('Only a maximum of ', I0, ' eigenvalues can be calculated')") Ntot
            STOP
      END IF

      CALL checkpoint(debug, str="Number of function evaluation points = ", val=Ntot)
      CALL checkpoint(debug, str="Number of first eigenvalues to compute = ", val=k_max)
      CALL checkpoint(debug, str="Omega = ", val=omega)
      CALL checkpoint(debug, str="")



      ! ------------------------------------------------
      ! creation and filling of tridiagonal elements
      ! ------------------------------------------------
      
      ! grid parameters
      Ltot = 4.D0*SQRT(k_max/omega)
      deltax = Ltot/(Ntot+1.D0)
      CALL checkpoint(debug, str="Total length of the symulated system = ", val=Ltot)
      CALL checkpoint(debug, str="Spacing between points = ", val=deltax)

      ! allocation of dynamic arrays
      ALLOCATE(xgrid(Ntot))
      ALLOCATE(diag(Ntot), upper_diag(Ntot-1))
      ALLOCATE(eigenvals(Ntot), eigenvects(Ntot,k_max))
      ALLOCATE(work(5*Ntot), iwork(5*Ntot), ifail(Ntot))

      ! fill arrays that define tridiagonal matrix
      xgrid(:) = (/ (-0.5D0*Ltot + ii*deltax, ii=1,Ntot) /)
      diag(:) = (/ (2.D0/(deltax**2) + (omega*(-0.5D0*Ltot + ii*deltax))**2, ii=1,Ntot) /)
      upper_diag(:) = -1.D0/(deltax**2)



      ! -----------------------------------------
      ! selected eigenvalues and eigenvectors
      ! -----------------------------------------

      ! create directory where to store results
      CALL SYSTEM('mkdir -p results')

      CALL DSTEVX('V', 'I', Ntot, diag, upper_diag, 0.D0, 1.D0, 1, k_max, &
                  2*DLAMCH('S'), eig_num, eigenvals, eigenvects, Ntot, &
                  work, iwork, ifail, info)
      CALL checkpoint(debug, str="Info value = ", val=info)



      ! function normalization
      DO ii = 1, k_max
            norm2 = (SUM(eigenvects(1:Ntot:2, ii)**2)*4.D0 + SUM(eigenvects(2:Ntot:2, ii)**2)*2.D0)*deltax/3.D0
            eigenvects(:,ii) = eigenvects(:,ii) / SQRT(norm2)
      END DO

      
      ! -----------------------------------------
      ! writing results
      ! -----------------------------------------

      ! write data into the output file 
      OPEN(unit=30, file=outfile)

            ! print eigenvectors dimensions (N, k_max)
            WRITE(30, "(2(I0, :, ' '))", advance="no") Ntot, k_max
            WRITE(30, *)
            WRITE(30, *)


            ! print points grid
            DO ii = 1, Ntot
                  WRITE(30, "(ES24.17, ' ')", advance="no") xgrid(ii)
            END DO
            WRITE(30, *)
            WRITE(30, *)

            ! print eigenvalues and then eigenvector
            DO ii = 1, k_max
                  
                  WRITE(30, "(ES24.17, ' ')") eigenvals(ii)
                  DO jj = 1, Ntot
                        WRITE(30, "(ES24.17, ' ')", advance="no") eigenvects(jj,ii)
                  END DO
                  WRITE(30, *)
                  WRITE(30, *)
            END DO
      CLOSE(unit=30)



      ! -----------------------------------------
      ! final deallocations
      ! -----------------------------------------

      DEALLOCATE(xgrid)
      DEALLOCATE(diag, upper_diag)
      DEALLOCATE(eigenvals, eigenvects)
      DEALLOCATE(work, iwork, ifail)
      CALL checkpoint(debug, str="Deallocations perfomed successfully")


      ! close 'errors.log'
      CLOSE(unit=10)
END PROGRAM HarmOsc_1D