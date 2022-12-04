! Program to simulate the time evolution of a wavefunction (in the ground state
! of the harmonic oscillator at t=0) subject to a quadratic potential which minimum
! point moves linearly in time along the positive direction.
! The parameters needed by the program are given in the command-line and the results are
! printed on file to produce the output plots.
! 
! Example of execution:
! ./TD_HarmOsc_1D.out 100 input/omega.dat input/Ttot.dat input/init_WF.dat
! will evolve the wavefunction from t=0 to t=Ttot through 100 time steps, using the value of omega found
! in 'input/omega.dat', the value of Ttot found in 'input/Ttot.dat' and the initialization wavefunction
! found in 'input/init_WF.dat'. Results are created in the 'results' folder.
!
! Author: Paolo Zinesi
!

! main program
PROGRAM TD_HarmOsc_1D
      USE checkpoint_mod
      USE, INTRINSIC :: iso_c_binding
      USE FFTW_mod
      USE Zwavefunc_mod
      IMPLICIT NONE
      LOGICAL :: debug = .TRUE.

      ! variables to read command-line arguments
      CHARACTER(LEN=100) :: arg, paramfile, initfile

      ! input information and parameters definition
      INTEGER :: Nx, Nt
      REAL*8 :: omega, Ttot, Ltot, xmin, deltat, deltax

      ! utility and loop variables
      REAL*8, PARAMETER :: PI = 4.D0*ATAN(1.D0)
      CHARACTER(LEN=100) :: str
      INTEGER :: print_idx, t_idx

      ! wavefunctions
      INTEGER :: Ninit, idx_minL, idx_minR
      REAL*8, DIMENSION(:), ALLOCATABLE :: psi0_x, psi0_xgrid
      TYPE(Zwavefunc) :: psi_x, psi_p
      DOUBLE COMPLEX, DIMENSION(:), ALLOCATABLE :: operator

      ! FFTW plans
      type(C_PTR) :: plan_direct, plan_inverse


      ! open 'errors.log'
      OPEN(unit=10, file='errors.log', access='APPEND')
      


      ! -----------------------------------------
      ! read command-line arguments
      ! -----------------------------------------

      IF (COMMAND_ARGUMENT_COUNT() .EQ. 4) THEN

            ! Nt: number of discretization points in time
            CALL GET_COMMAND_ARGUMENT(1, arg)
            READ(arg, *) Nt

            ! omega: filename where to read the value of omega
            ! gives an error if input file does not exist
            ! must be the same used to generate the initial wavefunction
            CALL GET_COMMAND_ARGUMENT(2, paramfile)
            OPEN(unit=20, file=TRIM(paramfile), status='OLD')
                  READ(20,*) omega
            CLOSE(unit=20)

            ! Ttot: filename where to read the value of T (maximum time)
            ! gives an error if input file does not exist
            CALL GET_COMMAND_ARGUMENT(3, paramfile)
            OPEN(unit=21, file=TRIM(paramfile), status='OLD')
                  READ(21,*) Ttot
            CLOSE(unit=21)

            ! initfile: filename where to read the initial wavefunction
            CALL GET_COMMAND_ARGUMENT(4, initfile)
            
      ELSE
            ! write error details into an error log file 
            PRINT *, "An error occurred, details on 'errors.log'"
            WRITE(10, "('Exactly 4 command-line arguments are expected:')")
            WRITE(10, "('Nt, omega_file, Ttot_file, init_file')")
            STOP
      END IF

      ! checks on positivity
      IF ((Nt .LE. 1) .OR. (omega .LE. 0) .OR. (Ttot .LE. 0)) THEN
            ! write error details into an error log file
            PRINT *, "An error occurred, details on 'errors.log'"
            WRITE(10, "('Variables `Nt`, `omega`, `Ttot` must be strictly positive')")
            STOP
      END IF

      CALL checkpoint(debug, str="Nt = ", val=Nt)
      CALL checkpoint(debug, str="Omega = ", val=omega)
      CALL checkpoint(debug, str="Ttot = ", val=Ttot)
      CALL checkpoint(debug, str="")



      ! -------------------------------------------------------
      ! read the ground-state eigenfunction from 'initfile'
      ! -------------------------------------------------------
      OPEN(unit=30, file=initfile)

            ! read length of array and allocate vector
            READ(30,*) Ninit

            ALLOCATE(psi0_x(Ninit), psi0_xgrid(Ninit))

            ! read grid on which the wavefunction is defined
            READ(30,*) psi0_xgrid

            ! read first eigenvalue (to be ignored)
            READ(30,*) psi0_x(1)

            ! read the first wavefunction
            READ(30,*) psi0_x

      CLOSE(unit=30)


      ! grid spacing of the initial wavefunction
      deltax = (psi0_xgrid(Ninit) - psi0_xgrid(1))/(Ninit-1)

      ! find indices of init WF to transfer to the total system 
      idx_minL = MAXLOC(psi0_x, MASK=psi0_x<(1.0D-10)*MAXVAL(psi0_x) .AND. (psi0_xgrid .LT. 0.D0), DIM=1)
      idx_minR = MAXLOC(psi0_x, MASK=psi0_x<(1.0D-10)*MAXVAL(psi0_x) .AND. (psi0_xgrid .GT. 0.D0), DIM=1)
      xmin = psi0_xgrid(idx_minL)


      ! Nx, Ltot depends directly on xmin
      Ltot = 1.D0 + 2.D0*ABS(xmin)
      Nx = CEILING(Ltot/deltax)

      CALL checkpoint(debug, str="Ltot = ", val=Ltot)
      CALL checkpoint(debug, str="Nx = ", val=Nx)
      CALL checkpoint(debug, str="Number of new points added to the grid = ", val=Nx - (idx_minR - idx_minL + 1))
      CALL checkpoint(debug, str="")



      ! -----------------------------------------
      ! definition of total wavefunctions
      ! -----------------------------------------

      ! create total wavefunctions in the x and p domain
      psi_x = Zwavefunc(length=Nx, need_fftw_alloc=.TRUE.)
      psi_p = Zwavefunc(length=Nx, need_fftw_alloc=.TRUE.)


      ! create FFTW plans before filling the arrays (this procedure might change the array values)
      plan_direct = fftw_plan_dft_1d(Nx, psi_x%elem_fftw,psi_p%elem_fftw, FFTW_FORWARD,FFTW_MEASURE)
      plan_inverse = fftw_plan_dft_1d(Nx, psi_p%elem_fftw,psi_x%elem_fftw, FFTW_BACKWARD,FFTW_MEASURE)
      CALL checkpoint(debug, str="FFTW plans created.")


      ! fill wavefunctions in the x domain
      CALL create_grid(psi_x, xmin=xmin, xmax=1.D0+ABS(xmin))
      psi_x%elem_fftw(1:(idx_minR - idx_minL + 1)) = psi0_x(idx_minL:idx_minR)
      psi_x%elem_fftw((idx_minR - idx_minL + 2):) = 0.D0

      ! normalize wavefunction
      psi_x%elem_fftw = psi_x%elem_fftw / SQRT(norm2_WF(psi_x))

      ! free memory used to store initial wavefunction
      DEALLOCATE(psi0_x, psi0_xgrid)

      ! fill total wavefunction in the p domain
      CALL create_grid(psi_p, xmin=0.D0, xmax=(2*PI*Nx)/Ltot)

      ! shift the grid of p in order to be compatible with both FFTW and physical
      ! conventions on the reciprocal-space array psi_p
      psi_p%grid(FLOOR((Nx - 1)/2.D0)+2:Nx) = psi_p%grid(FLOOR((Nx - 1)/2.D0)+2:Nx) - (2*PI*Nx)/Ltot



      ! -----------------------------------------
      ! time evolution
      ! -----------------------------------------

      CALL SYSTEM("mkdir -p results")
      CALL checkpoint(debug, str="Started computing time evolution")

      deltat = Ttot/Nt

      ! vector to transform wavefunctions
      ALLOCATE(operator(Nx))

      print_idx = 1
      DO t_idx = 1, Nt
            ! V/2 propagation
            operator = EXP(COMPLEX(0.D0,-1.D0)*deltat*(omega**2)*0.5D0*(psi_x%grid - (t_idx*1.D0)/Nt)**2)
            psi_x%elem_fftw = psi_x%elem_fftw * operator

            ! Fourier Transform
            CALL fftw_execute_dft(plan_direct, psi_x%elem_fftw, psi_p%elem_fftw)

            ! T propagation
            operator = EXP(COMPLEX(0.D0,-1.D0)*deltat*(psi_p%grid)**2)
            psi_p%elem_fftw = psi_p%elem_fftw * operator

            ! Inverse Fourier Transform
            CALL fftw_execute_dft(plan_inverse, psi_p%elem_fftw, psi_x%elem_fftw)
            psi_x%elem_fftw = psi_x%elem_fftw / Nx

            ! V/2 propagation
            operator = EXP(COMPLEX(0.D0,-1.D0)*deltat*(omega**2)*0.5D0*(psi_x%grid - (t_idx*1.D0)/Nt)**2)
            psi_x%elem_fftw = psi_x%elem_fftw * operator

            ! print wavefunctions on file
            IF(MOD(t_idx, 10) .EQ. 0) THEN
                  WRITE(str, "('results/outWF_', I0, '.dat')") t_idx
                  CALL writeWFFile(psi_x, unit=30+print_idx, file=TRIM(str), format="ES24.17")
                  print_idx = print_idx + 1

                  PRINT *, norm2_WF(psi_x, weights=psi_x%grid)
            END IF

      END DO

      CALL checkpoint(debug, str="Finished computing time evolution")


      ! free wavefunctions
      CALL freeWF(psi_x)
      CALL freeWF(psi_p)

      ! close 'errors.log'
      CLOSE(unit=10)

END PROGRAM TD_HarmOsc_1D