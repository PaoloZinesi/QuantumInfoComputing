! Program to simulate the time evolution of a wavefunction (in the ground state
! of the harmonic oscillator at t=0) subject to a quadratic potential which minimum
! point moves linearly in time along the positive direction.
! ........
! 
! Example of execution:
! ./TD_HarmOsc_1D
! will ....
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
      CHARACTER(LEN=100) :: arg, paramfile, initfile, outfile

      ! input information and parameters definition
      INTEGER :: Nx, Nt
      REAL*8 :: omega, Ttot, Ltot, xmin, deltat, deltax

      ! utility and loop variables
      REAL*8, PARAMETER :: PI = 4.D0*ATAN(1.D0)
      INTEGER :: t_idx

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

      IF (COMMAND_ARGUMENT_COUNT() .EQ. 5) THEN

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

            ! outifle: filename where to write results
            CALL GET_COMMAND_ARGUMENT(5, outfile)
            
      ELSE
            ! write error details into an error log file 
            PRINT *, "An error occurred, details on 'errors.log'"
            WRITE(10, "('Exactly 5 command-line arguments are expected:')")
            WRITE(10, "('Nt, omega_file, Ttot_file, init_file, out_file')")
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

      ! create and fill total wavefunction in the x domain
      psi_x = Zwavefunc(length=Nx, need_fftw_alloc=.TRUE.)
      CALL create_grid(psi_x, xmin=xmin, xmax=1.D0+ABS(xmin))
      psi_x%elem_fftw(1:(idx_minR - idx_minL + 1)) = psi0_x(idx_minL:idx_minR)
      psi_x%elem_fftw((idx_minR - idx_minL + 2):) = 0.D0

      ! normalize wavefunction
      psi_x%elem_fftw = psi_x%elem_fftw / SQRT(norm2_WF(psi_x))

      ! free memory used to store initial wavefunction
      DEALLOCATE(psi0_x, psi0_xgrid)

      ! create and fill total wavefunction in the p domain
      psi_p = Zwavefunc(length=Nx, need_fftw_alloc=.TRUE.)
      CALL create_grid(psi_p, xmin=0.D0, xmax=(2*PI*Nx)/Ltot)

      ! shift the grid of p in order to be compatible with both FFTW and physical
      ! conventions on the reciprocal-space array psi_p
      psi_p%grid(FLOOR((Nx - 1)/2.D0)+2:Nx) = psi_p%grid(FLOOR((Nx - 1)/2.D0)+2:Nx) - (2*PI*Nx)/Ltot



      ! -----------------------------------------
      ! time evolution
      ! -----------------------------------------

      ! vector to transform wavefunctions
      ALLOCATE(operator(Nx))




      ! WF declarations and allocations
      ! TYPE(Zwavefunc) :: psi_x, psi_p, psi_try
      ! psi_x = Zwavefunc(length=NN, need_fftw_alloc=.TRUE.)
      ! psi_p = Zwavefunc(length=NN, need_fftw_alloc=.TRUE.)
      ! psi_try = Zwavefunc(length=20, need_fftw_alloc=.TRUE.)

      ! ! create plans before filling the arrays
      ! CALL CPU_TIME(ti)
      ! plan_direct = fftw_plan_dft_1d(NN, psi_x%elem_fftw,psi_p%elem_fftw, FFTW_FORWARD,FFTW_MEASURE)
      ! plan_inverse = fftw_plan_dft_1d(NN, psi_p%elem_fftw,psi_x%elem_fftw, FFTW_BACKWARD,FFTW_MEASURE)
      ! CALL CPU_TIME(tf)
      ! CALL checkpoint(debug, str="Time [s] necessary to produce both plans =", val=tf-ti)

      ! ! fill array grids
      ! CALL create_grid(psi_x, xmin=0.D0, xmax=LL)
      ! CALL create_grid(psi_p, xmin=0.D0, xmax=(2*PI*NN)/LL)

      ! ! shift the grid of p in order to be compatible with both FFTW and physical
      ! ! conventions on the reciprocal-space array psi_p
      ! psi_p%grid(FLOOR((NN - 1)/2.D0)+2:NN) = psi_p%grid(FLOOR((NN - 1)/2.D0)+2:NN) - (2*PI*NN)/LL

      ! ! fill array
      ! sigmax = LL/20.D0
      ! psi_x%elem_fftw = (/ (  EXP(-0.5D0*((psi_x%grid(ii) - (LL/2.4D0))/sigmax)**2.D0) * &
      !                         EXP(COMPLEX(0,(2*PI/LL)*10*psi_x%grid(ii))), ii=1,NN) /)


      ! CALL checkpoint(debug, str="")


      ! ! write wavefunction into file
      ! CALL SYSTEM('mkdir -p results')
      ! CALL writeWFFile(psi_x, unit=10, file="results/init.dat", format="ES24.17")
      ! CALL checkpoint(debug, str="Norm^2 of psi_x =", val=norm2_WF(psi_x))
      ! CALL checkpoint(debug, str="")


      ! ! direct FFT
      ! CALL CPU_TIME(ti)
      ! CALL fftw_execute_dft(plan_direct, psi_x%elem_fftw, psi_p%elem_fftw)
      ! CALL CPU_TIME(tf)
      ! CALL checkpoint(debug, str="Time [s] necessary to execute direct plan =", val=tf-ti)

      ! CALL writeWFFile(psi_p, unit=20, file="results/FFT_init.dat", format="ES24.17")
      ! CALL checkpoint(debug, str="Norm^2 of psi_p =", val=norm2_WF(psi_p))
      ! CALL checkpoint(debug, str="")


      ! ! inverse FFT
      ! CALL CPU_TIME(ti)
      ! CALL fftw_execute_dft(plan_inverse, psi_p%elem_fftw, psi_x%elem_fftw)
      ! CALL CPU_TIME(tf)
      ! CALL checkpoint(debug, str="Time [s] necessary to execute inverse plan =", val=tf-ti)

      ! CALL writeWFFile(psi_x, unit=30, file="results/FFT_inverse.dat", format="ES24.17")
      ! CALL checkpoint(debug, str="Norm^2 of psi_x (after FT and inverse FT) =", val=norm2_WF(psi_x))
      ! CALL checkpoint(debug, str="")


      ! ! destroy plans and free WF memory
      ! CALL fftw_destroy_plan(plan_direct)
      ! CALL fftw_destroy_plan(plan_inverse)



      CALL freeWF(psi_x)
      CALL freeWF(psi_p)

      ! close 'errors.log'
      CLOSE(unit=10)

END PROGRAM TD_HarmOsc_1D