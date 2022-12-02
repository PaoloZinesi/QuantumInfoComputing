! Module to define a derived Zwavefunc type, which contains a double complex
! vector to store the wavefunction, a double array to store the grid on which the
! wavefunction is defined, and a variable to store the length of the wavefunction.
! The Zwavefunc needs to be allocated with an initializer and must be deallocated 
! using the dedicated subroutine 'freeWF'.
! 
! Author: Paolo Zinesi
!

MODULE Zwavefunc_mod
      USE checkpoint_mod
      USE, INTRINSIC :: iso_c_binding
      USE FFTW_mod
      IMPLICIT NONE

      ! enable or disable debugging for this module ONLY
      LOGICAL, PRIVATE :: debug = .TRUE.

      ! declaration of derived type
      TYPE Zwavefunc
            ! wavefunction length
            INTEGER :: len = 0

            ! wavefunction elements with normal allocation
            DOUBLE COMPLEX, DIMENSION(:), ALLOCATABLE :: elem

            ! elements with FFTW allocation
            COMPLEX(C_DOUBLE_COMPLEX), POINTER :: elem_fftw(:)
            TYPE(C_PTR) :: ptr
            LOGICAL :: useFFT = .FALSE.
            
            ! grid
            REAL*8, DIMENSION(:), ALLOCATABLE :: grid
            REAL*8 :: deltax = 0.D0

      END TYPE Zwavefunc


      ! initializer interface
      INTERFACE Zwavefunc
            MODULE PROCEDURE Zwavefunc_init
      END INTERFACE Zwavefunc


      CONTAINS
      


      ! This function returns a new Zwavefunc already allocated, with allocations
      ! specified in input. If not explicitly specified, the array is initialized
      ! in the usual fortran way. If the allocation for FFTW3 is requested, the
      ! vector is aligned in memory to ensure the best performances.
      !
      ! inputs:
      ! - length [integer]: number of elements to allocate
      ! - need_fftw_alloc [logical]: boolean value to switch on FFTW3 alignment
      !
      ! outputs:
      ! - WF [Zwavefunc]: initialized Zwavefunc
      ! 
      ! TODO: None
      ! 
      TYPE(Zwavefunc) FUNCTION Zwavefunc_init(length, need_fftw_alloc) RESULT(WF)
            IMPLICIT NONE

            INTEGER, INTENT(IN) :: length
            LOGICAL, INTENT(IN), OPTIONAL :: need_fftw_alloc


            ! check if length is positive
            IF(length .LE. 0) THEN
                  ! halt execution
                  PRINT *, "'Zwavefunc_init': Non-positive number given as wavefunction dimension"
                  STOP
            END IF


            ! allocate vector according to the desired allocation.
            ! unwanted allocations are set to null.
            IF(PRESENT(need_fftw_alloc) .AND. (need_fftw_alloc .EQV. .TRUE.)) THEN
                  ! when using FFTW, a special allocation is required to obtain the maximum
                  ! efficiency in the FFTW routines
                  WF%ptr = fftw_alloc_complex(int(length, C_SIZE_T))
                  CALL c_f_pointer(WF%ptr, WF%elem_fftw, [length])
                  WF%useFFT = .TRUE.

                  CALL checkpoint(debug, str="Vector allocated using FFTW allocation")
            ELSE
                  ! allocate in the usual FORTRAN way
                  ALLOCATE(WF%elem(length))
                  WF%useFFT = .FALSE.

                  ! set to null the C-type array and pointer
                  WF%ptr = c_null_ptr
                  WF%elem_fftw => NULL()
                  
                  CALL checkpoint(debug, str="Vector allocated using usual FORTRAN allocation")
            END IF


            WF%len = length
            ALLOCATE(WF%grid(length))

      END FUNCTION Zwavefunc_init



      ! This subroutine allocate and fills the grid of the Zwavefunc using
      ! the values of start point (first element of the array) and end of the grid
      ! (not included in grid). The number of points is the same as the vector length.
      ! The value of deltax is automatically computed and stored.
      !
      ! inputs:
      ! - WF [Zwavefunc]: input Zwavefunc to create the grid of
      ! - xmin [REAL*8]: first element of the grid
      ! - xmax [REAL*8]: last element of the grid (not included)
      !                              
      ! outputs: None
      ! 
      ! TODO: None
      ! 
      SUBROUTINE create_grid(WF, xmin, xmax)
            IMPLICIT NONE

            TYPE(Zwavefunc), INTENT(INOUT) :: WF
            REAL*8, INTENT(IN) :: xmin, xmax

            ! utility and loop variables
            INTEGER :: ii

            CALL checkpoint(debug, str="Entering 'create_grid' subroutine")


            ! check for validity of grid limits, xmin < xmax !
            IF(xmin .GE. xmax) THEN
                  PRINT "('`create_grid`: lower grid value is greater than the upper grid value')"
                  STOP
            END IF

            ! check if length is positive
            IF(WF%len .LE. 0) THEN
                  PRINT "('`create_grid`: nonpositive length of the grid')"
                  STOP
            END IF


            ! define grid spacing
            WF%deltax = (xmax-xmin)/(WF%len)

            ! define equally spaced grid
            WF%grid = (/ (xmin + ii * WF%deltax, ii = 0, WF%len - 1) /) 


            CALL checkpoint(debug, str="Exiting 'create_grid' subroutine")
      END SUBROUTINE create_grid



      ! This subroutine prints on a file a Zwavefunc WF.
      ! Unit, filename, and output format are given as input for flexibility.
      ! Only units different from 0,5,6 can be used to write into file (otherwise 
      ! collisions with default writing units arise).
      !
      ! inputs:
      ! - WF [Zwavefunc]: input Zwavefunc to write into a file
      ! - unit [integer]: identifier of the unit connected to the file
      ! - file [character(len=*)]: filename where matrix has to be written
      ! - format [character(len=*)]: string to specify how to format written numbers.
      !                              A suggestion is to use format="ES24.17".
      !                              
      ! outputs: None
      ! 
      ! TODO: None
      ! 
      SUBROUTINE writeWFFile(WF, unit, file, format)
            IMPLICIT NONE

            TYPE(Zwavefunc), INTENT(IN) :: WF
            INTEGER :: unit
            CHARACTER (LEN=*) :: file, format

            ! loop variable
            INTEGER :: ii

            CALL checkpoint(debug, str="Entering 'writeWFFile' subroutine")

            ! do not use units dedicated to standard fortran I/O (PRINT * for ex.)
            IF((unit .EQ. 0) .OR. (unit .EQ. 5) .OR. (unit .EQ. 6)) THEN
                  PRINT "('Not allowed to use unit', I2)", unit
                  PRINT "(A, ' has not been written')", file

            ELSE
                  ! write wavefunction into the file 
                  OPEN(unit=unit, file=file)

                        ! write the grid
                        DO ii = 1, WF%len
                              ! advance='no' to remain in the present line                        
                              WRITE(unit, "("//format//", ' ')", advance="no") WF%grid(ii)
                        END DO

                        WRITE(unit,*)
                        WRITE(unit,*)

                        IF(WF%useFFT .EQV. .FALSE.) THEN
                              ! write the function allocated using usual FORTRAN
                              DO ii = 1, WF%len                       
                                    WRITE(unit, "("//format//", ' ',"//format//")") WF%elem(ii)
                              END DO
                        ELSE
                              ! write the function allocated using FFTW3 alignment
                              DO ii = 1, WF%len
                                    WRITE(unit, "("//format//", ' ',"//format//")") WF%elem_fftw(ii)
                              END DO
                        END IF
                        

                  CLOSE(unit=unit)
            END IF

            CALL checkpoint(debug, str="Exiting 'writeWFFile' subroutine")
      END SUBROUTINE writeWFFile



      ! This subroutine frees the space allocated when creating the wavefunction.
      !
      ! inputs:
      ! - WF [Zwavefunc]: input Zwavefunc to deallocate
      !                              
      ! outputs: None
      ! 
      ! TODO: None
      ! 
      SUBROUTINE freeWF(WF)
            IMPLICIT NONE

            TYPE(Zwavefunc), INTENT(INOUT) :: WF

            ! deallocate wavefunction arrays
            IF(WF%useFFT .EQV. .FALSE.) THEN
                  DEALLOCATE(WF%elem)
            END IF
            
            CALL fftw_free(WF%ptr)

            ! deallocate grid
            DEALLOCATE(WF%grid)

            CALL checkpoint(debug, str="Deallocation performed successfully.")


      END SUBROUTINE freeWF



      ! This function returns the integral of the modulus squared of the wavefunction
      ! using Simpson's rule. The integral is computed according to the deltax parameter
      ! stored in Zwavefunc composite type.
      !
      ! inputs:
      ! - WF [Zwavefunc]: input Zwavefunc to compute the integral of
      !
      ! outputs:
      ! - norm2 [REAL*8]: integral of |WF|^2 in dx using Simpson's rule
      ! 
      ! TODO: add the possibility of computing a weighted integral (useful for average value)
      ! 
      REAL*8 FUNCTION norm2_WF(WF) RESULT(norm2)
            IMPLICIT NONE

            TYPE(Zwavefunc), INTENT(IN) :: WF
            INTEGER :: istart, iend


            ! check if length is positive
            IF(WF%len .LE. 0) THEN
                  PRINT "('`norm2_WF`: nonpositive length of the grid')"
                  STOP
            END IF

            ! check if wavefunction has deltax > 0
            IF(WF%deltax .LE. 0.D0) THEN
                  PRINT "('`norm2_WF`: invalid deltax')"
                  STOP
            END IF

            norm2 = 0.D0
            ! implement Simpson's rule
            IF(MOD(WF%len,2) .EQ. 1) THEN

                  istart = 1
                  iend = WF%len
                  ! odd case: a single sum is computed
                  IF(WF%useFFT .EQV. .FALSE.) THEN
                        norm2 =     norm2 + &
                                    (WF%deltax/3.D0)*(ABS(WF%elem(istart))**2 + ABS(WF%elem(iend))**2 + &
                                                4*SUM(ABS(WF%elem(istart+1:iend-1:2))**2) + &
                                                2*SUM(ABS(WF%elem(istart+2:iend-2:2))**2))
                  ELSE
                        norm2 =     norm2 + &
                                    (WF%deltax/3.D0)*(ABS(WF%elem_fftw(istart))**2 + ABS(WF%elem_fftw(iend))**2 + &
                                                4*SUM(ABS(WF%elem_fftw(istart+1:iend-1:2))**2) + &
                                                2*SUM(ABS(WF%elem_fftw(istart+2:iend-2:2))**2))
                  END IF
            ELSE
                  ! even case: two sums are computed and averaged together

                  ! first sum over [1 : N-1]
                  istart = 1
                  iend = WF%len - 1
                  IF(WF%useFFT .EQV. .FALSE.) THEN
                        norm2 =     norm2 + &
                                    0.5D0*(WF%deltax/3.D0)*(ABS(WF%elem(istart))**2 + ABS(WF%elem(iend))**2 + &
                                                4*SUM(ABS(WF%elem(istart+1:iend-1:2))**2) + &
                                                2*SUM(ABS(WF%elem(istart+2:iend-2:2))**2))
                  ELSE
                        norm2 =     norm2 + &
                                    0.5D0*(WF%deltax/3.D0)*(ABS(WF%elem_fftw(istart))**2 + ABS(WF%elem_fftw(iend))**2 + &
                                                4*SUM(ABS(WF%elem_fftw(istart+1:iend-1:2))**2) + &
                                                2*SUM(ABS(WF%elem_fftw(istart+2:iend-2:2))**2))
                  END IF

                  ! second sum over [2 : N]
                  istart = 2
                  iend = WF%len
                  ! even case: two sums are computed and averaged together
                  IF(WF%useFFT .EQV. .FALSE.) THEN
                        norm2 =     norm2 + &
                                    0.5D0*(WF%deltax/3.D0)*(ABS(WF%elem(istart))**2 + ABS(WF%elem(iend))**2 + &
                                                4*SUM(ABS(WF%elem(istart+1:iend-1:2))**2) + &
                                                2*SUM(ABS(WF%elem(istart+2:iend-2:2))**2))
                  ELSE
                        norm2 =     norm2 + &
                                    0.5D0*(WF%deltax/3.D0)*(ABS(WF%elem_fftw(istart))**2 + ABS(WF%elem_fftw(iend))**2 + &
                                                4*SUM(ABS(WF%elem_fftw(istart+1:iend-1:2))**2) + &
                                                2*SUM(ABS(WF%elem_fftw(istart+2:iend-2:2))**2))
                  END IF
            END IF

      END FUNCTION norm2_WF



      ! TODO:
      ! create a subroutine to read wavefunction from file


END MODULE Zwavefunc_mod