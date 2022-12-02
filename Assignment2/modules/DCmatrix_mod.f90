! Module to define a derived DCmatrix type, which contains a double complex
! matrix and a two-dimensional vector to store its dimensions. The DCmatrix needs
! to be allocated with an initializer and must be DEALLOCATED MANUALLY in the code.
! Function to compute the adjoint matrix and the matrix trace are provided, while
! a subroutine to write the matrix into a file is also provided.
! 
! Author: Paolo Zinesi
!

MODULE DCmatrix_mod
      USE checkpoint_mod
      IMPLICIT NONE

      ! enable or disable debugging for this module ONLY
      LOGICAL, PRIVATE :: debug = .TRUE.

      ! declaration of derived type
      TYPE DCmatrix
            ! matrix dimensions
            INTEGER, DIMENSION(2) :: dim  

            ! matrix elements
            DOUBLE COMPLEX, DIMENSION(:,:), ALLOCATABLE :: elem 
      END TYPE DCmatrix


      ! initializer interface
      INTERFACE DCmatrix
            MODULE PROCEDURE DCmatrix_init
      END INTERFACE DCmatrix

      ! adjoint operator interface
      INTERFACE OPERATOR(.ADJ.)
            MODULE PROCEDURE MatAdjoint
      END INTERFACE OPERATOR(.ADJ.)

      ! trace operator interface
      INTERFACE OPERATOR(.TR.)
            MODULE PROCEDURE MatTrace
      END INTERFACE OPERATOR(.TR.)



      CONTAINS


      ! This function returns a new DCmatrix already allocated, so that
      ! in the main program there is no need to manually allocate every DCmatrix.
      ! This initializer can be used to allocate both uninitialized matrices
      ! (if "val" is absent) and initialized matrices with a constant value ("val").
      !
      ! inputs:
      ! - rows [integer]: number of matrix rows
      ! - cols [integer]: number of matrix columns
      ! - val [double complex, optional]: value to fill matrix
      !
      ! outputs:
      ! - M [DCmatrix]: initialized DCmatrix
      ! 
      ! TODO:
      !     create a functionality to read a matrix from file (maybe using a dedicated
      !     subroutine) and initialize a new DCmatrix
      ! 
      TYPE(DCmatrix) FUNCTION DCmatrix_init(rows, cols, val) RESULT(M)
            IMPLICIT NONE

            INTEGER, INTENT(IN)                       :: rows, cols
            DOUBLE COMPLEX, INTENT(IN), OPTIONAL      :: val

            ! check if rows and columns are positive
            IF((rows .LE. 0) .OR. (cols .LE. 0)) THEN
                  ! halt execution
                  PRINT *, "'DCmatrix_init': Non-positive number given as matrix dimension"
                  STOP
            END IF

            ! initialization of rows and columns, and allocation
            M%dim = (/ rows, cols /)
            ALLOCATE(M%elem(M%dim(1),M%dim(2)))

            ! filling with "val" if present
            ! otherwise leaves it unspecified, but the user must fill the matrix itself
            IF(PRESENT(val)) THEN
                  M%elem = val
            ENDIF
            
            CALL checkpoint(debug, str="Created and allocated matrix")

      END FUNCTION DCmatrix_init



      ! This function returns a new DCmatrix already allocated that is the adjoint
      ! of a given input matrix. The initilizer previously defined is used.
      ! For this function the OPERATOR(.ADJ.) is defined through an interface.
      !
      ! inputs:
      ! - M [DCmatrix]: input DCmatrix to transpose
      !
      ! outputs:
      ! - M_Adj [DCmatrix]: transposed DCmatrix
      ! 
      ! TODO: None
      ! 
      TYPE(DCmatrix) FUNCTION MatAdjoint(M) RESULT(M_Adj)
            IMPLICIT NONE

            TYPE(DCmatrix), INTENT(IN) :: M

            ! utility and loop variables
            INTEGER :: ii, jj
            CHARACTER(LEN=1024) :: str

            CALL checkpoint(debug, str=NEW_LINE("")//"Entering 'MatAdjoint' function")

            ! print dimensions debug string
            WRITE (str, "('DIM(M)     =',I2,I2)") M%dim(1),M%dim(2)
            CALL checkpoint(debug, str=TRIM(str))


            ! allocate new matrix with transposed dimensions
            M_Adj = DCmatrix(rows=M%dim(2), cols=M%dim(1))



            ! print dimensions debug string
            WRITE (str, "('DIM(M_Adj) =',I2,I2)") M_Adj%dim(1),M_Adj%dim(2)
            CALL checkpoint(debug, str=TRIM(str))



            ! loop over rows and columns of original matrix
            DO ii = 1, M%dim(1)
                  DO jj = 1, M%dim(2)
                        ! assign values to the adjoint matrix
                        M_Adj%elem(jj,ii) = CONJG(M%elem(ii,jj))
                  END DO
            END DO
            
            CALL checkpoint(debug, str="Exiting 'MatAdjoint' function"//NEW_LINE(""))
      END FUNCTION MatAdjoint



      ! This function returns the trace of a square DCmatrix given as input.
      ! If the input is not a square DCmatrix, an error is printed and the 
      ! program is stopped.
      !
      ! inputs:
      ! - M [DCmatrix]: input square DCmatrix to compute the trace of
      !
      ! outputs:
      ! - trace [double complex]: trace of input matrix
      ! 
      ! TODO: None
      ! 
      DOUBLE COMPLEX FUNCTION MatTrace(M) RESULT(trace)
            IMPLICIT NONE

            TYPE(DCmatrix), INTENT(IN) :: M

            ! utility and loop variables
            INTEGER :: ii
            CHARACTER(LEN=1024) :: str
            trace = COMPLEX(0.D0, 0.D0)

            CALL checkpoint(debug, str=NEW_LINE("")//"Entering 'MatTrace' function")

            ! check dimensions
            IF(M%dim(1) .EQ. M%dim(2)) THEN

                  DO ii = 1, M%dim(1)
                        trace = trace + M%elem(ii,ii)
                  END DO

            ELSE
                  ! print dimensions debug string
                  WRITE (str, "('DIM(M):',I2,' !=',I2)") M%dim(1), M%dim(2)
                  CALL checkpoint(debug, str=TRIM(str))

                  ! halt execution
                  PRINT *, "'MatTrace': Trying to compute trace of a non-square matrix"
                  STOP
            END IF

            CALL checkpoint(debug, str="Exiting 'MatTrace' function"//NEW_LINE(""))

      END FUNCTION MatTrace


      ! This subroutine prints on a file a DCmatrix M in readable form.
      ! Unit, filename, and output format are given as input for flexibility.
      ! Only units different from 0,5,6 can be used to write into file (otherwise 
      ! collisions with default writing units arise).
      !
      ! inputs:
      ! - M [DCmatrix]: input DCmatrix to write into a file
      ! - unit [integer]: identifier of the unit connected to the file
      ! - file [character(len=*)]: filename where matrix has to be written
      ! - format [character(len=*)]: string to specify how to format written numbers.
      !                              A suggestion is to use format="(ES12.5,SP,ES13.5,' i   ')"
      !                              to have readable complex numbers in scientific notation
      !                              with 5 decimals and the 'i' at the end of the imag. part.
      !                              
      ! outputs: None
      ! 
      ! TODO: None
      ! 
      SUBROUTINE writeMatFile(M, unit, file, format)
            IMPLICIT NONE

            TYPE(DCmatrix), INTENT(IN) :: M
            INTEGER :: unit
            CHARACTER (LEN=*) :: file, format

            ! utility and loop variables
            INTEGER :: ii, jj

            CALL checkpoint(debug, str=NEW_LINE("")//"Entering 'writeMatFile' subroutine")

            ! do not use units dedicated to standard fortran I/O (PRINT * for ex.)
            IF((unit .EQ. 0) .OR. (unit .EQ. 5) .OR. (unit .EQ. 6)) THEN
                  PRINT "('Not allowed to use unit', I2)", unit
                  PRINT "(A, ' has not been written')", file

            ELSE
                  ! write data into the file 
                  OPEN(unit=unit, file=file)
                        DO ii = 1, M%dim(1)
                              DO jj = 1, M%dim(2)
                                    ! advance='no' to remain in the present line                        
                                    WRITE(unit, format, advance="no") M%elem(ii,jj)
                              END DO

                              ! new line
                              WRITE(unit,*)
                        END DO
                  CLOSE(unit=unit)
            END IF

            CALL checkpoint(debug, str="Exiting 'writeMatFile' subroutine"//NEW_LINE(""))
      END SUBROUTINE writeMatFile


      ! TODO:
      ! create a subroutine to check equality between two DCmatrices (with a given tolerance)


END MODULE DCmatrix_mod