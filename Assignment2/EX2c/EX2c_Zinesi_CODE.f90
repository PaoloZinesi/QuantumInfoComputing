! Program to test the DCmatrix_mod module
! This program takes as input from the command line the size of the DCmatrix
! to build, then it creates and fills the DCmatrix with an arbitrary procedure.
! Then it computes the adjoint matrix. All the matrices created are 
! printed in an output file and the trace is finally computed and returned.
! 
! Example of execution:
! ./a.out 20 20
! will create a 20x20 DCmatrix, calculate its adoint DCmatrix and its trace.
! The program is expected to halt at the trace computation if a non-square
! matrix is given as input.
!
! Author: Paolo Zinesi
!

! main program
PROGRAM test_DCmatrix
      USE checkpoint_mod
      USE DCmatrix_mod
      IMPLICIT NONE
      LOGICAL :: debug = .TRUE.

      ! utility and loop variables
      INTEGER :: ii, jj
      CHARACTER(LEN=1024) :: str

      ! variables to read command-line arguments
      CHARACTER(LEN=32) :: arg
      INTEGER :: arg_int
      INTEGER :: arg_idx

      ! DCmatrix variables
      INTEGER, DIMENSION(2) :: input_sizes
      TYPE(DCmatrix) :: M_in, M_in_Adj
      DOUBLE COMPLEX :: trace, trace_adj


      ! -----------------------------------------
      ! read command-line arguments
      ! -----------------------------------------

      IF (COMMAND_ARGUMENT_COUNT() .EQ. 2) THEN
            DO arg_idx = 1, 2

                  ! get argument number "arg_idx" and put it into "arg_int"
                  CALL GET_COMMAND_ARGUMENT(arg_idx, arg)
                  READ(arg, *) arg_int

                  input_sizes(arg_idx) = arg_int

            END DO
            
      ELSE
            PRINT *, "Exactly 2 command-line arguments are expected"
            STOP
      END IF
      CALL checkpoint(debug, str="Input rows = ", val=input_sizes(1))
      CALL checkpoint(debug, str="Input columns = ", val=input_sizes(2))



      ! -----------------------------------------
      ! creation and filling of matrices
      ! -----------------------------------------
      ! important checks on matrix dimensions
      ! positive dimensions
      DO ii = 1, 2
            IF(input_sizes(ii) .LE. 0) THEN
                  PRINT *, "Non-positive number given as matrix dimension"
                  STOP
            END IF
      END DO
      CALL checkpoint(debug, str="Matrix dimensions are positive")


      ! creation of the desired matrix
      ! M_in is now allocated
      M_in = DCmatrix(rows=input_sizes(1), cols=input_sizes(2))


      ! (arbitrary) filling of matrix
      DO ii = 1, M_in%dim(1)
            DO jj = 1, M_in%dim(2)
                  M_in%elem(ii,jj) = ii * jj * COMPLEX(2.0D-02, 1.0D-03)
            END DO
      END DO



      ! -----------------------------------------
      ! computing of the adjoint matrix
      ! -----------------------------------------
      M_in_Adj = .ADJ.(M_in)
      CALL checkpoint(debug, str="Adjoint matrix rows = ", val=M_in_Adj%dim(1))
      CALL checkpoint(debug, str="Adjoint matrix columns = ", val=M_in_Adj%dim(2))



      ! -----------------------------------------
      ! writing matrices on file
      ! -----------------------------------------
      ! write input matrix
      WRITE (str, "('M_',I0.0,'_',I0.0,'.dat')") M_in%dim(1), M_in%dim(2)
      CALL writeMatFile(M_in, unit=10, file=str, format ="(ES12.5,SP,ES13.5,' i   ')")

      ! write adjoint matrix
      WRITE (str, "('M_adj_',I0.0,'_',I0.0,'.dat')") M_in_Adj%dim(1), M_in_Adj%dim(2)
      CALL writeMatFile(M_in_Adj, unit=10, file=str, format ="(ES12.5,SP,ES13.5,' i   ')")



      ! -----------------------------------------
      ! traces computation
      ! -----------------------------------------
      trace = .TR.(M_in)
      trace_adj = .TR.(M_in_Adj)
      PRINT "('Trace input matrix =', ES12.5,SP,ES13.5,' i   ')", trace
      PRINT "('Trace adjoint matrix =', ES12.5,SP,ES13.5,' i   ')", trace_adj



      ! -----------------------------------------
      ! final deallocation of matrices
      ! -----------------------------------------      
      DEALLOCATE(M_in%elem)
      CALL checkpoint(debug, str="Deallocated input matrix")
      DEALLOCATE(M_in_Adj%elem)
      CALL checkpoint(debug, str="Deallocated adjoint matrix")

      
      STOP
END PROGRAM test_DCmatrix