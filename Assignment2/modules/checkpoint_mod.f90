! Module to define checkpoints around the code. A boolean variable ("debug")
! allows to switch on the debugging mode in all the main program by changing it
! only once. Each checkpoint may print strings and/or variables.
! 
! Author: Paolo Zinesi
!

MODULE checkpoint_mod
      IMPLICIT NONE

      ! This interface defines the generic checkpoint subroutine.
      ! The structure is similar for all the specific subroutines,
      ! the only changes are in the type of numeric variable ("val") to print
      !
      ! inputs:
      ! - debug [logical]: Boolean value to enable/disable debugging
      ! - str [character(len=*), optional]: optional string to print
      ! - val [different types, optional]: optional numerical value to print
      !
      ! outputs: None
      ! 
      ! TODO: None
      ! 
      INTERFACE checkpoint
            MODULE PROCEDURE STRcheckpoint
            MODULE PROCEDURE Dcheckpoint, Acheckpoint
            MODULE PROCEDURE Icheckpoint, Jcheckpoint, Kcheckpoint
            MODULE PROCEDURE Ccheckpoint, CDcheckpoint
      END INTERFACE checkpoint

      CONTAINS

      ! Only print the (optional) string. This subroutine is called
      ! whenever no numerical value needs to be printed. All the other subroutines
      ! must then have "val" as input.
      SUBROUTINE STRcheckpoint(debug, str)
            IMPLICIT NONE

            LOGICAL :: debug
            CHARACTER (LEN=*), OPTIONAL :: str


            IF (debug .EQV. .TRUE.) THEN
                  IF(PRESENT(str)) THEN
                        ! PRINT "(A)", str
                        WRITE(*, "(A)", advance="yes") str
                  ENDIF
            ENDIF
      END SUBROUTINE STRcheckpoint

      ! Subroutine to print non-optional REAL*8 val
      SUBROUTINE Dcheckpoint(debug, str, val)
            IMPLICIT NONE

            LOGICAL :: debug
            CHARACTER (LEN=*), OPTIONAL :: str
            REAL*8 :: val


            IF (debug .EQV. .TRUE.) THEN
                  IF(PRESENT(str)) THEN
                        WRITE(*, "(A)", advance="no") str
                  ENDIF
                  PRINT *, val
            ENDIF
      END SUBROUTINE Dcheckpoint

      ! Subroutine to print non-optional REAL*4 val
      SUBROUTINE Acheckpoint(debug, str, val)
            IMPLICIT NONE

            LOGICAL :: debug
            CHARACTER (LEN=*), OPTIONAL :: str
            REAL*4 :: val


            IF (debug .EQV. .TRUE.) THEN
                  IF(PRESENT(str)) THEN
                        WRITE(*, "(A)", advance="no") str
                  ENDIF
                  PRINT *, val
            ENDIF
      END SUBROUTINE Acheckpoint

      ! Subroutine to print non-optional INTEGER*2 val
      SUBROUTINE Icheckpoint(debug, str, val)
            IMPLICIT NONE

            LOGICAL :: debug
            CHARACTER (LEN=*), OPTIONAL :: str
            INTEGER*2 :: val


            IF (debug .EQV. .TRUE.) THEN
                  IF(PRESENT(str)) THEN
                        WRITE(*, "(A)", advance="no") str
                  ENDIF
                  WRITE (*,"(I0.0)") val
            ENDIF
      END SUBROUTINE Icheckpoint

      ! Subroutine to print non-optional INTEGER*4 val
      SUBROUTINE Jcheckpoint(debug, str, val)
            IMPLICIT NONE

            LOGICAL :: debug
            CHARACTER (LEN=*), OPTIONAL :: str
            INTEGER*4 :: val


            IF (debug .EQV. .TRUE.) THEN
                  IF(PRESENT(str)) THEN
                        WRITE(*, "(A)", advance="no") str
                  ENDIF
                  WRITE (*,"(I0.0)") val
            ENDIF
      END SUBROUTINE Jcheckpoint

      ! Subroutine to print non-optional INTEGER*8 val
      SUBROUTINE Kcheckpoint(debug, str, val)
            IMPLICIT NONE

            LOGICAL :: debug
            CHARACTER (LEN=*), OPTIONAL :: str
            INTEGER*8 :: val


            IF (debug .EQV. .TRUE.) THEN
                  IF(PRESENT(str)) THEN
                        WRITE(*, "(A)", advance="no") str
                  ENDIF
                  WRITE (*,"(I0.0)") val
            ENDIF
      END SUBROUTINE Kcheckpoint

      ! Subroutine to print non-optional COMPLEX val
      SUBROUTINE Ccheckpoint(debug, str, val)
            IMPLICIT NONE

            LOGICAL :: debug
            CHARACTER (LEN=*), OPTIONAL :: str
            COMPLEX :: val


            IF (debug .EQV. .TRUE.) THEN
                  IF(PRESENT(str)) THEN
                        WRITE(*, "(A)", advance="no") str
                  ENDIF
                  PRINT *, val
            ENDIF
      END SUBROUTINE Ccheckpoint

      ! Subroutine to print non-optional DOUBLE COMPLEX val
      SUBROUTINE CDcheckpoint(debug, str, val)
            IMPLICIT NONE

            LOGICAL :: debug
            CHARACTER (LEN=*), OPTIONAL :: str
            DOUBLE COMPLEX :: val


            IF (debug .EQV. .TRUE.) THEN
                  IF(PRESENT(str)) THEN
                        WRITE(*, "(A)", advance="no") str
                  ENDIF
                  PRINT *, val
            ENDIF
      END SUBROUTINE CDcheckpoint

END MODULE checkpoint_mod