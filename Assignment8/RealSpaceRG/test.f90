! Program to test...
! 
! Example of execution:
! ./test_RG ....
! will ...
!
! Author: Paolo Zinesi
!

! main program
PROGRAM test_RG
      USE checkpoint_mod
      USE ManyBodyUtils_mod
      IMPLICIT NONE
      !LOGICAL :: debug = .TRUE.

      INTEGER, PARAMETER :: NN = 3
      REAL :: ti, tf

      REAL*8, DIMENSION(:,:), ALLOCATABLE :: HN
      REAL*8, DIMENSION(2,2) :: sigmaZ
      REAL*8, DIMENSION(4,4) :: sigmaXsigmaX
      ALLOCATE(HN(2**NN,2**NN))

      ! sigmax sigmax
      sigmaXsigmaX = 0.D0
      CALL fillSigmaXSigmaX(sigmaXsigmaX, 1, 2)

      ! sigmaz
      sigmaZ = 0.D0
      CALL fillSigmaZ(sigmaZ, 1, 1)


      ! fill all sigmax sigmax interactions
      HN = 0.D0
      CALL CPU_TIME(ti)
      CALL generalTensorProduct(HN, sigmaXsigmaX, sigmaZ, 2, 1)
      CALL CPU_TIME(tf)
      PRINT "('XXZ computation with generalTensorProduct = ', ES15.7)", tf-ti


      CALL writeMatFile(HN, unit=10, file="HN.dat", format="(ES24.16)")

      DEALLOCATE(HN)


END PROGRAM test_RG