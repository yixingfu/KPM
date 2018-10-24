! Last Modified=
! Last Modified=

! 2D complex to complex
Complex ::  X_2D(32,100)
Complex ::  X(3200)
Equivalence (X_2D, X)
type(DFTI_DESCRIPTOR), POINTER :: My_Desc1_Handle, My_Desc2_Handle
Integer :: Status, L(2)
!...put input data into X_2D(j,k), Y_2D(j,k), 1<=j=32,1<=k<=100
L(1) = 32
L(2) = 100
!...the transform is a 32-by-100

X_2D = 2
 
! Perform a complex to complex transform
Status = DftiCreateDescriptor( My_Desc1_Handle, DFTI_SINGLE,&
          DFTI_COMPLEX, 2, L)
Status = DftiCommitDescriptor( My_Desc1_Handle)
Status = DftiComputeForward( My_Desc1_Handle, X)
Status = DftiFreeDescriptor(My_Desc1_Handle)
! result is given by X_2D(j,k), 1<=j<=32, 1<=k<=100
        open(15,file="test.txt")
        write(15,*)X_2D
        close(15)
        end program main
