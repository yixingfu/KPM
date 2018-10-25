! Adapted From MKL GUIDE=Tue 23 Oct 2018 8:07:36 PM DST
! Last Modified=Tue 23 Oct 2018 10:25:42 PM DST
! This is a subroutine that does complex to complex FFT
        subroutine FFT_3D(Nx, Ny, Nz, Ntot, X_in, X_out)
            integer,intent(in)::Nx,Ny,Nz,Ntot
            complex*16,dimension(Ntot),intent(in)::X_in
            complex*16,dimension(Ntot),intent(out)::X_out
            type(DFTI_DESCRIPTOR), POINTER :: My_Desc_Handle
            Integer :: Status
            Integer, dimension(3)::L
            
            L(1) = Nx
            L(2) = Ny
            L(3) = Nz

            X_out = X_in 
            ! Perform a complex to complex transform
            Status = DftiCreateDescriptor( My_Desc_Handle, DFTI_DOUBLE,&
                DFTI_COMPLEX, 3, L)
            Status = DftiCommitDescriptor( My_Desc_Handle)
            Status = DftiComputeForward( My_Desc_Handle, X_out)
            Status = DftiFreeDescriptor(My_Desc_Handle)
        end subroutine FFT_3D


