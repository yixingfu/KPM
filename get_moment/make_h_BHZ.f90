! Created=Tue 12 Dec 2017 03:28:22 PM STD
! Last Modified=Fri 28 Sep 2018 01:36:16 AM DST
      !This file creates H
      !The matrix is stored as CSR(A,col,rp)
      !

      allocate(A(NNZ),col(NNZ),rp(N+1))
      if (fixedTwist) then
          Twist = OrigTwist*pi/real(L)
      else

          allocate(TwistAll(num_procs*3*seq_rep))
          call random_number(TwistAll)
          Twist=TwistAll(my_id*3+1:my_id*3+3)
          write(*,*) 'twist:',Twist
          deallocate(TwistAll)

          Twist = (Twist)*pi/L! 0 to pi??
      endif

      ! BHZ Definition
        ! f is for c(r+ex)\dagger c(r). col = row - 1
      if (BHZ) then
          xf=-BHZ_SPIN * III*pauli_x + pauli_z
          xb= BHZ_SPIN * III*pauli_x + pauli_z

          yf=-III*pauli_y + pauli_z
          yb= III*pauli_y + pauli_z

          zf=0d0
          zb=0d0

          zz=(BHZ_M-2d0)*pauli_z
      endif

      ! Add Twist
      txf = xf*cdexp(dcmplx(0d0,Twist(1)))
      txb = xb*cdexp(dcmplx(0d0,-Twist(1)))
      tyf = yf*cdexp(dcmplx(0d0,Twist(2)))
      tyb = yb*cdexp(dcmplx(0d0,-Twist(2)))
      tzf = zf*cdexp(dcmplx(0d0,Twist(3)))
      tzb = zb*cdexp(dcmplx(0d0,-Twist(3)))
      if (my_id .eq. 0) then
          write(*,*) 'BHZ:',BHZ
          write(*,*) 'BHZ spin:',BHZ_SPIN
          write(*,*) 'BHZ M:',BHZ_M
          write(*,*) '-',txf
          write(*,*) '-',txb
          write(*,*) '-',tyf
          write(*,*) '-',tyb
          write(*,*) '-',tzf
          write(*,*) '-',tzb
      endif





      if (D.eq.2) then
        include "2d_QP_prepare_BHZ.f90"

          rp = 0
          col = 0
          col_ind = 1
          rp_ind = 1
          eps_ind = 1

          do j=1,L! y
          do i=1,L! x
          do s = 0,1
          s_ = 1-s
          ! set row pointer and disorder term
          rp(rp_ind) = PerSite*rp_ind-(PerSite-1)
          col(col_ind) = rp_ind
          rp_ind = rp_ind+1
          A(col_ind) = eps(eps_ind) + zz(s,s)! zz for BHZ
          col_ind = col_ind+1

          ! x forward
          t_tmp = t0 + Trnd*random2D(i-0.5d0,j+0d0,P,Q)&
              + TQP*quasiperiodic(i-0.5d0,j+0d0,P,Q,phase)
          i_=modulo(i-2,L)+1
          ind_r = xys2i(i_,j,s_,L)
          col(col_ind) = ind_r
          A(col_ind) = txf(s,s_)*t_tmp*open_bc(i,i_,L,OPEN_BC_x)
!        write(*,*)A(col_ind)
          col_ind = col_ind+1
          ind_r = xys2i(i_,j,s,L)
          col(col_ind) = ind_r
          A(col_ind) = txf(s,s)*t_tmp*open_bc(i,i_,L,OPEN_BC_x)
!        write(*,*)A(col_ind)
          col_ind = col_ind+1

          ! x backward
          t_tmp = t0 + Trnd*random2D(i+0.5d0,j+0d0,P,Q)&
              + TQP*quasiperiodic(i+0.5d0,j+0d0,P,Q,phase)
          i_=modulo(i,L)+1
          ind_r = xys2i(i_,j,s_,L)
          col(col_ind) = ind_r
          A(col_ind) = txb(s,s_)*t_tmp*open_bc(i,i_,L,OPEN_BC_x)
!        write(*,*)A(col_ind)
          col_ind = col_ind+1
          ind_r = xys2i(i_,j,s,L)
          col(col_ind) = ind_r
          A(col_ind) = txb(s,s)*t_tmp*open_bc(i,i_,L,OPEN_BC_x)
!        write(*,*)A(col_ind)
          col_ind = col_ind+1

          ! y forward
          t_tmp = t0 + Trnd*random2D(i+0.0d0,j-0.5d0,P,Q)&
              + TQP*quasiperiodic(i+0.0d0,j-0.5d0,P,Q,phase)
          j_=modulo(j-2,L)+1
          ind_r = xys2i(i,j_,s_,L)
          col(col_ind) = ind_r
          A(col_ind) = tyf(s,s_)*t_tmp*open_bc(j,j_,L,OPEN_BC_y)
!        write(*,*)A(col_ind)
          col_ind = col_ind+1
          ind_r = xys2i(i,j_,s,L)
          col(col_ind) = ind_r
          A(col_ind) = tyf(s,s)*t_tmp*open_bc(j,j_,L,OPEN_BC_y)
!        write(*,*)A(col_ind)
          col_ind = col_ind+1

          ! y backward
          t_tmp = t0 + Trnd*random2D(i+0.0d0,j+0.5d0,P,Q)&
              + TQP*quasiperiodic(i+0.0d0,j+0.5d0,P,Q,phase)
          j_=modulo(j,L)+1
          ind_r = xys2i(i,j_,s_,L)
          col(col_ind) = ind_r
          A(col_ind) = tyb(s,s_)*t_tmp*open_bc(j,j_,L,OPEN_BC_y)
!        write(*,*)A(col_ind)
          col_ind = col_ind+1
          ind_r = xys2i(i,j_,s,L)
          col(col_ind) = ind_r
          A(col_ind) = tyb(s,s)*t_tmp*open_bc(j,j_,L,OPEN_BC_y)
!        write(*,*)A(col_ind)
          col_ind = col_ind+1

          End do
          eps_ind = eps_ind+1
          End do
          End do
          rp(rp_ind) = NNZ+1

          deallocate(eps)
      else if (D.eq.3) then
          write(*,*) "BHZ for 3D not implemented!"
          ! not differentiating QP or not
          rp = 0
          col = 0
          col_ind = 1
          rp_ind = 1
          allocate(eps(L*L*L))
          eps = 0
!                call ResetRandSeed(my_id*7)
!                call random_number(eps)
!                eps = W*2.0d0*(eps-0.5d0)
!                idum=-my_id*17
          call random3D(eps,W,my_id)
!                write(*,*)eps
          eps_ind = 1
          do k=1,L!z
          do j=1,L!y 
          do i=1,L!x
          do s = 0,1
          s_ = 1-s
          ! set row pointer and disorder term
          rp(rp_ind) = 7*rp_ind-6
          col(col_ind) = rp_ind
          rp_ind = rp_ind+1
          A(col_ind) = eps(eps_ind)
          col_ind = col_ind+1

          ! x forward
          i_ = modulo(i,L)+1
          ind_r = xyzs2i(i_,j,k,s_,L)
          col(col_ind) = ind_r
          A(col_ind) = txf(s,s_)
          col_ind = col_ind+1

          ! x backward
          i_ = modulo(i-2,L)+1
          ind_r = xyzs2i(i_,j,k,s_,L)
          col(col_ind) = ind_r
          A(col_ind) = txb(s,s_)
          col_ind = col_ind+1

          ! y forward
          j_ = modulo(j,L)+1
          ind_r = xyzs2i(i,j_,k,s_,L)
          col(col_ind) = ind_r
          A(col_ind) = tyf(s,s_)
          col_ind = col_ind+1

          ! y backward
          j_ = modulo(j-2,L)+1
          ind_r = xyzs2i(i,j_,k,s_,L)
          col(col_ind) = ind_r
          A(col_ind) = tyb(s,s_)
          col_ind = col_ind+1

          ! caution! z has s->s, not s->s_
          ! z forward
          k_ = modulo(k,L)+1
          ind_r = xyzs2i(i,j,k_,s,L)
          col(col_ind) = ind_r
          A(col_ind) = tzf(s,s_)
          col_ind = col_ind+1

          ! z backward
          k_ = modulo(k-2,L)+1
          ind_r = xyzs2i(i,j,k_,s,L)
          col(col_ind) = ind_r
          A(col_ind) = tzb(s,s_)
          col_ind = col_ind+1

          End do
          eps_ind = eps_ind + 1
          End do
          End do
          End do
          rp(rp_ind) = NNZ+1
          deallocate(eps)

      endif
!      if (my_id.eq.0) then
!          open(111,file="MATRIX.txt",status="replace",&
!              form="unformatted",access="stream")
!          write(111) N
!          write(111) NNZ
!          write(*,*)"NNZ",NNZ
!          write(111) A, col, rp
!          close(111)
!      endif



! -------------------------EIGENVALUE
      include "get_eigenvalue.f90"
! --------------------EIGEnVALUE end

      ! normalization same for both

      call LanczosBound(N,NNZ,A,rp,col,1000,Emax,Emin)
      norm_a = (Emax-Emin)/(2d0-0.2d0)
      norm_b = (Emax+Emin)/2
        if (Set_norm_a .ne. 0) then
                norm_a = Set_norm_a
                norm_b = 0d0
        endif
      call rescale_to_1(N,NNZ,A,rp,col,norm_a,norm_b)


