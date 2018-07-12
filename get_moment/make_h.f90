! Created=Tue 12 Dec 2017 03:28:22 PM STD
! Last Modified=Sat 23 Jun 2018 02:46:03 PM EDT
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
          twist=twist*2d0-1d0
          write(*,*) 'twist:',Twist
          deallocate(TwistAll)

          Twist = (Twist)*pi/L! 0 to pi??
      endif



      ! Add Twist
      txf = xf*cdexp(dcmplx(0d0,Twist(1)))
      txb = xb*cdexp(dcmplx(0d0,-Twist(1)))
      tyf = yf*cdexp(dcmplx(0d0,Twist(2)))
      tyb = yb*cdexp(dcmplx(0d0,-Twist(2)))
      tzf = zf*cdexp(dcmplx(0d0,Twist(3)))
      tzb = zb*cdexp(dcmplx(0d0,-Twist(3)))
      if (my_id .eq. 0) then
          write(*,*) '-',txf
          write(*,*) '-',txb
          write(*,*) '-',tyf
          write(*,*) '-',tyb
          write(*,*) '-',tzf
          write(*,*) '-',tzb
      endif






      if (D.eq.2) then
          allocate(phase_all(D,L*L))

          if (LIMIT_CORRELATION) then
              include "2d_QP_prepare_piecewise.f90"
          else
              include "2d_QP_prepare.f90"
              phase_all(1,:) = phase(1)
              phase_all(2,:) = phase(2)
          endif
          ! initiate indices
          eps_ind = 1
          rp = 0
          col = 0
          col_ind = 1
          rp_ind = 1

          do j=1,L! y
          do i=1,L! x
          do s = 0,1
          s_ = 1-s
          ! set row pointer and disorder term
          rp(rp_ind) = PerSite*rp_ind-(PerSite-1)
          col(col_ind) = rp_ind
          rp_ind = rp_ind+1
          A(col_ind) = eps(eps_ind)
          col_ind = col_ind+1

          ! x forward
          i_=modulo(i,L)+1
          ind_r = xys2i(i_,j,s_,L)
          t_tmp = t0 + Trnd*random2D(i-0.5d0,j+0d0,P,Q)&
              + TQP*quasiperiodic(i-0.5d0,j+0d0,P,Q,&
              phase)
          col(col_ind) = ind_r
          A(col_ind) = txf(s,s_)*t_tmp*open_bc(i,i_,L,OPEN_BC_x)
          col_ind = col_ind+1

          ! x backward
          i_=modulo(i-2,L)+1
          ind_r = xys2i(i_,j,s_,L)
          t_tmp = t0 + Trnd*random2D(i+0.5d0,j+0d0,P,Q)&
              + TQP*quasiperiodic(i+0.5d0,j+0d0,P,Q,phase)
          col(col_ind) = ind_r
          A(col_ind) = txb(s,s_)*t_tmp*open_bc(i,i_,L,OPEN_BC_x)
          col_ind = col_ind+1

          ! y forward
          j_=modulo(j,L)+1
          ind_r = xys2i(i,j_,s_,L)
          t_tmp = t0 + Trnd*random2D(i+0d0,j-0.5d0,P,Q)&
              + TQP*quasiperiodic(i+0d0,j-0.5d0,P,Q,&
              phase)
          col(col_ind) = ind_r
          A(col_ind) = tyf(s,s_)*t_tmp*open_bc(j,j_,L,OPEN_BC_y)
          col_ind = col_ind+1

          ! y backward
          j_=modulo(j-2,L)+1
          ind_r = xys2i(i,j_,s_,L)
          t_tmp = t0 + Trnd*random2D(i+0d0,j+0.5d0,P,Q)&
              + TQP*quasiperiodic(i+0d0,j+0.5d0,P,Q,phase)
          col(col_ind) = ind_r
          A(col_ind) = tyb(s,s_)*t_tmp*open_bc(j,j_,L,OPEN_BC_y)
          col_ind = col_ind+1

          End do
          eps_ind = eps_ind+1
          End do
          End do
          rp(rp_ind) = NNZ+1

          deallocate(eps,phase_all)
      else if (D.eq.3) then
          ! not differentiating QP or not
          allocate(phase_all(D,L*L*L))
        write(*,*) '3d'
          if (LIMIT_CORRELATION) then
              include "3d_QP_prepare_piecewise.f90"
          else
              include "3d_QP_prepare.f90"
              phase_all(1,:) = phase(1)
              phase_all(2,:) = phase(2)
          endif
          ! initiate indices
          eps_ind = 1
          rp = 0
          col = 0
          col_ind = 1
          rp_ind = 1

          do k=1,L!z
          do j=1,L!y 
          do i=1,L!x
          do s = 0,1
          s_ = 1-s
          ! set row pointer and disorder term
          rp(rp_ind) = PerSite*rp_ind-(PerSite-1)
          col(col_ind) = rp_ind
          rp_ind = rp_ind+1
          A(col_ind) = eps(eps_ind)
          col_ind = col_ind+1

          ! x forward
          i_=modulo(i,L)+1
          ind_r = xyzs2i(i_,j,k,s_,L)
          col(col_ind) = ind_r
          t_tmp = t0 ! + HOPPING
          A(col_ind) = txf(s,s_)*t_tmp*open_bc(i,i_,L,OPEN_BC_x)
          col_ind = col_ind+1

          ! x backward
          i_=modulo(i-2,L)+1
          ind_r = xyzs2i(i_,j,k,s_,L)
          col(col_ind) = ind_r
          t_tmp = t0 ! + HOPPING
          A(col_ind) = txb(s,s_)*t_tmp*open_bc(i,i_,L,OPEN_BC_x)
          col_ind = col_ind+1

          ! y forward
          j_=modulo(j,L)+1
          ind_r = xyzs2i(i,j_,k,s_,L)
          t_tmp = t0 ! + HOPPING
          col(col_ind) = ind_r
          A(col_ind) = tyf(s,s_)*t_tmp*open_bc(j,j_,L,OPEN_BC_y)
          col_ind = col_ind+1

          ! y backward
          j_=modulo(j-2,L)+1
          ind_r = xyzs2i(i,j_,k,s_,L)
          t_tmp = t0 ! + HOPPING
          col(col_ind) = ind_r
          A(col_ind) = tyb(s,s_)*t_tmp*open_bc(j,j_,L,OPEN_BC_y)
          col_ind = col_ind+1

          ! caution! z has s->s, not s->s_
          ! z forward
          k_=modulo(k,L)+1
          ind_r = xyzs2i(i,j,k_,s,L)
          t_tmp = t0 ! + HOPPING
          col(col_ind) = ind_r
          A(col_ind) = tzf(s,s)*t_tmp*open_bc(k,k_,L,OPEN_BC_z)
          col_ind = col_ind+1

          ! z backward
          k_=modulo(k-2,L)+1
          ind_r = xyzs2i(i,j,k_,s,L)
          t_tmp = t0 ! + HOPPING
          col(col_ind) = ind_r
          A(col_ind) = tzb(s,s)*t_tmp*open_bc(k,k_,L,OPEN_BC_z)
          col_ind = col_ind+1

          End do
          eps_ind = eps_ind + 1
          End do
          End do
          End do
          rp(rp_ind) = NNZ+1
          deallocate(eps,phase_all)

      endif
      if (my_id.eq.0) then
          open(111,file="MATRIX.txt",status="replace",&
              form="unformatted",access="stream")
          write(111) N
          write(111) NNZ
          write(*,*)"NNZ",NNZ
          write(111) A, col, rp
          close(111)
      endif



! -------------------------EIGENVALUE
      include "get_eigenvalue.f90"
! --------------------EIGEnVALUE end

      ! normalization same for both

      call LanczosBound(N,NNZ,A,rp,col,1000,Emax,Emin)
      norm_a = (Emax-Emin)/(2d0-0.2d0)
      norm_b = (Emax+Emin)/2
      call rescale_to_1(N,NNZ,A,rp,col,norm_a,norm_b)


