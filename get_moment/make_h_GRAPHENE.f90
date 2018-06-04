! Branched=Thu 31 May 2018 07:42:14 PM DST
! Last Modified=Thu 31 May 2018 11:02:36 PM DST
      !This file creates H for graphene
      ! 
      !The matrix is stored as CSR(A,col,rp)
      ! See Fradkins 16.3.2
      ! Set t = 1/2
      ! Quasi periodic parameter W
      ! 2D only

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

      ! term between c^\dagger and c
      ! Twist
      tx = 0.5d0*cdexp(dcmplx(0d0,Twist(1)))
      tx_ = 0.5d0*cdexp(dcmplx(0d0,-Twist(1)))
      ty = 0.5d0*cdexp(dcmplx(0d0,Twist(2)))
      ty_ = 0.5d0*cdexp(dcmplx(0d0,-Twist(2)))
      

      if (my_id .eq. 0) then
          write(*,*)'tx - ',tx
          write(*,*)'tx_ - ',tx_
          write(*,*)'ty - ',ty
          write(*,*)'ty_ - ',ty_
      endif






      if (D.eq.2) then
!          allocate(phase_all(D,L*L))
          include "2d_QP_prepare.f90"
          ! reset before get going
          rp = 0
          col = 0
          col_ind = 1
          rp_ind = 1
          eps_ind = 1
          do j=1,L! y
          do i=1,L! x
          ! set row pointer and disorder term (almost same for every
          ! model)
          rp(rp_ind) = PerSite*rp_ind-(PerSite-1)
          col(col_ind) = rp_ind
          rp_ind = rp_ind+1
          A(col_ind) = eps(eps_ind)
          col_ind = col_ind+1

          ! x forward: c^\dagger(row), c(col). col get varied
          ! QP when Trnd or TQP turns on, it modifies the 
          ! real part of t
          t_tmp = tx + Trnd*random2D(i+0.5d0,j+0d0,P,Q)&
              + TQP*quasiperiodic(i+0.5d0,j+0d0,P,Q,phase)
          i_=modulo(i,L)+1
          ind_r = xy2i(i_,j,L)
          col(col_ind) = ind_r
          A(col_ind) = eiAx(modulo(i,2))*t_tmp*open_bc(i,i_,L,OPEN_BC_x)
          col_ind = col_ind+1

          ! x backward
          t_tmp = tx_ + Trnd*random2D(i-0.5d0,j+0d0,P,Q)&
              + TQP*quasiperiodic(i-0.5d0,j+0d0,P,Q,phase)
          i_=modulo(i-2,L)+1
          ind_r = xy2i(i_,j,L)
          col(col_ind) = ind_r
          A(col_ind) = eiAx(modulo(i_,2))*t_tmp*open_bc(i,i_,L,OPEN_BC_x)
          col_ind = col_ind+1

          ! y forward
          t_tmp = ty + Trnd*random2D(i+0.0d0,j+0.5d0,P,Q)&
              + TQP*quasiperiodic(i+0.0d0,j+0.5d0,P,Q,phase)
          j_=modulo(j,L)+1
          ind_r = xy2i(i,j_,L)
          col(col_ind) = ind_r
          A(col_ind) = eiAy(modulo(i,2))*t_tmp*open_bc(j,j_,L,OPEN_BC_y)
          col_ind = col_ind+1

          ! y backward
          t_tmp = ty_ + Trnd*random2D(i+0.0d0,j-0.5d0,P,Q)&
              + TQP*quasiperiodic(i+0.0d0,j-0.5d0,P,Q,phase)
          j_=modulo(j-2,L)+1
          ind_r = xy2i(i,j_,L)
          col(col_ind) = ind_r
          A(col_ind) = eiAy(modulo(i,2))*t_tmp*open_bc(j,j_,L,OPEN_BC_y)
          col_ind = col_ind+1

          eps_ind = eps_ind+1
          End do
          End do
          rp(rp_ind) = NNZ+1

          deallocate(eps)
      else if (D.eq.3) then
          write(*,*) "NO 3D for PI FLUX"
      endif

! FOR TESTING: SAVE MATRIX
!      if (my_id.eq.0) then
!          open(111,file="MATRIX.txt",status="replace",&
!              form="unformatted",access="stream")
!          write(111) N
!          write(111) NNZ
!          write(*,*)"NNZ",NNZ
!          write(111) A, col, rp
!          close(111)
!      endif
!


! -------------------------EIGENVALUE
      include "get_eigenvalue.f90"
! --------------------EIGEnVALUE end

      ! normalization same for both

      call LanczosBound(N,NNZ,A,rp,col,1000,Emax,Emin)
      norm_a = (Emax-Emin)/(2d0-0.2d0)
      norm_b = (Emax+Emin)/2
      call rescale_to_1(N,NNZ,A,rp,col,norm_a,norm_b)


