! Branched=Thu 31 May 2018 07:42:14 PM DST
! Last Modified=Wed 06 Jun 2018 05:45:04 PM DST
      !This file creates H for graphene
      ! 
      !The matrix is stored as CSR(A,col,rp)
      ! See Fradkins 16.3.2
      ! Set t = 1/2
      ! Quasi periodic parameter W
      ! 2D only
      ! Twist not implemented

      allocate(A(NNZ),col(NNZ),rp(N+1))
      ! NO TWIST
!      if (fixedTwist) then
!          Twist = OrigTwist*pi/real(L)
!      else
!          allocate(TwistAll(num_procs*3*seq_rep))
!          call random_number(TwistAll)
!          Twist=TwistAll(my_id*3+1:my_id*3+3)
!          write(*,*) 'twist:',Twist
!          deallocate(TwistAll)
!          Twist = (Twist)*pi/L! 0 to pi??
!      endif
      Twist =  0

      t = 0.5d0 ! the only term




      if (D.eq.2) then
!          allocate(phase_all(D,L*L))
          include "2d_QP_prepare_GRAPHENE.f90" 

          ! reset before get going
          rp = 0
          col = 0
          col_ind = 1
          rp_ind = 1
          eps_ind = 1
          do j=1,L! y
          do i=1,L! x
          do s=0,1
          ! set row pointer and disorder term (almost same for every
          ! model)
          rp(rp_ind) = PerSite*rp_ind-(PerSite-1)
          col(col_ind) = rp_ind
          rp_ind = rp_ind+1
          A(col_ind) = eps(eps_ind)
          col_ind = col_ind+1

          xx=i*ix+j*jx+s*AB
          yy=i*iy+j*jy
          
          do NNi=1,3
          
          ! NN: c^\dagger(row), c(col). col get varied
          i_ = safe_mod(i+NNx(s,NNi),L)
          j_ = safe_mod(j+NNy(s,NNi),L)
!          write(*,*)i,j,i_,j_
          xx_ = i_*ix+j_*jx+(1-s)*AB ! A to B, B to A
          yy_ = i_*iy+j_*jy
          col(col_ind) = xys2i(i_,j_,1-s,L)

          ! QP when Trnd or TQP turns on, it modifies the 
          ! real part of t
          ! take midway between itself and its neighbor
          t_tmp = t + Trnd*random2D((xx+xx_)/2d0,(yy+yy_)/2d0,P,Q)&
              + TQP*quasiperiodic((xx+xx_)/2d0,(yy+yy_)/2d0,P,Q,phase)

          A(col_ind) = t_tmp!*open_bc(i,i_,L,OPEN_BC_x)
          col_ind = col_ind+1
          enddo

          eps_ind = eps_ind+1
          end do


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


