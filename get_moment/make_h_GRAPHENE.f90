! Branched=Thu 31 May 2018 07:42:14 PM DST
! Last Modified=Fri 27 Jul 2018 10:42:46 AM DST
      !This file creates H for graphene
      ! 
      !The matrix is stored as CSR(A,col,rp)
      ! See Fradkins 16.3.2
      ! Set t = 1/2
      ! Quasi periodic parameter W
      ! 2D only
      ! Twist testing
        ! Note: xx,yy,ix,jx,iy,jy not used 
        ! over HERE
        ! unless we want to do hopping model


      allocate(A(NNZ),col(NNZ),rp(N+1))

      include "2d_QP_prepare_GRAPHENE.f90" 

      ! Get Twist x,y
      if (fixedTwist) then
          Twist = OrigTwist*pi/real(L)
      else
          ! we need 3 phases, not 3 twists. 
          ! here twist has dimension 3 to keep up
          ! with other cases
          allocate(TwistAll(num_procs*3*seq_rep))
          call random_number(TwistAll)
          Twist=TwistAll(my_id*3+1:my_id*3+3)
          write(*,*) 'twist:',Twist
          deallocate(TwistAll)
          Twist = (Twist)*2d0*pi/L! 0 to pi??
      endif

      t = 0.5d0 ! the only term

      texp_theta = t
      ! A row B col

      ! NOTE: twist has nothing to do with
      ! the QP potential.
      ! This is multiplying the DUAL.
      ! trust it for now.
      ! IF we want to change lattice shape, 
      ! then these need to be modified. not for now
      texp_theta(0,1) = (2d0/3d0)*Twist(1)+(-1d0/3d0)*Twist(2)
      texp_theta(0,2) = texp_theta(0,1) - Twist(1)
      texp_theta(0,3) = texp_theta(0,2) + Twist(2)
      texp_theta(1,:) = -texp_theta(0,:)
      ! t exp(i *) 
      texp_theta = t*cdexp(III*texp_theta)
      write(*,*)texp_theta




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

      !xx=i*ix+j*jx+s*ABx
      !yy=i*iy+j*jy+s*ABy
      HC_r_ref = (L*HC_a1+L*HC_a2+HC_b)/2d0

      do NNi=1,3

      ! NN: c^\dagger(row), c(col). col get varied
      i_ = safe_mod(i+NNx(s,NNi),L)
      j_ = safe_mod(j+NNy(s,NNi),L)
!          write(*,*)i,j,i_,j_
      !xx_ = i_*ix+j_*jx+(1-s)*ABx ! A to B, B to A
      !yy_ = i_*iy+j_*jy+(1-s)*ABy
      col(col_ind) = xys2i(i_,j_,1-s,L)

      ! QP when Trnd or TQP turns on, it modifies the 
      ! real part of t
      ! take midway between itself and its neighbor
!          t_tmp = t 
      ! Hopping off for now
      !+ Trnd*random2D((xx+xx_)/2d0,(yy+yy_)/2d0,P,Q)&
      !  + TQP*quasiperiodic((xx+xx_)/2d0,(yy+yy_)/2d0,P,Q,phase)
      HC_r = i*HC_a1+j*HC_a2+s*HC_b
      HC_r_ = i_*HC_a1+j_*HC_a2+(1-s)*HC_b

      A(col_ind) = texp_theta(s,NNi) & 
          * set_shape(HC_r,HC_r_,HC_r_ref)
      col_ind = col_ind+1
      enddo

      eps_ind = eps_ind+1
      end do


      End do
      End do
      rp(rp_ind) = NNZ+1

      deallocate(eps)

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


