! prepare for QP
! This generates correlated potential term
! eps has a term on each site
! For GRAPHENE
! Fully correlated version
! also set parameter for NNs

      phase = 0
      if (QP) then
          Wrnd = 0d0
          WQP = W
          if (RandPhase) then
              allocate(TwistAll(num_procs*3*seq_rep))
              call random_number(TwistAll)
              phase = TwistAll(my_id*3+1:my_id*3+3)
              deallocate(TwistAll)
          else 
              phase = inputPhase
          endif

          phase = phase*2d0*pi
          P = 2.0d0*pi*fibonacci(fiboN-fiboM,FiboBasis)&
                /fibonacci(fiboN,FiboBasis)
          Q = P
          if (commC .ne. 0) then
              ! if set commensurate C, then Q=P=2pi/C
              P=2.0d0*pi/commC
              Q=P
          endif

          ! offset Q
        if (setQ .ne. 0) then
          P = setQ
          Q = P
        write(*,*) "Q set to" , Q
        endif
        if (Qn .ne. 0) then
          P = 2d0*pi*Qn/L
          Q = P
          write(*,*) "---Q set to ",Q
        endif

          if (my_id.eq.0) then
              write(*,*)"Q=" , Q,"; M=",fiboM
          endif
      else 
          WQP = 0d0
          Wrnd = W
      End if
      TQP = 0d0
      Trnd = 0d0
      t0 = 1d0
      if (RandType.eq.RandHOP) then
          WQP = 0d0
          Wrnd = 0d0
          t0 = dsqrt(1d0-W**2)
          if (QP) then
              TQP = W
              Trnd = 0d0
          else 
              TQP = 0d0
              Trnd = W
          endif
      endif

      ! L*L lattice sites, 2 particles each
      ! parameters for honeycomb shape
     if (HONEYCOMB_BASIS.eq.HC_LAT) then
! Along lattice vectors
      jx=0d0
      jy=1d0
      ix=1d0
      iy=0d0
      ABy=-1d0/3d0
      ABx=2d0/3d0
! ~~~~~~~~~~~~~~~~
  else if (HONEYCOMB_BASIS.eq.HC_XY) then
! Along x,y
      jx = 0d0
      jy = dsqrt(3d0)
      ix = 1.5d0
      iy = dsqrt(3d0)/2d0
      ABx=1d0
      ABy=0d0
! ~~~~~~~~~~~~~~~~
  else if (HONEYCOMB_BASIS.eq.HC_set) then
!!!! EXPLANATION: We are setting matrix T
! such that (i, j) = T (delta x,y), where
! delta x,y represent the direction we are 
! setting QP. These directions are normalized
! so that it respect finite size periodicity.
! Each lattice point has its i,j,
! and we want to convert to x,y defined on
! those basis. If choosing this option,
! make sure all the 6 parameters are done
! correctly
      jx = HC_jx_in
      jy = HC_jy_in
      ix = HC_ix_in
      iy = HC_iy_in
        ABx = 2d0/3d0*ix - 1d0/3d0*jx
        ABy = 2d0/3d0*iy - 1d0/3d0*jy
  else if (HONEYCOMB_BASIS.eq.HC_set_theta) then
      ! This is giving a choice that is a rotation
      ! of angle theta, and conforming to the 
      ! requiement discussed above.

        ! counterclockwise theta
        a1 = 1d0
        b2 = 1d0
        a2 = -2d0*dtan(HC_theta_in)/(dsqrt(3d0)+dtan(HC_theta_in))
        b1 =  2d0*dtan(HC_theta_in)/(dsqrt(3d0)-dtan(HC_theta_in))
        HC_denom = a1*b2-a2*b1
        ix = b2/HC_denom
        iy = -b1/HC_denom
        jx = -a2/HC_denom
        jy = a1/HC_denom
        ABx = 2d0/3d0*ix - 1d0/3d0*jx
        ABy = 2d0/3d0*iy - 1d0/3d0*jy
  endif

! Along x,y, as parameters
      real_jx = 0d0
      real_jy = dsqrt(3d0)
      real_ix = 1.5d0
      real_iy = dsqrt(3d0)/2d0
      real_ABx=1d0
      real_ABy=0d0
! ~~~~~~~~~~~~~~~~

      allocate(eps(L*L*2)) 
      eps_ind = 1
      if (RandType.eq.RandPOT) then
          do j=1,L
          do i=1,L
          xx = i*ix+j*jx
          yy = i*iy+j*jy
          ! site A: + 0
          eps(eps_ind) = &
              WQP*quasiperiodic_HC(xx,yy,P,Q,phase)&
              + Wrnd*random2D(xx,yy,P,Q)
          eps_ind=eps_ind+1

          ! site B: + 1
          eps(eps_ind) = &
              WQP*quasiperiodic_HC(xx+ABx,yy+ABy,P,Q,phase)&
              + Wrnd*random2D(xx+ABx,yy+ABy,P,Q)
          eps_ind=eps_ind+1
          enddo
          enddo
          ! normalize
          eps = eps - sum(eps)/real(2*L*L)
      else if (RandType.eq.RandHOP) then
          eps = 0d0
      endif


      ! parameter
      ! 0 is A, left site.
      ! clockwise from same index A-B
      ! left bottom to right top
      NNx(0,1) = 0
      NNy(0,1) = 0
      NNx(0,2) = -1
      NNy(0,2) = 0
      NNx(0,3) = -1
      NNy(0,3) = 1

      NNx(1,1) = 0
      NNy(1,1) = 0
      NNx(1,2) = 1
      NNy(1,2) = 0
      NNx(1,3) = 1
      NNy(1,3) = -1

