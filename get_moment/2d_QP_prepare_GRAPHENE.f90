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
          R = P
          if (commC .ne. 0) then
              ! if set commensurate C, then Q=P=2pi/C
              P=2.0d0*pi/commC
              Q=P
              R=P
          endif

          ! offset Q
        if (setQ .ne. 0) then
          P = setQ
          Q = P
          R = P
        write(*,*) "Q set to" , Q
        endif
        if (Qn .ne. 0) then
          P = 2d0*pi*Qn/L
          Q = P
          R = P
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
      HC_a1(1) = 3d0/2d0 ! eastnorth
      HC_a1(2) = dsqrt(3d0)/2d0
      HC_a2(1) = 0d0 ! north
      HC_a2(2) = dsqrt(3d0)
      HC_b(1) = 1d0 ! east
      HC_b(2) = 0d0


     if (HONEYCOMB_BASIS.eq.HC_LAT) then
         ! Along lattice vectors
         call dual_2d(HC_a1,HC_a2,HC_d1,HC_d2)
         HC_d3 = - (HC_d1+HC_d2)
         ! ~~~~~~~~~~~~~~~~
     else if (HONEYCOMB_BASIS.eq.HC_XY) then
         ! Along x,y
         write(*,*)"THIS is doing something funny because of normalization"
         HC_d1 = HC_a1
         HC_d1(2) = 0d0
         HC_d2 = HC_a2
         HC_d2(1) = 0d0
         HC_d3 = 0d0
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
         !!!!!
         !!!!! SEE NOTES
        if (my_id .eq. 0) then
            write(*,*) "HC: delta = "
            write(*,*) HC_d1
            write(*,*) HC_d2
            write(*,*) HC_d3
        endif
  else if (HONEYCOMB_BASIS.eq.HC_set_theta) then
      ! This is giving a choice that is a rotation
      ! of angle theta, and conforming to the 
      ! requiement discussed above.
        write(*,*) "HC w/ theta: THIS IS NOT WORKING"

        ! counterclockwise theta
!        a1 = 1d0
!        b2 = 1d0
!        a2 = -2d0*dtan(HC_theta_in)/(dsqrt(3d0)+dtan(HC_theta_in))
!        b1 =  2d0*dtan(HC_theta_in)/(dsqrt(3d0)-dtan(HC_theta_in))
!        HC_denom = a1*b2-a2*b1
!        ix = b2/HC_denom
!        iy = -b1/HC_denom
!        jx = -a2/HC_denom
!        jy = a1/HC_denom
!        ABx = 2d0/3d0*ix - 1d0/3d0*jx
!        ABy = 2d0/3d0*iy - 1d0/3d0*jy
        HC_d1 = 0
        HC_d2 = 0
        HC_d3 = 0
  endif

! Along x,y, as parameters
!      real_jx = 0d0
!      real_jy = dsqrt(3d0)
!      real_ix = 1.5d0
!      real_iy = dsqrt(3d0)/2d0
!      real_ABx=1d0
!      real_ABy=0d0
! ~~~~~~~~~~~~~~~~

      allocate(eps(L*L*2)) 
      eps_ind = 1
      if (RandType.eq.RandPOT) then
          do j=1,L
          do i=1,L
          HC_r = i*HC_a1+j*HC_a2

          ! site A: + 0
          xx = dot_product(HC_r,HC_d1)
          yy = dot_product(HC_r,HC_d2)
          zz = dot_product(HC_r,HC_d3)
          eps(eps_ind) = &
              WQP*quasiperiodic_HC(xx,yy,zz,P,Q,R,phase)&
              + Wrnd*random2D(xx,yy,P,Q)
          eps_ind=eps_ind+1

          ! site B: + 1
          HC_r = HC_r + HC_b
          xx = dot_product(HC_r,HC_d1)
          yy = dot_product(HC_r,HC_d2)
          zz = dot_product(HC_r,HC_d3)
          eps(eps_ind) = &
              WQP*quasiperiodic_HC(xx,yy,zz,P,Q,R,phase)&
              + Wrnd*random2D(xx,yy,P,Q)
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
      ! they are acting on i,j
      ! and need to flip AB
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

