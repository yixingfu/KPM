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
              phase = TwistAll(my_id*3+1:my_id*3+2)
              deallocate(TwistAll)
          else 
              phase = inputPhase
          endif

          phase = phase*2d0*pi
          P = 2.0d0*pi*fibonacci(fiboN-fiboM)/fibonacci(fiboN)
          Q = P
          if (commC .ne. 0) then
              ! if set commensurate C, then Q=P=2pi/C
              P=2.0d0*pi/commC
              Q=P
          endif

          ! offset Q
          P = P + setQ
          Q = P

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
      jx = 0d0
      jy = dsqrt(3d0)
      ix = 1.5d0
      iy = dsqrt(3d0)/2d0

      allocate(eps(L*L*2)) 
      eps_ind = 1
      if (RandType.eq.RandPOT) then
          do j=1,L
          do i=1,L
          xx = i*ix+j*jx
          yy = i*iy+j*jy
          ! site A: + 0
          eps(eps_ind) = &
              WQP*quasiperiodic(xx,yy,P,Q,phase)&
              + Wrnd*random2D(xx,yy,P,Q)
          eps_ind=eps_ind+1

          ! site B: + 1
          eps(eps_ind) = &
              WQP*quasiperiodic(xx+1d0,yy,P,Q,phase)&
              + Wrnd*random2D(xx+1d0,yy,P,Q)
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

