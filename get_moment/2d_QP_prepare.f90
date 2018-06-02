! prepare for QP
! This generates correlated potential term
! eps has a term on each site
! Fully correlated version

      phase = 0
      if (QP) then
          Wrnd = 0d0
          WQP = W
          if (RandPhase) then
              call random_number(phase)
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


      allocate(eps(L*L)) 
      eps_ind = 1
      if (RandType.eq.RandPOT) then
          do j=1,L
          do i=1,L

          eps(eps_ind) = &
              WQP*quasiperiodic(i+0d0,j+0d0,P,Q,phase)&
              + Wrnd*random2D(i+0d0,j+0d0,P,Q)
          eps_ind=eps_ind+1
          enddo
          enddo
          eps = eps - sum(eps)/real(L*L)
      else if (RandType.eq.RandHOP) then
          eps = 0d0
      endif
      phase_all(1,:) = phase(1)
      phase_all(2,:) = phase(2)


