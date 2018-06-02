! prepare for QP
! This generates correlated potential term
! eps has a term on each site
! partially correlated version
      allocate(phase_vals(D,PIECE))
      ! setting parameters
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
      call random_number(phase_vals)
      phase_vals =  phase_vals * 2d0 * pi
      ! reset 
      phase_all = -1

      ! place phase_vals in random locations
      do i=1,PIECE
 20      call random_number(temp_i_real)
         temp_i = ceiling(temp_i_real*L*L)
         if (phase_all(1,temp_i).eq.-1) then
             phase_all(:,temp_i) = phase_vals(:,i)
         else
             goto 20
         endif
      enddo

      ! Now that the seeds are set, we grow them
      ! to get entire phase (phase_all)

      ! for each iteration, assign  a random starting point
      ! At each point, check if it is defined
      ! if defined, then extend it at a random direction
      ! When one location gets updated, also update the total location
      ! number count. End when all locations populated.
      set_count = 0
      total_count = L**D

      while (set_count .lt. total_count) do
        call random(temp_i_real)
        temp_start = ceiling(temp_i_real*total_count)
        do temp_i=1,total_count
          if (phase_all(1,i).eq.-1) CYCLE
          i = modulo(temp_i+temp_start,total_count) + 1
          call random_number(temp_i_real)
          if (temp_i_real .lt. 0.25) then
                i_ = i_left_2d(i,L)
          else if (temp_i_real.lt.0.5) then
                i_ = i_right_2d(i,L)
          else if (temp_i_real.lt.0.75) then
                i_ = i_up_2d(i,L)
          else
                i_ = i_down_2d(i,L)
          endif
          if (phase_all(1,i_).eq.-1) then
              phase_all(:,i_)=phase_all(:,i)
              set_count=set_count+1
          endif
        enddo
      enddo
      write(*,*)"phase for all"
      write(*,*)phase_all

      

      ! Populating eps
      allocate(eps(L*L))  !!2D
      eps_ind = 1
      if (RandType.eq.RandPOT) then
          do j=1,L
          do i=1,L

          eps(eps_ind) = &
              WQP*quasiperiodic(i+0d0,j+0d0,P,Q,phase_all(eps_ind))
          eps_ind=eps_ind+1
          enddo
          enddo
          eps = eps - sum(eps)/real(L*L)
      else if (RandType.eq.RandHOP) then
          eps = 0d0
      endif

        deallocate(phase_vals,phase_all)
