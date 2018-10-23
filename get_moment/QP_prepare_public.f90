! The part of preparing QP parameters, the
! same for all models


      phase = 0
      if (DISORDER_TYPE.eq.DISORDER_QP) then
          WQP = W
          Wrnd = 0d0
      else if (DISORDER_TYPE.eq.DISORDER_RAND) then
          WQP = 0d0
          Wrnd = W
      else if (DISORDER_TYPE.eq.DISORDER_BOTH) then
          WQP = W
          Wrnd = W2
      endif
          if (RandPhase) then
              allocate(TwistAll(num_procs*3*seq_rep))
              call random_number(TwistAll)
              phase = TwistAll(my_id*3+1:my_id*3+D)
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

          ! set Q
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
          write(*,*) "---Q set to ",Q,",Qn=",Qn
        endif

          if (my_id.eq.0) then
              write(*,*)"Q=" , Q,"; M=",fiboM
          endif





      TQP = 0d0
      Trnd = 0d0
      t0 = 1d0
      if (RandType.eq.RandHOP) then
          WQP = 0d0
          Wrnd = 0d0
          t0 = dsqrt(1d0-W**2)
          if (DISORDER_TYPE.eq.DISORDER_QP) then
              TQP = W
              Trnd = 0d0
          else if (DISORDER_TYPE.eq.DISORDER_RAND) then
              TQP = 0d0
              Trnd = W
          else if (DISORDER_TYPE.eq.DISORDER_BOTH) then
              TQP = 0d0
              Trnd =  0d0
              write(*,*) "Hopping not implemented for QP+RAND"
          endif
      endif


