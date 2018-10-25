! Branched=Wed 05 Sep 2018 03:15:40 PM DST
! Last Modified=Wed 05 Sep 2018 03:45:49 PM DST
      !This file creates H for Long Range Hopping 1D
      !The matrix is stored as CSR(A,col,rp)
      !Momentum space 
      ! clean: H(k) = sgn(cos(k))|cos(k)|^\sigma
      allocate(A(NNZ),col(NNZ),rp(N+1))
      if (fixedTwist) then
          Twist = OrigTwist*pi/real(L)
      else
          allocate(TwistAll(num_procs*3*seq_rep))
          call random_number(TwistAll)
          Twist=TwistAll(my_id*3+1:my_id*3+3)
          twist=twist*2d0-1d0
          deallocate(TwistAll)
          Twist = (Twist)*pi/L! 0 to pi??
      endif


      ! phase
      phase = 0
      if (DISORDER_TYPE.eq.DISORDER_QP) then
          
          if (RandPhase) then
              allocate(TwistAll(num_procs*3*seq_rep))
              call random_number(TwistAll)
              phase = TwistAll(my_id*3+1:my_id*3+2)
              deallocate(TwistAll)
          else 
              phase = inputPhase
          endif
          phase = phase*2d0*pi
          P = 2.0d0*pi*fibonacci(fiboN-fiboM,FiboBasis)&
                /fibonacci(fiboN,FiboBasis)
          Q = P
          iQ = fibonacci(fiboN-fiboM,FiboBasis)

      else 
          write(*,*) "QP makes it sparse. Non QP not implemented"
      endif

          ! initiate indices
          rp = 0
          col = 0
          col_ind = 1
          rp_ind = 1

          U_fwd = 0.5d0*W*(dcos(phase(1))-III*dsin(phase(1)))
          U_bwd = 0.5d0*W*(dcos(phase(1))+III*dsin(phase(1)))


          do i=1,L! momentum 0 to 2pi for 1 to L
          ! set row pointer and self term
          rp(rp_ind) = PerSite*rp_ind-(PerSite-1)
          col(col_ind) = rp_ind
          k_i = (i-0.5d0)/L*2d0*pi+Twist(1)
          A(col_ind) = sign(abs(dcos(k_i))**LRH_sigma,dcos(k_i))
          rp_ind = rp_ind+1
          col_ind = col_ind+1

          !  forward
          i_=modulo(i+iQ-1,L)+1
          col(col_ind) = i_
          A(col_ind) = U_fwd
          col_ind = col_ind+1

          ! backward
          i_=modulo(i-iQ-1,L)+1
          col(col_ind) = i_
          A(col_ind) = U_bwd
          col_ind = col_ind+1

          End do
          rp(rp_ind) = NNZ+1

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
        if (Set_norm_a .ne. 0) then
                norm_a = Set_norm_a
                norm_b = 0d0
        endif
      call rescale_to_1(N,NNZ,A,rp,col,norm_a,norm_b)


