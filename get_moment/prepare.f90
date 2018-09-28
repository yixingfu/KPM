! Created =Tue 12 Dec 2017 03:10:19 PM STD
! Last Modified=Thu 20 Sep 2018 05:20:50 PM DST
! 

      ! This file prepares a few derived parameters from input file
      BHZ = .false.
      PIFLUX =  .false.
      GRAPHENE = .false.
      LRH1D = .false.
      if (MODEL_TYPE.eq.TYPE_BHZ) then
          BHZ=.true.
          PerSite=1+4*D
          N = 2*(L**D)
          ! Potential + spin up/down fwd/bwd 
      else if (MODEL_TYPE.eq.TYPE_PI_SM) then
          PIFLUX=.true.
          PerSite=1+2*D
          N = (L**D)
          ! Potential +  fwd/bwd no spin
          ! This is the version w/o t'
      else if (MODEL_TYPE.eq.TYPE_GRAPHENE) then
          GRAPHENE=.true.
          PerSite=1+3
          N = 2*(L**D)
      else if (MODEL_TYPE.eq.TYPE_LRH1D) then
          LRH1D=.true.
          PerSite=1+2 ! back and forward, momentum space
          N = L! no spin
          if (D.ne.1) then
              write(*,*) "Error: D=1 only for LRH model"
          endif
      else
          PerSite=1+2*D
          N = 2*(L**D)
          ! Potential + spin up/down fwd/bwd matching
      endif

      NNZ = PerSite*N ! fwd & bwd each site per dim + disorder
        if (GRAPHENE) then
            JNNZ=3*N
        else if (LRH1D) then
            JNNZ=N
            write(*,*) "not implemented - LRH 1D cond"
        else if (BHZ) then
            JNNZ = 4*N
            write(*,*)"J has 4 nonzero elements for each  point"
        else 
      JNNZ = (2)*N ! fwd & bwd each site @ x
        endif
!      write(*,*)D,N,L,NNZ


        XYNNZ = N

      EigValTot = 0
      EigValLancTot = 0


