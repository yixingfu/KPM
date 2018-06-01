! Created =Tue 12 Dec 2017 03:10:19 PM STD
! Last Modified=Fri 01 Jun 2018 02:59:33 PM DST
! 

      ! This file prepares a few derived parameters from input file
      BHZ = .false.
      PIFLUX =  .false.
      GRAPHENE=.false.
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
      else
          PerSite=1+2*D
          N = 2*(L**D)
          ! Potential + spin up/down fwd/bwd matching
      endif
      
      NNZ = PerSite*N ! fwd & bwd each site per dim + disorder

      JNNZ = (2)*N ! fwd & bwd each site @ x
!      write(*,*)D,N,L,NNZ


      EigValTot = 0
      EigValLancTot = 0


