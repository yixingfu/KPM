! Created=Tue 12 Dec 2017 03:19:48 PM STD
! Last Modified=Sat 27 Oct 2018 03:20:44 PM DST
      ! read inputfile name from command line
      call getarg(1,inputfile)
      ! read inputs D,L,Nc,W,QP?, from input file
      open(11, file=trim(inputfile))
      read(11,*) D,L,Nc,W,W2,Rep,REALIZATION0
        write(*,*)"1"
      read(11,*) DISORDER_TYPE,fixedTwist, slowOPTCOND
        write(*,*)"2"
      read(11,*) task, RandType, Scramble
        write(*,*)"3"
      read(11,*) OrigTwist
        write(*,*)"4"
      read(11,*) outputfile
        write(*,*)"5"
!!        read(11,*) Inherit, SaveAll
      read(11,*) ExactSpectrum, ExactStates, ExactIPR, EigvalCount,&
          BinCount,LSEmax, LanczosLS
        write(*,*)"6"
      read(11,*) RandPhase, inputPhase
        write(*,*)"7"
      read(11,*) fiboM,commC,setQ,Qn
        write(*,*)"8"
      read(11,*) MODEL_TYPE
        write(*,*)"9"
      read(11,*) BHZ_SPIN,BHZ_M
        write(*,*)"10"
      read(11,*) OPEN_BC_x,OPEN_BC_y,OPEN_BC_z
        write(*,*)"11"
      read(11,*) LIMIT_CORRELATION,PIECE
        write(*,*)"12"
      ! direction of velocity calculation
      read(11,*) HONEYCOMB_BASIS,HC_Jx,HC_Jy 
        write(*,*)"13"
      ! direction of QP
      read(11,*) HC_d1, HC_d2, HC_d3, HC_d3_abs
        write(*,*)"14"
      ! extra rotation -- not working for now!!!
      read(11,*) HC_theta_in
        write(*,*)"15"
      read(11,*) FLAKE_SHAPE
        write(*,*)"16"
      read(11,*) FiboBasis
        write(*,*)"17"
      read(11,*) LRH_sigma
        write(*,*)"18"
      read(11,*) CondTensor, Dir_a, Dir_b
        write(*,*)"19"
      read(11,*) Set_Norm_a
        write(*,*)"20"
        read(11,*) IQHE_B
        write(*,*) "21"
        read(11,*) SB
      close(11)
        allocate(LanczosEigs(seq_rep*8),MinEigs(seq_rep*15))
      if ((DISORDER_TYPE.eq.DISORDER_QP).or. &
                (DISORDER_TYPE.eq.DISORDER_BOTH)) then
          fiboN = inv_fibo(L,FiboBasis)
      endif

!         check
!         write(*,*) D,L,Nc,W,QP,outputfile



