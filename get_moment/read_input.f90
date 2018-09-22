! Created=Tue 12 Dec 2017 03:19:48 PM STD
! Last Modified=Sat 22 Sep 2018 03:35:35 PM DST
      ! read inputfile name from command line
      call getarg(1,inputfile)
      ! read inputs D,L,Nc,W,QP?, from input file
      open(11, file=trim(inputfile))
      read(11,*) D,L,Nc,W,Rep,REALIZATION0
      read(11,*) QP,fixedTwist, slowOPTCOND
      read(11,*) task, RandType, Scramble
      read(11,*) OrigTwist
      read(11,*) outputfile
!!        read(11,*) Inherit, SaveAll
      read(11,*) ExactSpectrum, ExactStates
      read(11,*) RandPhase, inputPhase
      read(11,*) fiboM,commC,setQ,Qn
      read(11,*) MODEL_TYPE
      read(11,*) BHZ_SPIN,BHZ_M
      read(11,*) OPEN_BC_x,OPEN_BC_y,OPEN_BC_z
      read(11,*) LIMIT_CORRELATION,PIECE
      ! direction of velocity calculation
      read(11,*) HONEYCOMB_BASIS,HC_Jx,HC_Jy 
      ! direction of QP
      read(11,*) HC_d1, HC_d2, HC_d3, HC_d3_abs
      ! extra rotation -- not working for now!!!
      read(11,*) HC_theta_in
      read(11,*) FLAKE_SHAPE
      read(11,*) FiboBasis
      read(11,*) LRH_sigma
      read(11,*) CondTensor, Dir_a, Dir_b
      close(11)

      if (QP) then
          fiboN = inv_fibo(L,FiboBasis)
      endif

!         check
!         write(*,*) D,L,Nc,W,QP,outputfile



