! Created=Wed 13 Dec 2017 01:15:17 PM STD
! Last Modified=Fri 21 Sep 2018 11:42:45 PM DST
      ! This file computes the Jxy moment
      ! for arbitrary direction
      ! dir_a and dir_b
      !  we look for \sigma_ab
      ! calculating \mu_{m,n}

      if (Dir_a .eq.DIR_X) then
          RaA = XA
      else if (Dir_a .eq. DIR_Y) then
          RaA = YA
      endif

      if (Dir_b .eq. DIR_X) then
          RbA = XA
      else if (Dir_b .eq. DIR_Y) then
          RbA = YA
      endif

      

      allocate(mu2d_tot(0:Nc-1,0:Nc-1),mu2d2_tot(0:Nc-1,0:Nc-1))
      allocate(mu2d_avg(0:Nc-1,0:Nc-1),mu2d2_avg(0:Nc-1,0:Nc-1))
      allocate(mu2d(0:Nc-1,0:Nc-1))
      allocate(psi0R(N),psi0(N),psi1(N),psi_tmp(N))
      allocate(psi0_out(N),psi1_out(N),psi_tmp_out(N))
      allocate(psi_in(N),psi_p_in(N),psi_pp_in(N))
      allocate(psi_out(N),psi_p_out(N),psi_pp_out(N))
      allocate(psi_all_out(0:Nc-1,N))
! those marked with in, multiplies J at last
      !those marked with out, multiplies J at first
!<out|in>=<0|JaHHHHJbHHHHH|0>

      mu2d_tot = 0
      mu2d2_tot = 0


      ! Repetition for KPM
      do i=1,Rep
      mu2d = 0

      ! set \psi(0)
      call ResetRandSeed()
      call random_number(psi0R)
      psi0 = zexp(dcmplx(0d0,2.0d0*pi*psi0R))

      ! \psi(0)_out = Ja\psi(0)
!      call CSRmultVc16(N,NNZ,JA,Jrp,Jcol,psi0,psi0_out)
      call op_commutator(N,&
          XYNNZ,RaA,XYrp,XYcol,&
          NNZ,  A,  rp,  col,&
          psi0,JaPsi)
      TnPsi = psi0
      TmJaPsi = JaPsi ! m=0
      
      !  Ja Tm(H) Jb\psi(0), m=0
      call op_commutator(N,&
          XYNNZ,RbA,XYrp,XYcol,&
          NNZ,  A,  rp,  col,&
          TmJaPsi,JbTmJaPsi)

       
      ! \mu_xx_{0,0} = \psi(0)_out \cdot \psi(0)_in
      psi_all_out(0,:) = JbTmJaPsi
      mu2d(0,0) = - dot_product(JbTmJaPsi,TnPsi)

      ! m = m+1
      TmpJaPsi = TmJaPsi
      call CSRmultVc16(N,NNZ,A,rp,col, &
          TmpJaPsi,TmJaPsi)
      call op_commutator(N,&
          XYNNZ,RbA,XYrp,XYcol,&
          NNZ,  A,  rp,  col,&
          TmJaPsi,JbTmJaPsi)
! Let's set Tmpp to be 0
      ! n = n+1
      TnpPsi = TnPsi
      call CSRmultVc16(N,NNZ,A,rp,col, &
          TnpPsi,TnPsi)

       
      ! \psi(1) = H\psi(0)
      call CSRmultVc16(N,NNZ,A,rp,col,psi0,psi1)
      call CSRmultVc16(N,JNNZ,JA,Jrp,Jcol,psi1,psi_tmp)
      ! ################
      ! CONSTRUCTION ZONE
      ! ################
      call 
      mu2d(0,1) = dot_product(psi0_out,psi_tmp)
      mu2d(1,0) = mu2d(0,1)

      call CSRmultVc16(N,NNZ,A,rp,col,psi0_out,psi1_out)
      mu2d(1,1) = dot_product(psi1_out,psi_tmp)

      psi_pp_in = psi0
      psi_p_in  = psi1
      ! now the central part: do the loops and find mu2d
      ! two ways: 1. slow; 2. fast
      ! what we need is mu2d
        
          ! here is the fast way (takes more memory)
          ! first: construct all psi and save
          ! then: find the expectation one by one.
          psi_pp_out = psi0_out
          psi_p_out = psi1_out
          psi_all_out(0,:) = psi_pp_out
          psi_all_out(1,:) = psi_p_out


          do j=2,Nc-1
          call CSRmultVc16(N,NNZ,A,rp,col,psi_p_in,psi_tmp)
          psi_in = 2d0*psi_tmp-psi_pp_in
          call CSRmultVc16(N,JNNZ,JA,Jrp,Jcol,psi_in,psi_tmp)
          psi_pp_in = psi_p_in
          psi_p_in  = psi_in

          call CSRmultVc16(N,NNZ,A,rp,col,psi_p_out,psi_tmp_out)
          psi_out = 2d0*psi_tmp_out-psi_pp_out
          psi_pp_out = psi_p_out
          psi_p_out  = psi_out
          psi_all_out(j,:) = psi_out

          do k=0,j
          mu2d(j,k) = dot_product(psi_all_out(k,:),psi_tmp)
          mu2d(k,j) = mu2d(j,k)
          End do
          enddo





      mu2d_tot = mu2d_tot + mu2d
      mu2d2_tot = mu2d2_tot + mu2d*mu2d
      enddo
      mu2d_avg = mu2d_tot/(Rep*N)
      mu2d2_avg = mu2d2_tot/(Rep*N)
      deallocate(mu2d_tot,mu2d2_tot,mu2d)
      deallocate(psi0R,psi0,psi1,psi_tmp)
      deallocate(psi0_out,psi1_out,psi_tmp_out)
      deallocate(psi_in,psi_p_in,psi_pp_in)
      deallocate(psi_out,psi_p_out,psi_pp_out)
      deallocate(psi_all_out)
