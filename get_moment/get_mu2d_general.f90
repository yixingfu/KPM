! Created=Wed 13 Dec 2017 01:15:17 PM STD
! Last Modified=Sat 22 Sep 2018 03:40:03 PM DST
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
      allocate(psi0R(N),psi0(N))
      allocate(JaPsi(N),JbTmJaPsi(N))
      allocate(TmJaPsi(N),TmpJaPsi(N),TmppJaPsi(N))
      allocate(TnPsi(N),TnpPsi(N),TnppPsi(N))
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
      TnPsi = psi0

      ! Ja\psi(0)
      call op_commutator(N,&
          XYNNZ,RaA,XYrp,XYcol,&
          NNZ,  A,  rp,  col,&
          psi0,JaPsi)
      TmJaPsi = JaPsi ! m=0
      
      !  Jb Tm(H) Ja\psi(0), m=0
      call op_commutator(N,&
          XYNNZ,RbA,XYrp,XYcol,&
          NNZ,  A,  rp,  col,&
          TmJaPsi,JbTmJaPsi)
       
      ! \mu_xx_{0,0} = \psi(0)_out \cdot \psi(0)_in
      psi_all_out(0,:) = JbTmJaPsi
      mu2d(0,0) = - dot_product(JbTmJaPsi,TnPsi)

      ! increment in m
      TmppJaPsi = 0d0
      TmpJaPsi = TmJaPsi
      do cond_m = 1,Nc-1
      call op_chebyshev(N,NNZ,A,rp,col,&
          TmppJaPsi,TmpJaPsi,TmJaPsi)
      TmppJaPsi = TmpJaPsi
      TmpJaPsi = TmJaPsi

      ! close and save
      call op_commutator(N,&
          XYNNZ,RbA,XYrp,XYcol,&
          NNZ,  A,  rp,  col,&
          TmJaPsi,JbTmJaPsi)
      psi_all_out(cond_m,:) = JbTmJaPsi
      enddo

      do cond_m = 0,Nc-1
        mu2d(cond_m,0) = - dot_product(psi_all_out(cond_m,:),TnPsi)
      enddo

      ! increment in n
      TnppPsi = 0d0
      TnpPsi = TnPsi
      do cond_n = 1,Nc-1
      call op_chebyshev(N,NNZ,A,rp,col,&
          TnppPsi,TnpPsi,TnPsi)
      TnppPsi = TnpPsi
      TnpPsi = TnPsi

      ! calculate mu
      do cond_m = 0,Nc-1
        mu2d(cond_m,cond_n) = - dot_product(psi_all_out(cond_m,:),TnPsi)
      enddo
      enddo





      mu2d_tot = mu2d_tot + mu2d
      mu2d2_tot = mu2d2_tot + mu2d*mu2d
      enddo ! all KPM random repeat

      mu2d_avg = mu2d_tot/(Rep*N)
      mu2d2_avg = mu2d2_tot/(Rep*N)
      deallocate(mu2d_tot,mu2d2_tot,mu2d)
      deallocate(psi0R,psi0,psi1,psi_tmp)
      deallocate(psi0_out,psi1_out,psi_tmp_out)
      deallocate(psi_in,psi_p_in,psi_pp_in)
      deallocate(psi_out,psi_p_out,psi_pp_out)
      deallocate(psi_all_out)
