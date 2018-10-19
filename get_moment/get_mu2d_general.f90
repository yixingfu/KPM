! Created=Wed 13 Dec 2017 01:15:17 PM STD
! Last Modified=Fri 28 Sep 2018 01:57:05 AM DST
      ! This file computes the Jxy moment
      ! for arbitrary direction
      ! dir_a and dir_b
      !  we look for \sigma_ab
      ! calculating \mu_{m,n}
      idum=-(my_id+1)*64

        ! DEBUG
!        open(37,file="ALL.dat",&
!                status="replace",&
!                form="unformatted",access="stream")
!        write(37) N,Nc
!
        ! DEBUG
        
        allocate(JaA(JNNZ),JbA(JNNZ))
        allocate(Jarp(N+1),Jbrp(N+1))
        allocate(Jacol(JNNZ),Jbcol(JNNZ))
      if (Dir_a .eq.DIR_X) then
          JaA = JxA
          Jarp = Jxrp
          Jacol = Jxcol
      else if (Dir_a .eq. DIR_Y) then
          JaA = JyA
          Jarp = Jyrp
          Jacol = Jycol
      endif

      if (Dir_b .eq. DIR_X) then
          JbA = JxA
          Jbrp = Jxrp
          Jbcol = Jxcol
      else if (Dir_b .eq. DIR_Y) then
          JbA = JyA
          Jbrp = Jyrp
          Jbcol = Jycol
      endif


!       do i=1,N
!          do j=rp(i),rp(i+1)-1
!               write(*,"(2I5,2F11.8)") i,col(j),real(A(j)),imag(A(j))
!          enddo
!       enddo
!       write(*,*) "END of H"
!
!       do i=1,N
!          do j=Jarp(i),Jarp(i+1)-1
!          write(*,"(2I5,2F11.8)") i,Jacol(j),real(JaA(j)),imag(JaA(j))
!          enddo
!       enddo
!       write(*,*) "END of Ja"
!       do i=1,N
!          do j=Jbrp(i),Jbrp(i+1)-1
!          write(*,"(2I5,2F11.8)") i,Jbcol(j),real(JbA(j)),imag(JbA(j))
!          enddo
!       enddo
!       write(*,*) "END of Jb"
!
     

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

      mu2d_tot = 0d0
      mu2d2_tot = 0d0


      ! Repetition for KPM
      do i=1,Rep
      mu2d = 0d0

      ! set \psi(0)
!      call ResetRandSeed()
!      call random_number(psi0R)
!      psi0 = zexp(dcmplx(0d0,2.0d0*pi*psi0R))
      do ipsi=1,N
      psirand = ran2(idum)*pi*2d0
      psi0(ipsi) = dcos(psirand)+III*dsin(psirand)
      enddo


!        write(37) real(psi0), imag(psi0)
      TnPsi = psi0

      ! Ja\psi(0)
!      call op_commutator(N,&
!          XYNNZ,RaA,XYrp,XYcol,&
!          NNZ,  A,  rp,  col,&
!          psi0,JaPsi)
      call CSRmultVc16(N,JNNZ,JaA,Jarp,Jacol,psi0,JaPsi)
      TmJaPsi = JaPsi ! m=0
      
      !  Jb Tm(H) Ja\psi(0), m=0
!      call op_commutator(N,&
!          XYNNZ,RbA,XYrp,XYcol,&
!          NNZ,  A,  rp,  col,&
!          TmJaPsi,JbTmJaPsi)
      call CSRmultVc16(N,JNNZ,JbA,Jbrp,Jbcol,TmJaPsi,JbTmJaPsi)
       
      ! \mu_xx_{0,0} = \psi(0)_out \cdot \psi(0)_in
      psi_all_out(0,:) = JbTmJaPsi
      mu2d(0,0) = - dot_product(JbTmJaPsi,TnPsi)

      ! increment in m. For initial, T_{-1} = xT_0
      TmpJaPsi = TmJaPsi
      call CSRmultVc16(N,NNZ,A,rp,col,TmpJaPsi,TmppJaPsi)
      ! the rest then recursive definition.
      do cond_m = 1,Nc-1
      call op_chebyshev(N,NNZ,A,rp,col,&
          TmppJaPsi,TmpJaPsi,TmJaPsi)
      TmppJaPsi = TmpJaPsi
      TmpJaPsi = TmJaPsi

      ! close and save
!      call op_commutator(N,&
!          XYNNZ,RbA,XYrp,XYcol,&
!          NNZ,  A,  rp,  col,&
!          TmJaPsi,JbTmJaPsi)
      call CSRmultVc16(N,JNNZ,JbA,Jbrp,Jbcol,TmJaPsi,JbTmJaPsi)
      psi_all_out(cond_m,:) = JbTmJaPsi
      enddo

      do cond_m = 0,Nc-1
        mu2d(cond_m,0) = - dot_product(psi_all_out(cond_m,:),TnPsi)
!        write(*,*)cond_m,0,mu2d(cond_m,0)
      enddo

      ! increment in n. Similar treatment as m
      TnpPsi = TnPsi
      call CSRmultVc16(N,NNZ,A,rp,col,TnpPsi,TnppPsi)
      ! and the rest follows recursion
      do cond_n = 1,Nc-1
      call op_chebyshev(N,NNZ,A,rp,col,&
          TnppPsi,TnpPsi,TnPsi)
      TnppPsi = TnpPsi
      TnpPsi = TnPsi

      ! calculate mu
      do cond_m = 0,Nc-1
        mu2d(cond_m,cond_n) = - dot_product(psi_all_out(cond_m,:),TnPsi)
!        write(*,*)cond_m,cond_n,mu2d(cond_m,cond_n)
      enddo
      enddo



!        write(37)real(mu2d),imag(mu2d)

      mu2d_tot = mu2d_tot + mu2d/N
!        write(*,*)mu2d_tot(0:4,0:4)/real(i)
!        write(*,*)"----"
      mu2d2_tot = mu2d2_tot + mu2d*mu2d
      enddo ! all KPM random repeat

      mu2d_avg = mu2d_tot/(Rep*N)
      mu2d2_avg = mu2d2_tot/(Rep*N)
      deallocate(mu2d_tot,mu2d2_tot,mu2d)
      deallocate(psi0R,psi0)
      deallocate(psi_all_out)
      deallocate(JaPsi,JbTmJaPsi)
      deallocate(TmJaPsi,TmpJaPsi,TmppJaPsi)
      deallocate(TnPsi,TnpPsi,TnppPsi)
        deallocate(JaA,JbA,Jarp,Jbrp,Jacol,Jbcol)

!        close(37)
