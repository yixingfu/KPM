! Branched=Wed 14 Nov 2018 02:58:10 PM STD
! Last Modified=Thu 15 Nov 2018 01:07:58 AM STD
      ! This file computes mu, that comes from (J^2 \rho)
        if (my_id.eq.0) then
            write(*,*) "2D only - mu with J2"
!        write(*,*)"J"
!        call CSR_print(N,JNNZ,JxA,Jxrp,Jxcol)
        write(*,*)"H"
        call CSR_print(N,NNZ,A,rp,col)
        endif
      allocate(mu_tot(0:Nc-1),mu2_tot(0:Nc-1))
      allocate(mu_avg(0:Nc-1),mu2_avg(0:Nc-1))
      allocate(mu_J2_tot(0:Nc-1),mu_J22_tot(0:Nc-1))
      allocate(mu_J2_avg(0:Nc-1),mu_J22_avg(0:Nc-1))
      allocate(mu(0:Nc-1))
      allocate(mu_J2(0:Nc-1))
      allocate(psi0R(N),psi0(N),psi1(N),psi_tmp(N))
      allocate(psi0J(N),psi1J(N))
      allocate(psi_j(N),psi_j_p(N),psi_j_pp(N))
      allocate(psi_jJ(N),psi_j_pJ(N),psi_j_ppJ(N))
      mu_tot = 0
      mu2_tot = 0
!        call ResetRandSeed()
      idum=-(rlz_id+1)*64
!      randomtest1 = ran2(idum)
!      randomtest2 = ran2(idum)
!      write(*,*) randomtest1,randomtest2

      do i=1,Rep
      ! for each repetition
      ! first, create random vector psi0(random phase)
      do ipsi=1,N
      psirand = ran2(idum)*pi*2d0
      psi0(ipsi) = dcos(psirand)+III*dsin(psirand)
      enddo
      ! psi0 gives mu0
      mu = 0
      mu_J2 = 0
      mu(0) = dot_product(psi0,psi0)
        
!      call CSRmultVc16(N,JNNZ,JxA,Jxrp,Jxcol,psi0,psi_tmp)
      call CSRmultVc16(N,JNNZ,JxA,Jxrp,Jxcol,psi0,psi0J)
!      call CSRmultVc16(N,JNNZ,JyA,Jyrp,Jycol,psi0,psi_tmp)
!      call CSRmultVc16(N,JNNZ,JyA,Jyrp,Jycol,psi_tmp,psi0J2y)
      mu_J2(0) = dot_product(psi0,psi0J)

      ! Then calculate psi1
      call CSRmultVc16(N,NNZ,A,rp,col,psi0,psi1)
      call CSRmultVc16(N,NNZ,A,rp,col,psi0J,psi1J)
      ! psi1 gives mu1
      mu(1) = dot_product(psi0,psi1)
      mu_J2(1) = dot_product(psi0,psi1J)

      psi_j_pp = psi0
      psi_j_p  = psi1
      psi_j_ppJ = psi0J
      psi_j_pJ  = psi1J
      j0 = 2

      do j=j0,Nc-1
      ! \psi_j = 2H \psi_{j-1} - \psi_{j-2}
      call CSRmultVc16(N,NNZ,A,rp,col,psi_j_p,psi_tmp)
      psi_j = 2d0*psi_tmp-psi_j_pp
      call CSRmultVc16(N,NNZ,A,rp,col,psi_j_pJ,psi_tmp)
      psi_jJ = 2d0*psi_tmp-psi_j_ppJ

      ! mu(j)
      mu(j) = dot_product(psi0,psi_j)
      mu_J2(j) = dot_product(psi0,psi_jJ)

      ! renew psi_j_pp and psi_j_p
      psi_j_pp = psi_j_p
      psi_j_p = psi_j
      psi_j_ppJ = psi_j_pJ
      psi_j_pJ = psi_jJ
      enddo
      mu = mu/N
      mu_J2 = mu_J2/N
      if (task .eq. RHODER) then
          write(*,*) "THIS IS ABANDONED"
      endif

!        write(*,*)mu(0:5)
      mu_tot = mu_tot+mu
      mu_J2_tot = mu_J2_tot+mu_J2
        write(*,*)mu_tot(0:4)/real(i)
        write(*,*)mu_J2_tot(0:4)/real(i)
        write(*,*)"-----"
      mu2_tot = mu2_tot+mu**2
      mu_J22_tot = mu_J22_tot+mu_J2**2
      enddo
      mu_avg = mu_tot/(Rep)
      mu_J2_avg = mu_J2_tot/(Rep)
      mu2_avg = mu2_tot/(Rep)
      mu_J22_avg = mu_J22_tot/(Rep)
      ! divide by N as required by DoS
      ! This is because the delta function rep of rho has 1/N,
      ! dividing by total # of states.
      deallocate(mu_tot,mu2_tot,mu)
      deallocate(mu_J2_tot,mu_J22_tot,mu_J2)
      deallocate(psi0R,psi0,psi1)
      deallocate(psi_j,psi_j_p,psi_j_pp,psi_tmp)
      deallocate(psi_jJ,psi_j_pJ,psi_j_ppJ)
      deallocate(psi0J,psi1J)
