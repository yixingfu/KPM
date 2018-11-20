! Created=Wed 13 Dec 2017 03:05:25 PM STD
! Last Modified=Fri 07 Sep 2018 11:02:22 AM DST
      ! 
      allocate(rho(1:Ntilde))
      allocate(rho_J2(1:Ntilde))
      ! it need to be scaled here, and scaled back later
      allocate(Egrid_t(1:Ntilde))
      Egrid_t = (Egrid-norm_b)/norm_a
        rho = 0
        rho_J2 = 0


      if (useFFT) then
          write(*,*) "Not implemented yet"
      else 
          if (Der .eq. 0) then
              do i=1,Ntilde
              !1/(pi sqrt(1-x^2)) (g0 mu0+2 sum (g_n mu_n Tn(x)))
              Ei = Egrid_t(i)
              rhoi = gJ(0)*mu_avg(0)
              rho_J2i = gJ(0)*mu_J2_avg(0)
              do j=1,ForceNc-1
              rhoi = rhoi+&
                  2d0*gJ(j)*mu_avg(j)*ChebyT(j,Ei)
              rho_J2i = rho_J2i+&
                  2d0*gJ(j)*mu_J2_avg(j)*ChebyT(j,Ei)
              End do
                rhoi = rhoi/(norm_a*pi*dsqrt(1d0-Ei*Ei))
              rho(i) = rhoi
!              rho_J2(i) = rho_J2i/(norm_a*pi*dsqrt(1d0-Ei*Ei))
                rho_J2i=rho_J2i/(norm_a*pi*dsqrt(1d0-Ei*Ei))
              rho_J2(i) = rho_J2i/(rhoi+0.000000001)
        write(*,*)Ei,rhoi,rho_J2i,rho_J2i/rhoi
              End do
          else if (Der .eq. 1) then
              do i=1,Ntilde
              Ei = Egrid_t(i)
              rhoi = gJ(0)*mu_avg(0)*(1d0/(1d0-Ei*Ei))
              ! T_0(x)=1,U_-1(x)=0
              do j=1,ForceNc-1
              rhoi = rhoi+&
                  2d0*gJ(j)*mu_avg(j)*(ChebyT(j,Ei)/(1d0-Ei*Ei) &
                  + ChebyU(j-1,Ei)*j)
              enddo
              rho(i) = rhoi/(norm_a*pi*dsqrt(1d0-Ei*Ei))
              enddo
          else if (Der .eq. 2) then
              write(*,*) "The Entire Array is E=0!"
              rhoi = D2rhomu0(norm_a,norm_b,0d0)&
                  * gJ(0)*mu_avg(0)
              do j=1,ForceNc-1
              rhoi = rhoi+D2rho(j,norm_a,norm_b,0d0) &
                  *2d0* gJ(j)*mu_avg(j)
              enddo

!              rhoi = gJ(0)*mu_avg(0)
!              do j=1,ForceNc-1
!                rhoi = rhoi + &
!                    2d0*dcos(j*pi/2d0)*gJ(j)*mu_avg(j) + &
!                    2d0*j*j*dcos(pi*(j-2d0)/2d0)*gJ(j)*mu_avg(j)
!              enddo
!              rhoi = rhoi/(norm_a**3*pi)
              do i=1,Ntilde
              rho(i) = rhoi
              enddo
          else if (Der .eq. 4) then
              write(*,*) "The Entire Array is E=0!"
              rhoi = D4rhomu0(norm_a,norm_b,0d0)&
                  * gJ(0)*mu_avg(0)
              do j=1,ForceNc-1
              rhoi = rhoi+D4rho(j,norm_a,norm_b,0d0) &
                  * 2d0 * gJ(j) * mu_avg(j)
              enddo
              do i=1,Ntilde
              rho(i) = rhoi
              enddo

          endif

      endif
      rho_tot = rho_tot+(rho)
      rho_J2_tot = rho_J2_tot+(rho_J2)
      rho2_tot = rho2_tot+rho*rho
      rho_J22_tot = rho_J22_tot+rho_J2*rho_J2
      deallocate(rho,rho_J2,Egrid_t)
      deallocate(mu_avg,mu2_avg)! no longer used
      deallocate(mu_J2_avg,mu_J22_avg)! no longer used

