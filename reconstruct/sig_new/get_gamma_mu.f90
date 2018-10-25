! Cretaed=Sun 23 Sep 2018 10:04:27 PM DST
! Last Modified=Sun 23 Sep 2018 11:10:50 PM DST
      ! This file computes \sum \Gamma_{m,n} \mu_{m,n}
      Gamma_mu = 0d0
      do eps_ind=1,Noutput
      eps = eps_grid(eps_ind)
      Gamma_mu_tmp = 0d0
        write(*,*)"----",eps
      do m=0,Nc-1
      do n=0,m
!        write(*,*) m,n,gamma_mn(m,n,eps),mu_tilde(m,n)
      Gamma_mu_tmp = Gamma_mu_tmp & 
          + real(gamma_mn(m,n,eps) * mu_tilde(m,n))
      enddo
      enddo
      Gamma_mu(eps_ind) = Gamma_mu_tmp
      enddo
