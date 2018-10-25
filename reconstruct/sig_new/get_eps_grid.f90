! Created=Sun 23 Sep 2018 10:22:47 PM DST
! Last Modified=Sun 23 Sep 2018 11:02:49 PM DST
      ! This file creates the grid to be used for eps
      allocate(eps_grid(Noutput))
      allocate(gamma_mu(Noutput))
      eps_grid = 0d0
      do i=1,Noutput
        eps_grid(i) = 0.99d0*(i/real(Noutput)*2d0-1d0)
      enddo

