! Last Modified=Sun 23 Sep 2018 11:24:23 PM DST
      ! get \tilde{\mu} by applying kernel
      do m=0,Nc-1
      do n=0,Nc-1
      mu_tilde(m,n) = mu_avg(m,n)*hm(m)*hm(n)*gJ(m)*gJ(n)
        write(*,*)m,n,mu_tilde(m,n)
      enddo
      enddo

