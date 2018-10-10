! Last Modified=Sun 23 Sep 2018 11:24:23 PM DST
      ! get \tilde{\mu} by applying kernel
        open(52,file=trim(inputfile_k)//"mu.txt",status="replace")
      do m=0,Nc-1
      do n=0,Nc-1
      mu_tilde(m,n) = mu_avg(m,n)*hm(m)*hm(n)*gJ(m)*gJ(n)
        write(52,*)m,n,mu_avg(m,n)
      enddo
      enddo
        close(52)

