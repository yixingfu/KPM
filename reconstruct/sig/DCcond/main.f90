

      program main
          integer::m,n,pts,i
          real*8::beta,mu,int_gamma
         write(*,*)"int_gamma_mn(m,n,mu,beta,pts)"
         read(*,*)m,n,mu,beta,pts
        call cpu_time(t1)
          do i=1,10000
         int_gamma = int_gamma_mn(m,n,mu,beta,pts)
         enddo
         call cpu_time(t2)
         write(*,*)(t2-t1)/real(10000)
        
        contains
            include "GammaMN.f90"
            include "Chebyshev.f90"
            include "StatMech.f90"
      end program main
