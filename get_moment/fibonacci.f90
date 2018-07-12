      integer function fibonacci(N,FiboBasis_in)
          integer::N,i
          integer::fi,fi_1,fi_2
          integer,dimension(2),optional::FiboBasis_in
          integer,dimension(2)::FiboBasis
          if (present(FiboBasis_in)) then
                FiboBasis = FiboBasis_in
          else
                FiboBasis(1) = 0
                FiboBasis(2) = 1
          endif
          if (N .eq. 0) then
              fibonacci = FiboBasis(1)
          else if (N .eq. 1) then
              fibonacci = FiboBasis(2)
          else
              fi_2 = 0
              fi_1 = 1
              do i=2,N
              fi=fi_1+fi_2
              fi_2=fi_1
              fi_1=fi
              enddo
              fibonacci = fi
          endif
          return

      end function fibonacci

      integer function inv_fibo(L,FiboBasis_in)
                integer::N,L
                integer,dimension(2),optional::FiboBasis_in
                integer,dimension(2)::FiboBasis
          if (present(FiboBasis_in)) then
                FiboBasis = FiboBasis_in
          else
                FiboBasis(1) = 0
                FiboBasis(2) = 1
          endif
                do N=1,L
                    if (fibonacci(N,FiboBasis).eq.L) then
                        inv_fibo=N
                        return
                    endif
                enddo
      end function inv_fibo
