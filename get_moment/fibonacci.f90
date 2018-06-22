      integer function fibonacci(N)
          integer::N,i
          integer::fi,fi_1,fi_2
          if (N .eq. 0) then
              fibonacci = 0
          else if (N .eq. 1) then
              fibonacci = 1
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

      integer function inv_fibo(L)
                integer::N,L
                do N=1,L
                    if (fibonacci(N).eq.L) then
                        inv_fibo=N
                        return
                    endif
                enddo
      end function inv_fibo
