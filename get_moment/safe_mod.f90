! Created=Wed 06 Jun 2018 04:37:07 PM DST
! Last Modified=Wed 06 Jun 2018 04:38:43 PM DST
      ! Gives a modulo that ends in 1 to N
      integer function safe_mod(i,L)
          integer,intent(in)::i,L
          safe_mod=mod(i+L-1,L)+1
      end function safe_mod
