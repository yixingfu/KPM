! Last Modified=Wed 30 May 2018 11:08:59 PM DST
! returns 0 if need to be cut
      integer function open_bc(i,i_,L,open_bc_logical)
          integer::i,i_,L,open_bc_logical
          open_bc=1-open_bc_logical*abs(i-i_)/(L-1)
      end function open_bc
