! Last Modified=Wed 30 May 2018 11:08:59 PM DST
! returns 0 if need to be cut
      integer function open_bc(i,i_,L,open_bc_logical)
          integer::i,i_,L,open_bc_logical
          open_bc=1-open_bc_logical*abs(i-i_)/(L-1)
      end function open_bc

      integer function cut_bc_below(i,j,L)
              !!!! Set left bottom cornor to 0
          integer::i,j,L,temp
          temp=((i+j) .le. (L/2))
          cut_bc_below = 1-temp
      end function cut_bc_below


      integer function cut_bc_above(i,j,L)
              !!!! Set right top cornor to 0
          integer::i,j,L,temp
          temp=((i+j) .ge. (3*L/2))
          cut_bc_above=1-temp
      end function cut_bc_above
