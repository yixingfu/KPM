! Last Modified=Thu 25 Oct 2018 01:10:35 AM EDT
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

      integer function set_shape(HC_r,HC_r_,HC_r_ref)
          real*8,dimension(2)::HC_r,HC_r_,HC_r_ref
          integer::temp,temp2
          real*8,parameter::threshold=6400
          real*8::radius
          ! circle with radius 80

          temp = 1
          radius = dot_product(HC_r-HC_r_ref,HC_r-HC_r_ref)
          temp2 = (abs(radius).lt.threshold)
          temp = temp * temp2
          radius = dot_product(HC_r_-HC_r_ref,HC_r_-HC_r_ref)
          temp2 = (abs(radius).lt.threshold)
          temp = temp * temp2
        write(*,*)HC_r,temp

          set_shape = temp
      end function set_shape
