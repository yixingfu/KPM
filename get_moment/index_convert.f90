! created=Tue 12 Dec 2017 03:55:22 PM STD
! Last Modified=Sat 02 Jun 2018 03:05:21 PM DST
      ! changing indices between x,y,z,s to single index
      ! s=0 for spin up, s=1 for spin down
      integer*8 function xyzs2i(x,y,z,s,L)
          integer,intent(in)::x,y,z
          integer,intent(in)::L,s
          xyzs2i = 2_8*((z-1)*L**2_8+(y-1)*L+(x-1))+1+s
          return
      end function xyzs2i
      ! changing indices between x,y,s to single index
      ! s=0 for spin up, s=1 for spin down
      integer*8 function xys2i(x,y,s,L)
          ! For graphene
          ! s = 0: A site (bottom right of hex)
          ! s = 1: B site (bottom left of hex)
          integer,intent(in)::x,y
          integer,intent(in)::L,s
          xys2i = 2_8*((y-1)*L+(x-1))+1+s
          return
      end function xys2i
      integer*8 function xy2i(x,y,L)
          integer,intent(in)::x,y,L
          xy2i = 1_8*((y-1)*L+(x-1))+1
          return
      end function xy2i
      integer*8 function i_up_2d(i,L)
          integer,intent(in)::i,L
          i_up_2d=modulo(i+L*L-L-1,L*L)+1_8
      end function i_up_2d
      integer*8 function i_down_2d(i,L)
          integer,intent(in)::i,L
          i_up_2d=modulo(i+L-1,L*L)+1_8
      end function i_down_2d
      integer*8 function i_left_2d(i,L)
          integer,intent(in)::i,L
          i_up_2d=i-1_8+(((i-1)/L)-((i-2)/L))*L
      end function i_left_2d
      integer*8 function i_right_2d(i,L)
          integer,intent(in)::i,L
          i_up_2d=i+1_8-(((i)/L)-((i-1)/L))*L
      end function i_right_2d
