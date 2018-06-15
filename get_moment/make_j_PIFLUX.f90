! Created=Wed 13 Dec 2017 11:51:36 PM STD
! Last Modified=Thu 14 Jun 2018 01:03:36 PM DST
      ! This file makes J operator for Piflux, 2D

          ! not differentiating QP or not
          if (my_id.eq.0) then
              write(*,*) "2D,J,piflux, testing"
          endif
          ! reset
          Jrp = 0
          Jcol= 0
          col_ind = 1
          rp_ind = 1
          ! txf*(-i),txb*(i), same as the ready use 3D case
!          Jtxf = -tx ! forward is i_= i+1
!          Jtxb = -tx_ !!!DOING it directly.

          do j=1,L!y
          do i=1,L!x
          Jrp(rp_ind) = 2*rp_ind-1
          rp_ind = rp_ind+1

          ! x_forward
          i_=modulo(i,L)+1
          ind_r = xy2i(i_,j,L)
          Jcol(col_ind) = ind_r
          JA(col_ind) = -tx
          col_ind = col_ind+1

          ! x backward
          i_=modulo(i-2,L)+1
          ind_r = xy2i(i_,j,L)
          Jcol(col_ind) = ind_r
          JA(col_ind) = -tx_
          col_ind = col_ind+1

          End do
          End do
          Jrp(rp_ind) = JNNZ+1

