! Created=Wed 13 Dec 2017 11:51:36 PM STD
! Last Modified=Thu 15 Nov 2018 01:00:52 AM STD
      ! This file makes J operator for SOC, 2D

          ! not differentiating QP or not
          if (my_id.eq.0) then
              write(*,*) "2D,J,seem to work"
          endif
          ! reset
          Jxrp = 0
          Jxcol= 0
          col_ind = 1
          rp_ind = 1
          ! txf*(-i),txb*(i), same as the ready use 3D case
          Jtxf = -txf*III
          Jtxb = txb*III

          do j=1,L!y
          do i=1,L!x
          do s=0,1
          s_ = 1-s
          Jxrp(rp_ind) = 2*rp_ind-1
          rp_ind = rp_ind+1

          ! x_forward
          i_=modulo(i,L)+1
          ind_r = xys2i(i_,j,s_,L)
          Jxcol(col_ind) = ind_r
          JxA(col_ind) = Jtxf(s,s_)
          col_ind = col_ind+1

          ! x backward
          i_=modulo(i-2,L)+1
          ind_r = xys2i(i_,j,s_,L)
          Jxcol(col_ind) = ind_r
          JxA(col_ind) = Jtxb(s,s_)
          col_ind = col_ind+1

          End do
          End do
          End do
          Jxrp(rp_ind) = JNNZ+1


          ! reset
          Jyrp = 0
          Jycol= 0
          col_ind = 1
          rp_ind = 1
          ! txf*(-i),txb*(i), same as the ready use 3D case
          Jtyf = -tyf*III
          Jtyb = tyb*III

          do j=1,L!y
          do i=1,L!x
          do s=0,1
          s_ = 1-s
          Jyrp(rp_ind) = 2*rp_ind-1
          rp_ind = rp_ind+1

          ! y_forward
          j_=modulo(j,L)+1
          ind_r = xys2i(i,j_,s_,L)
          Jycol(col_ind) = ind_r
          JyA(col_ind) = Jtyf(s,s_)
          col_ind = col_ind+1

          ! y backward
          j_=modulo(j-2,L)+1
          ind_r = xys2i(i,j_,s_,L)
          Jycol(col_ind) = ind_r
          JyA(col_ind) = Jtyb(s,s_)
          col_ind = col_ind+1

          End do
          End do
          End do
          Jyrp(rp_ind) = JNNZ+1

