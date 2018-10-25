! Created=Wed 13 Dec 2017 11:51:36 PM STD
! Last Modified=Tue 23 Oct 2018 12:08:47 PM DST
      ! This file makes J operator for SOC, 2D

          ! not differentiating QP or not
          if (my_id.eq.0) then
              write(*,*) "SELF DUAL, conductivity not implemented"
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

