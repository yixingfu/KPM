! Created=Wed 13 Dec 2017 11:51:36 PM STD
! Last Modified=Wed 19 Sep 2018 11:26:07 AM DST
      ! This file makes J operator for BHZ, 2D

          ! not differentiating QP or not
          if (my_id.eq.0) then
              write(*,*) "2D,J,BHZ,not ready"
          endif
          ! reset
          Jrp = 0
          Jcol= 0
          col_ind = 1
          rp_ind = 1
          ! txf*(-i),txb*(i), same as the ready use 3D case
          Jtxf = -txf*III
          Jtxb = txb*III

          do j=1,L!y
          do i=1,L!x
          do s=0,1
          s_ = 1-s
          Jrp(rp_ind) = 2*2*rp_ind-1
          ! each real space point has 4
          rp_ind = rp_ind+1

          ! x_forward
          i_=modulo(i,L)+1
          ind_r = xys2i(i_,j,s_,L)
          Jcol(col_ind) = ind_r
          JA(col_ind) = Jtxf(s,s_)*open_bc(i,i_,L,OPEN_BC_x)
          col_ind = col_ind+1
          ind_r =  xys2i(i_,j,s,L)
          Jcol(col_ind) = ind_r
          JA(col_ind) = Jtxf(s,s)*open_bc(i,i_,L,OPEN_BC_x)
          col_ind = col_ind+1


          ! x backward
          i_=modulo(i-2,L)+1
          ind_r = xys2i(i_,j,s_,L)
          Jcol(col_ind) = ind_r
          JA(col_ind) = Jtxb(s,s_)*open_bc(i,i_,L,OPEN_BC_x)
          col_ind = col_ind+1
          ind_r = xys2i(i_,j,s,L)
          col(col_ind) = ind_r
          JA(col_ind) = Jtxb(s,s)*open_bc(i,i_,L,OPEN_BC_x)
          col_ind = col_ind+1

          End do
          End do
          End do
          Jrp(rp_ind) = JNNZ+1

