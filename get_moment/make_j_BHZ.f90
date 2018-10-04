! Created=Wed 13 Dec 2017 11:51:36 PM STD
! Last Modified=Fri 28 Sep 2018 01:47:09 AM DST
      ! This file makes Jx and Jy operator for BHZ, 2D

          ! not differentiating QP or not
          if (my_id.eq.0) then
              write(*,*) "2D,J,BHZ, working"
          endif
          ! reset
          Jxrp = 0
          Jxcol= 0
          col_ind = 1
          rp_ind = 1
          ! txf*(-i),txb*(i), same as the ready use 3D case
          Jtxf = -txf*III
          Jtxb = txb*III
        write(*,*)"Jtxf"
        write(*,*)Jtxf
        write(*,*)"Jtxb"
        write(*,*)Jtxb
        write(*,*)"----"

          do j=1,L!y
          do i=1,L!x
          do s=0,1
          s_ = 1-s
          Jxrp(rp_ind) = 2*2*rp_ind-3
          ! each real space point has 4: bwd/fwd,same/diff spin
          rp_ind = rp_ind+1

          ! x_forward
          i_=modulo(i-2,L)+1
          ind_r = xys2i(i_,j,s_,L)
          Jxcol(col_ind) = ind_r
          JxA(col_ind) = Jtxf(s,s_)*open_bc(i,i_,L,OPEN_BC_x)
          col_ind = col_ind+1
          ind_r =  xys2i(i_,j,s,L)
          Jxcol(col_ind) = ind_r
          JxA(col_ind) = Jtxf(s,s)*open_bc(i,i_,L,OPEN_BC_x)
          col_ind = col_ind+1


          ! x backward
          i_=modulo(i,L)+1
          ind_r = xys2i(i_,j,s_,L)
          Jxcol(col_ind) = ind_r
          JxA(col_ind) = Jtxb(s,s_)*open_bc(i,i_,L,OPEN_BC_x)
          col_ind = col_ind+1
          ind_r = xys2i(i_,j,s,L)
          Jxcol(col_ind) = ind_r
          JxA(col_ind) = Jtxb(s,s)*open_bc(i,i_,L,OPEN_BC_x)
          col_ind = col_ind+1

          End do
          End do
          End do
          Jxrp(rp_ind) = JNNZ+1


          ! Now Jy

          ! reset
          Jyrp = 0
          Jycol= 0
          col_ind = 1
          rp_ind = 1
          ! txf*(-i),txb*(i), same as the ready use 3D case
          Jtyf = -tyf*III
          Jtyb = tyb*III
        write(*,*)"Jtyf"
        write(*,*)Jtyf
        write(*,*)"Jtyb"
        write(*,*)Jtyb
        write(*,*)"----"

          do j=1,L!y
          do i=1,L!x
          do s=0,1
          s_ = 1-s
          Jyrp(rp_ind) = 2*2*rp_ind-3
          ! each real space point has 4: bwd/fwd,same/diff spin
          rp_ind = rp_ind+1

          ! y_forward
          j_=modulo(j-2,L)+1
          ind_r = xys2i(i,j_,s_,L)
          Jycol(col_ind) = ind_r
          JyA(col_ind) = Jtyf(s,s_)*open_bc(j,j_,L,OPEN_BC_y)
          col_ind = col_ind+1
          ind_r =  xys2i(i,j_,s,L)
          Jycol(col_ind) = ind_r
          JyA(col_ind) = Jtyf(s,s)*open_bc(j,j_,L,OPEN_BC_y)
          col_ind = col_ind+1


          ! y backward
          j_=modulo(j,L)+1
          ind_r = xys2i(i,j_,s_,L)
          Jycol(col_ind) = ind_r
          JyA(col_ind) = Jtyb(s,s_)*open_bc(j,j_,L,OPEN_BC_y)
          col_ind = col_ind+1
          ind_r = xys2i(i,j_,s,L)
          Jycol(col_ind) = ind_r
          JyA(col_ind) = Jtyb(s,s)*open_bc(j,j_,L,OPEN_BC_y)
          col_ind = col_ind+1

          End do
          End do
          End do
          Jyrp(rp_ind) = JNNZ+1


