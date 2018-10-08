! Created=Wed 13 Dec 2017 11:51:36 PM STD
! Last Modified=Thu 04 Oct 2018 05:45:15 PM DST
      ! This file makes Jx and Jy operator for BHZ, 2D

          ! not differentiating QP or not
          if (my_id.eq.0) then
              write(*,*) "2D,J,IQHE_SQ, working"
          endif
          ! reset
          Jxrp = 0
          Jxcol= 0
          col_ind = 1
          rp_ind = 1


          do j=1,L!y
          do i=1,L!x
          Jxrp(rp_ind) = 2*rp_ind-1
          ! each real space point has 2: bwd/fwd
          rp_ind = rp_ind+1

          ! x_forward
          i_=modulo(i-2,L)+1
          ind_r = xy2i(i_,j,L)
          Jxcol(col_ind) = ind_r
          t_tmp = t0 + Trnd*random2D(i+0.5d0,j+0d0,P,Q)&
              + TQP*quasiperiodic(i+0.5d0,j+0d0,P,Q,phase)
          JxA(col_ind) = III * cdexp(dcmplx(0d0,Twist(1)))&
              *t_tmp*open_bc(i,i_,L,OPEN_BC_x)
          col_ind = col_ind+1


          ! x backward
          i_=modulo(i,L)+1
          ind_r = xy2i(i_,j,L)
          Jxcol(col_ind) = ind_r
          t_tmp = t0 + Trnd*random2D(i+0.5d0,j+0d0,P,Q)&
              + TQP*quasiperiodic(i+0.5d0,j+0d0,P,Q,phase)
          JxA(col_ind) = - III * cdexp(dcmplx(0d0,-Twist(1)))&
              *t_tmp*open_bc(i,i_,L,OPEN_BC_x)

          col_ind = col_ind+1

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

          do j=1,L!y
          do i=1,L!x
          Jyrp(rp_ind) = 2*rp_ind-1
          ! each real space point has 4: bwd/fwd
          rp_ind = rp_ind+1

          ! y_forward
          j_=modulo(j-2,L)+1
          ind_r = xy2i(i,j_,L)
          Jycol(col_ind) = ind_r
          t_tmp = t0 + Trnd*random2D(i+0.0d0,j-0.5d0,P,Q)&
              + TQP*quasiperiodic(i+0.0d0,j-0.5d0,P,Q,phase)
          JyA(col_ind) = III &
              *cdexp(dcmplx(0d0, Twist(2)+2d0*pi*i*IQHE_B)) &
              *t_tmp*open_bc(j,j_,L,OPEN_BC_y)
          col_ind = col_ind+1


          ! y backward
          j_=modulo(j,L)+1
          ind_r = xy2i(i,j_,L)
          Jycol(col_ind) = ind_r
          t_tmp = t0 + Trnd*random2D(i+0.0d0,j+0.5d0,P,Q)&
              + TQP*quasiperiodic(i+0.0d0,j+0.5d0,P,Q,phase)
          JyA(col_ind) = - III &
              *cdexp(dcmplx(0d0,-Twist(2)-2d0*pi*i*IQHE_B)) &
              *t_tmp*open_bc(j,j_,L,OPEN_BC_y)
          col_ind = col_ind+1

          End do
          End do
          Jyrp(rp_ind) = JNNZ+1


