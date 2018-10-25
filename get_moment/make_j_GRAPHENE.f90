! Created=Wed 13 Dec 2017 11:51:36 PM STD
! Last Modified=Fri 28 Sep 2018 01:52:13 AM DST
      ! This file makes J operator for graphene, 2D

          ! not differentiating QP or not
          if (my_id.eq.0) then
              write(*,*) "2D,J,honeycomb. Not ready - - "
          endif
          ! reset
          Jxrp = 0
          Jxcol= 0
          col_ind = 1
          rp_ind = 1
          ! txf*(-i),txb*(i), same as the ready use 3D case
        Jtexp_theta = III*texp_theta
        Jtexp_theta(0,:) = Jtexp_theta(0,:)&
                *((NNx(0,:)*real_ix+&
                NNy(0,:)*real_jx+real_ABx)*HC_Jx &
                +(NNx(0,:)*real_iy+&
                NNy(0,:)*real_jy+real_ABy)*HC_Jy)
        Jtexp_theta(1,:) = conjg(Jtexp_theta(0,:))
        write(*,*)'--',Jtexp_theta(:,1)
        write(*,*)'--',Jtexp_theta(:,2)
        write(*,*)'--',Jtexp_theta(:,3)


          do j=1,L!y
          do i=1,L!x
          do s=0,1
          s_ = 1-s
          Jxrp(rp_ind) = 3*rp_ind-1
          rp_ind = rp_ind+1

          do NNi=1,3
      i_ = safe_mod(i+NNx(s,NNi),L)
      j_ = safe_mod(j+NNy(s,NNi),L)

          ind_r = xys2i(i_,j_,s_,L)
          Jxcol(col_ind) = ind_r
          JxA(col_ind) = Jtexp_theta(s,NNi)
          col_ind = col_ind+1

          enddo
          End do
          End do
          End do
          Jxrp(rp_ind) = JNNZ+1

