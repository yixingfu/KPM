! Created=Wed 13 Dec 2017 11:51:36 PM STD
! Last Modified=Mon 08 Oct 2018 02:59:25 PM DST
      ! This file makes J operator

      if (D.eq.2) then
          write(*,*) "Testing Conductivity Feature for 2D"
          if (RandType.eq.RandHOP) then
              write(*,*) "Hopping conductivity not implemented"
          else
              if (BHZ) then
                  include "make_j_BHZ.f90"
              else if (PIFLUX) then
                  include "make_j_PIFLUX.f90"
              else if (GRAPHENE) then
                  include "make_j_GRAPHENE.f90"
              else if (IQHE_SQ) then
                  include "make_j_IQHE_SQ.f90"
              else if (IQHE_SOC) then
                  include "make_j_IQHE_SOC.f90"
              else
                  include "make_j_SOC.f90"
              endif
          endif
      else if (D.eq.3) then
          ! not differentiating QP or not
!          write(*,*) "3D,J"
          Jxrp = 0
          Jxcol= 0
          col_ind = 1
          rp_ind = 1
          ! txf*(-i),txb*(i), that's it.
          Jtxf = -txf*III
          Jtxb = txb*III

          do k=1,L!z
          do j=1,L!y
          do i=1,L!x
          do s=0,1
          s_ = 1-s
          Jxrp(rp_ind) = 2*rp_ind-1
          rp_ind = rp_ind+1

          ! x_forward
          ind_r = xyzs2i(modulo(i-2,L)+1,j,k,s_,L)
          Jxcol(col_ind) = ind_r
          JxA(col_ind) = Jtxf(s,s_)
          col_ind = col_ind+1
!        write(*,*) rp_ind-1,',',Jxcol(col_ind-1),&
!        ',',real(JxA(col_ind-1)),',',imag(JxA(col_ind-1))

          ! x backward
          ind_r = xyzs2i(modulo(i,L)+1,j,k,s_,L)
          Jxcol(col_ind) = ind_r
          JxA(col_ind) = Jtxb(s,s_)
          col_ind = col_ind+1
!        write(*,*) rp_ind-1,',',Jxcol(col_ind-1),&
!        ',',real(JxA(col_ind-1)),',',imag(JxA(col_ind-1))

          End do
          End do
          End do
          End do
          Jxrp(rp_ind) = JNNZ+1
      endif

