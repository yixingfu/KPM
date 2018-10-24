! Last Created=Tue 23 Oct 2018 06:32:55 PM DST
! Last Modified=Wed 24 Oct 2018 12:48:46 AM DST
      ! This file computes IPR based on the ED states saved in H_dense
      ! as the first index (column) vectors. IPRx is direct, while IPRk
      ! is done after a fourier transformation
      if (seq_i .eq. 0) then
        write(*,*) "IPR CALCULATION"
          allocate(IPRx_all(N),IPRk_all(N),IPR_E(N))
          IPRx_all = 0d0
          IPRk_all = 0d0
          IPRcount = 0
      endif
        allocate(IPRx(N),IPRk(N),psi_x(N),psi_k(N))
      do i=1,N
        IPR_E(i) = EigVal(i)
        psi_x = H_dense(:,i)

        ! real space
        IPRx(i) = sum(zabs(psi_x)**4d0)&
            /(sum(zabs(psi_x)**2d0)**2d0)

        ! momentum space
        if (TWOSPIN3D) then
            allocate(psi_x_up(L*L*L),psi_x_down(L*L*L))
            allocate(psi_k_up(L*L*L),psi_k_down(L*L*L))
            psi_x_up = psi_x(1:N:2)
            psi_x_down = psi_x(2:N:2)
            call FFT_3D(L,L,L,L*L*L,psi_x_up,psi_k_up)
            call FFT_3D(L,L,L,L*L*L,psi_x_down,psi_k_down)
            psi_k(1:N:2) = psi_k_up
            psi_k(2:N:2) = psi_k_down
            deallocate(psi_x_up,psi_x_down,psi_k_up,psi_k_down)
        else 
            write(*,*) "not implemented"
        endif
        IPRk(i) = sum(zabs(psi_k)**4d0)&
            /(sum(zabs(psi_k)**2d0)**2d0)

      end do
      IPRx_all = IPRx_all + IPRx
      IPRk_all = IPRk_all + IPRk
      IPRcount = IPRcount + 1
      deallocate(IPRx,IPRk,psi_x,psi_k)
