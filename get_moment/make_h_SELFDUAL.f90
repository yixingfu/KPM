! Created=Tue 12 Dec 2017 03:28:22 PM STD
! Last Modified=Tue 23 Oct 2018 12:49:11 PM DST
      !This file creates H for self dual model
      !The matrix is stored as CSR(A,col,rp)
      !
      allocate(A(NNZ),col(NNZ),rp(N+1))
      if (fixedTwist) then
          Twist = OrigTwist*pi/real(L)
      else

          allocate(TwistAll(num_procs*3*seq_rep))
          call random_number(TwistAll)
          Twist=TwistAll(my_id*3+1:my_id*3+3)
          twist=twist*2d0-1d0
          write(*,*) 'twist:',Twist
          deallocate(TwistAll)

          Twist = (Twist)*pi/L! 0 to pi??
      endif



      ! Add Twist
      txf = xf*cdexp(dcmplx(0d0,Twist(1)))
      txb = xb*cdexp(dcmplx(0d0,-Twist(1)))
      tyf = yf*cdexp(dcmplx(0d0,Twist(2)))
      tyb = yb*cdexp(dcmplx(0d0,-Twist(2)))
      tzf = zf*cdexp(dcmplx(0d0,Twist(3)))
      tzb = zb*cdexp(dcmplx(0d0,-Twist(3)))
      if (my_id .eq. 0) then
          write(*,*) '-',txf
          write(*,*) '-',txb
          write(*,*) '-',tyf
          write(*,*) '-',tyb
          write(*,*) '-',tzf
          write(*,*) '-',tzb
      endif






      if (D.eq.3) then

          include "QP_prepare.f90"
          ! initiate indices
          eps_ind = 1
          rp = 0
          col = 0
          col_ind = 1
          rp_ind = 1
          
          do k=1,L! z
          do j=1,L! y
          do i=1,L! x
          do s = 0,1
          s_ = 1-s

          ! set row pointer
          rp(rp_ind) = PerSite*rp_ind-(PerSite-1)
          rp_ind = rp_ind+1

          !!! Disorder term
          ! same spin
          col(col_ind) = xyzs2i(i,j,k,s,L) ! = i,j,k,s,L
          A(col_ind) = WQP * ( &
              (dcos(Q*i+phase(1))*pauli_z(s,s)) + &
              (dcos(P*j+phase(2))*pauli_x(s,s)) + &
              (dcos(R*k+phase(3))*pauli_y(s,s)))
          col_ind = col_ind+1
          
          ! opposite spin
          col(col_ind) = xyzs2i(i,j,k,s_,L)
          A(col_ind) = WQP * ( &
              (dcos(Q*i+phase(1))*pauli_z(s,s_)) + &
              (dcos(P*j+phase(2))*pauli_x(s,s_)) + &
              (dcos(R*k+phase(3))*pauli_y(s,s_)))
          col_ind = col_ind+1

          !!! Hopping term
          ! (x,x+1) (xf)
          i_=modulo(i,L)+1
          ind_r = xyzs2i(i_,j,k,s_,L)
          t_tmp = t0 
          col(col_ind) = ind_r
          A(col_ind) = txf(s,s_)*t_tmp*open_bc(i,i_,L,OPEN_BC_x)
          col_ind = col_ind+1

          ! (x,x-1) (xb)
          i_=modulo(i-2,L)+1
          ind_r = xyzs2i(i_,j,k,s_,L)
          t_tmp = t0 
          col(col_ind) = ind_r
          A(col_ind) = txb(s,s_)*t_tmp*open_bc(i,i_,L,OPEN_BC_x)
          col_ind = col_ind+1

          ! (y,y+1) (yf)
          j_=modulo(j,L)+1
          ind_r = xyzs2i(i,j_,k,s_,L)
          t_tmp = t0 
          col(col_ind) = ind_r
          A(col_ind) = tyf(s,s_)*t_tmp*open_bc(j,j_,L,OPEN_BC_y)
          col_ind = col_ind+1

          ! (y,y-1) (yb)
          j_=modulo(j-2,L)+1
          ind_r = xyzs2i(i,j_,k,s_,L)
          t_tmp = t0 
          col(col_ind) = ind_r
          A(col_ind) = tyb(s,s_)*t_tmp*open_bc(j,j_,L,OPEN_BC_y)
          col_ind = col_ind+1

          ! (z,z+1) (zf)
          k_=modulo(k,L)+1
          ind_r = xyzs2i(i,j,k_,s,L)
          t_tmp = t0 
          col(col_ind) = ind_r
          A(col_ind) = tzf(s,s)*t_tmp*open_bc(k,k_,L,OPEN_BC_z)
          col_ind = col_ind+1

          ! (z,z-1) (zb)
          k_=modulo(k-2,L)+1
          ind_r = xyzs2i(i,j,k_,s,L)
          t_tmp = t0 
          col(col_ind) = ind_r
          A(col_ind) = tzf(s,s)*t_tmp*open_bc(k,k_,L,OPEN_BC_z)
          col_ind = col_ind+1


          End do
          End do
          End do
          rp(rp_ind) = NNZ+1

          deallocate(eps,phase_all)
      else if (D.eq.2) then
          write(*,*) "not implemented - self dual 2D"
          

      endif
      if (my_id.eq.0) then
          open(111,file="MATRIX.txt",status="replace",&
              form="unformatted",access="stream")
          write(111) N
          write(111) NNZ
          write(*,*)"NNZ",NNZ
          write(111) A, col, rp
          close(111)
      endif



! -------------------------EIGENVALUE
      include "get_eigenvalue.f90"
! --------------------EIGEnVALUE end

      ! normalization same for both

      call LanczosBound(N,NNZ,A,rp,col,1000,Emax,Emin)
      norm_a = (Emax-Emin)/(2d0-0.2d0)
      norm_b = (Emax+Emin)/2
        if (Set_norm_a .ne. 0) then
                norm_a = Set_norm_a
                norm_b = 0d0
        endif
      call rescale_to_1(N,NNZ,A,rp,col,norm_a,norm_b)


