! prepare for QP  -- 3D version
! This generates correlated potential term
! eps has a term on each site
! Fully correlated version

      include "QP_prepare_public.f90"
      


      allocate(eps(L*L*L)) 
      eps_ind = 1

        ! if fully random, do it now
        call random3D(eps,Wrnd,my_id)
        ! no matter this or that, then go through
        ! QP.  saving trouble at the price of 
        ! a bit of CPU time
        ! this is already normalized

      if (RandType.eq.RandPOT) then
          do k=1,L
          do j=1,L
          do i=1,L
          eps(eps_ind) = eps(eps_ind)+&
              WQP*quasiperiodic3D(i+0d0,j+0d0,k+0d0,P,Q,R,phase)
          eps_ind=eps_ind+1
          enddo
          enddo
          enddo
          eps = eps - sum(eps)/real(L*L*L)
      else if (RandType.eq.RandHOP) then
          eps = 0d0
      endif


