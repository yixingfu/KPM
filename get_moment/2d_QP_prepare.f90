! prepare for QP
! This generates correlated potential term
! eps has a term on each site
! Fully correlated version

      include "QP_prepare_public.f90"

      allocate(eps(L*L)) 
      eps_ind = 1
      if (RandType.eq.RandPOT) then
          do j=1,L
          do i=1,L

          eps(eps_ind) = &
              WQP*quasiperiodic(i+0d0,j+0d0,P,Q,phase)&
              + Wrnd*random2D(i+0d0,j+0d0,P,Q)
          eps_ind=eps_ind+1
          enddo
          enddo
          eps = eps - sum(eps)/real(L*L)

        if (Scramble) then
                write(*,*) "scrambled"
                call ScrambleArray(eps)
        endif
      else if (RandType.eq.RandHOP) then
          eps = 0d0
      endif


