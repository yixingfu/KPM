! Created=Tue 12 Dec 2017 06:05:12 PM STD
! Last Modified=Thu 12 Jul 2018 10:41:12 AM DST
      ! some random or quasiperiodic number routine

      real*8 function quasiperiodic(i,j,P,Q,phase)
          real*8::i,j
          real*8::P,Q
          real*8,dimension(2)::phase
          quasiperiodic = dcos(P*i+phase(1))+&
              dcos(Q*j+phase(2))
          return
      end function quasiperiodic

      real*8 function quasiperiodic_HC(i,j,P,Q,phase)
          real*8::i,j
          real*8::P,Q
          real*8,dimension(3)::phase
          quasiperiodic_HC = dcos(P*i+phase(1))+&
              dcos(Q*j+phase(2))+&
              dcos(-dsqrt(3d0)/2d0*Q*(j+i)+phase(3))
          return
      end function quasiperiodic_HC

      real*8 function quasiperiodic3D(i,j,k,P,Q,R,phase)
          real*8::i,j,k
          real*8::P,Q,R
          real*8,dimension(3)::phase
          quasiperiodic3D = dcos(P*i+phase(1))+&
              dcos(Q*j+phase(2))+&
              dcos(R*k+phase(3))
          return
      end function quasiperiodic3D


      real*8 function random2D(i,j,P,Q)
          ! This is the random that fits 2D case best
          real*8::i,j
          real*8::P,Q
          real*8,dimension(2)::phase
          real*8,parameter::pi=3.1415926535897932384626433832795d0
          call random_number(phase)
          phase = phase*2d0*pi
          random2D = dcos(P*i+phase(1))+&
              dcos(Q*j+phase(2))
          return
      end function random2D

      subroutine random3D(eps,W,my_id)
          real*8,dimension(:),intent(inout)::eps
          real*8,intent(in)::W
          real*8,dimension(:),allocatable::U1,U2
          integer::N,i
          integer,intent(in)::my_id
          real*8,parameter::pi=3.1415926535897932384626433832795d0
          N = size(eps)
          allocate(U1(N),U2(N))
          do i=0,2*my_id+3
          call random_number(U1)
          call random_number(U2)
          enddo
          eps = dsqrt(-2d0*dlog(U1))*dcos(2d0*pi*U2)*W
          eps = eps - sum(eps)/N!shift
          !write(*,*)eps
          deallocate(U1,U2)
          return
      end subroutine random3D
      subroutine ResetRandSeed(SeedOffset)
          integer :: N,i
          integer,allocatable,dimension(:)::seed
          integer::un
          integer,optional,intent(in)::SeedOffset
          call random_seed(size=N)
          allocate(seed(N))
          if (present(SeedOffset)) then
              seed = SeedOffset*(/ (i,i=1,n) /)
          else
              open(newunit=un, file="/dev/urandom", access="stream", &
                  form="unformatted", action="read",&
                  status="old")
              read(un) seed
              close(un)
          endif
          call random_seed(put=seed)
      End subroutine ResetRandSeed



