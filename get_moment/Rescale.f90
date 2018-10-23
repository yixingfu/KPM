! Created=Mon 29 Jan 2018 02:45:12 PM EST
! Last Modified=Thu 10 May 2018 03:05:35 PM DST
! This file gives rescale back and forth.


      subroutine rescale_to_1(N,NNZ,A,rp,col,norm_a,norm_b)
          integer*8,intent(in)::N,NNZ
          integer*8,dimension(NNZ),intent(in)::col
          integer*8,dimension(N+1),intent(in)::rp
          complex*16,dimension(NNZ),intent(inout)::A
          real*8,intent(in)::norm_a,norm_b
          integer*8::i,j
        complex*16::Arpi

        do i=1,N
                do j=rp(i),rp(i+1)-1
        write(*,'(i4.4,a,i4.4,a,F11.9,a,F11.9)')i,',',col(j),&
                ',',real(A(j)),',',imag(A(j))
                enddo
        enddo

          do i=1,N
        Arpi = A(rp(i))
          A(rp(i)) = Arpi-norm_b*(abs(Arpi).gt.0.000001d0)
          enddo
          A = A/norm_a
          write(*,*)'normalization factors:',norm_a,norm_b

          return
      end subroutine rescale_to_1

      subroutine rescale_back(N,NNZ,A,rp,col,norm_a,norm_b)
          integer*8,intent(in)::N,NNZ
          integer*8,dimension(NNZ),intent(in)::col
          integer*8,dimension(N+1),intent(in)::rp
          complex*16,dimension(NNZ),intent(inout)::A
          real*8,intent(in)::norm_a,norm_b
          integer*8::i
          A = A*norm_a
          do i=1,N
          A(rp(i)) = A(rp(i))+norm_b
          enddo
          return
      end subroutine rescale_back
