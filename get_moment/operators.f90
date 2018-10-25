! Created=Thu 20 Sep 2018 05:36:14 PM DST
! Last Modified=Fri 21 Sep 2018 11:39:53 PM DST
      ! this is a subroutine that applies operators
      ! to a wave vector

      subroutine op_commutator(N,NNZa,Aa,rpa,cola,NNZb,Ab,rpb,colb,&
              v_in,v_out)
          integer*8,intent(in)::N,NNZa,NNZb
          integer*8,dimension(NNZ),intent(in)::cola,colb
          integer*8,dimension(N+1),intent(in)::rpa,rpb
          complex*16,dimension(NNZ),intent(in)::Aa,Ab
          complex*16,dimension(N),intent(in)::v_in
          complex*16,dimension(N),intent(out)::v_out
          complex*16,dimension(N)::v_tmp1,v_tmp2,v_tmp3,v_tmp4

          call CSRmultVc16(N,NNZb,Ab,rpb,colb,v_in,v_tmp1)
          call CSRmultVc16(N,NNZa,Aa,rpa,cola,v_tmp1,v_tmp2)


          call CSRmultVc16(N,NNZa,Aa,rpa,cola,v_in,v_tmp3)
          call CSRmultVc16(N,NNZb,Ab,rpb,colb,v_tmp3,v_tmp4)

          v_out = v_tmp2 - v_tmp4
 
        ! v_out = [a,b]v_in
      end subroutine op_commutator

      subroutine op_chebyshev(N,NNZ,A,rp,col,v_pp,v_p,v)
! Warning: this works from T2. 
! If want to start from T1,  set T-1 as x
          integer*8,intent(in)::N,NNZ
          integer*8,dimension(NNZ),intent(in)::col
          integer*8,dimension(N+1),intent(in)::rp
          complex*16,dimension(NNZ),intent(in)::A
          complex*16,dimension(N),intent(in)::v_pp,v_p
          complex*16,dimension(N),intent(out)::v
          complex*16,dimension(N)::v_tmp

        call CSRmultVc16(N,NNZ,A,rp,col,v_p,v_tmp)
        v = 2d0*v_tmp - v_pp

      end subroutine op_chebyshev

