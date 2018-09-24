!Created=Wed 13 Dec 2017 07:12:25 PM STD
!Last Modified=Sun 23 Sep 2018 05:07:19 PM DST
      real*8 function ChebyT(n,x)
          real*8::x
          integer::n
          ChebyT = dcos(n*dacos(x))
          return
      end function ChebyT

      real*8 function ChebyU(n,x)
          real*8::x
          integer::n
          ChebyU = dsin((n+1)*dacos(x))/dsin(dacos(x))
          return
      end function ChebyU


