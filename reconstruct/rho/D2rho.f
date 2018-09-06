          FUNCTION D2rho(n,a,b,x)
          INTEGER :: n
          Integer, parameter :: idp = kind(1.d0)
          REAL(idp) :: a,b,x,D2rho,ChebyU,ChebyT

          D2rho=((a**2 + 2*(b - x)**2)*ChebyT(n,(-b + x)/a))/&
      (a**2 - (b - x)**2)**2.5 + &
      (2*n*(-b + x)*ChebyU(-1 + n,(-b + x)/a))/&
      (a*(a**2 - (b - x)**2)**1.5) + &
      (n*(a*n*ChebyU(-2 + n,(-b + x)/a) + &
      (-1 + n)*(b - x)*ChebyU(-1 + n,(-b + x)/a)))/&
      (a*(a**2 - (b - x)**2)**1.5)
          RETURN
          END FUNCTION D2rho



          FUNCTION D2rhomu0(a,b,x)
              Integer, parameter :: idp = kind(1.d0)
              REAL(idp) :: a,b,x,D2rhomu0

              D2rhomu0=(3*(-b + x)**2)/(a**2-(-b + x)**2)**2.5&
     +(a**2 - (-b + x)**2)**(-1.5)
              RETURN
          END FUNCTION D2rhomu0
