! This is the 2rd derivative to be inserted to get_rho
       rhoi= (3*norm_b**2*(gJ(0)*mu(0) + &
       Sum(2*ChebyshevT(j,-(norm_b/norm_a)) &
       * gJ(j)*mu(j), List(j,1,-1 + ForceNc))))/(norm_a**5 &
       * (1 - norm_b**2/norm_a**2)**2.5*Pi) &
       +  (gJ(0)*mu(0) + Sum(2*ChebyshevT(j,-(norm_b/norm_a))&
       * gJ(j)*mu(j), List(j,1,-1 + ForceNc)))/(norm_a**3&
       * (1 - norm_b**2/norm_a**2)**1.5*Pi) &
       -  (2*norm_b*Sum((2*j*ChebyshevU(-1 + j,&
       - (norm_b/norm_a))*gJ(j)*mu(j))/norm_a, &
       List(j,1,-1 + ForceNc)))/(norm_a**3*(1 - norm_b**2&
       / norm_a**2)**1.5*Pi) + Sum((2*j*(-(j*&
       ChebyshevU(-2 + j,-(norm_b/norm_a))) -&
        (norm_b*(-1 + j)*ChebyshevU(-1 + j,&
        - (norm_b/norm_a)))/norm_a)*gJ(nnn)*mu(nnn)) &
        / (norm_a**2*(-1 + norm_b**2/norm_a**2)),&
        List(nnn,1,-1 + ForceNc))/ &
        (norm_a*Sqrt(1 - norm_b**2/norm_a**2)*Pi)
