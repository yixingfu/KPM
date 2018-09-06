! This is the 2rd derivative to be inserted to get_rho
       rhoi= (3*norm_b**2*(gJ(0)*mu_avg(0) + &
       Sum(2*ChebyT(j,-(norm_b/norm_a)) &
       * gJ(j)*mu_avg(j), List(j,1,-1 + ForceNc))))/(norm_a**5 &
       * (1 - norm_b**2/norm_a**2)**2.5*Pi) &
       +  (gJ(0)*mu_avg(0) + Sum(2*ChebyT(j,-(norm_b/norm_a))&
       * gJ(j)*mu_avg(j), List(j,1,-1 + ForceNc)))/(norm_a**3&
       * (1 - norm_b**2/norm_a**2)**1.5*Pi) &
       -  (2*norm_b*Sum((2*j*ChebyU(-1 + j,&
       - (norm_b/norm_a))*gJ(j)*mu_avg(j))/norm_a, &
       List(j,1,-1 + ForceNc)))/(norm_a**3*(1 - norm_b**2&
       / norm_a**2)**1.5*Pi) + Sum((2*j*(-(j*&
       ChebyU(-2 + j,-(norm_b/norm_a))) -&
        (norm_b*(-1 + j)*ChebyU(-1 + j,&
        - (norm_b/norm_a)))/norm_a)*gJ(j)*mu_avg(j)) &
        / (norm_a**2*(-1 + norm_b**2/norm_a**2)),&
        List(j,1,-1 + ForceNc))/ &
        (norm_a*Sqrt(1 - norm_b**2/norm_a**2)*Pi)
