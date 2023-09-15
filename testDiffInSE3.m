function testDiffInSE3()

pose1 = [rodrigues([0.1 0.2 0.3]) [0.1 0.2 0.3]';0 0 0 1];
pose2 = [rodrigues([0.1 0.2 -0.3]) [0.3 -0.2 0.3]';0 0 0 1];

[dr, dp] = DiffInSE3(pose1, pose2);

pose11 = [rodrigues(dr) dp;0 0 0 1]*pose2;

end
function [dr, dp] = DiffInSE3(T1, T2)


  R = T1(1:3,1:3) * T2(1:3,1:3)';
  t = T1(1:3,4) - R * T2(1:3,4);
  dr = rodrigues(R);
  dp = JrInv(-dr) * t;


end
function J = JrInv(phi) 
  EPSILON = 1e-6;
  EPSILONSQRT = sqrt(EPSILON);

  J = eye(3);

  phi_norm2 = sum(phi.^2);
  phi_hat = SkewSymMat(phi);
  phi_hat2 = phi_hat * phi_hat;

  J = J + phi_hat / 2;
  if (phi_norm2 > EPSILON) 
    phi_norm = sqrt(phi_norm2);

    

    if (phi_norm < 3.141592653 - EPSILONSQRT) 
      
      J = J + phi_hat2 * (1 / phi_norm2 - (1 + cos(phi_norm)) / (2 * phi_norm * sin(phi_norm)));
     else 
      
      J = J + phi_hat2 / (3.141592653 * 3.141592653);
    end
    
   else 
    
    J = J + phi_hat2 / 12;
  
  end
  
end