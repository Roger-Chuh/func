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