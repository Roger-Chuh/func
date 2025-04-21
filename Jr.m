function J = Jr(phi)
EPSILON = 1e-10;

J = eye(3);

phi_norm2 = sum(phi.^2);
phi_hat = SkewSymMat(phi);
phi_hat2 = phi_hat * phi_hat;

if (phi_norm2 > EPSILON)
    phi_norm = sqrt(phi_norm2);
    phi_norm3 = phi_norm2 * phi_norm;
    
    J = J - phi_hat * (1 - cos(phi_norm)) / phi_norm2;
    J = J + phi_hat2 * (phi_norm - sin(phi_norm)) / phi_norm3;
else
    % sin and cos Taylor expansion around 0
    J = J - phi_hat / 2;
    J = J + phi_hat2 / 6;
end
end