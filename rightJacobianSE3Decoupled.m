function J = rightJacobianSE3Decoupled(phi)
J = zeros(6, 6);
J(1:3,1:3) = Jr(phi(1:3));
J(4:6,4:6) = rodrigues(-phi(1:3));
end