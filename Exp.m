function result =  Exp(w, v)


omega = w;
theta_sq = sum(omega.^2);

if (theta_sq < 1e-10)
    theta = 0;
else
    theta = sqrt(theta_sq);
end
so3 = rodrigues(omega);
Omega = SkewSymMat(omega);
Omega_sq = Omega * Omega;
V = zeros(3, 3);

if (theta < 1e-5)
    V = so3;
    
else
    theta_sq = theta * theta;
    V = (eye(3) + ((1) - cos(theta)) / (theta_sq)*Omega + (theta - sin(theta)) / (theta_sq * theta) * Omega_sq);
end

tran = V * v;
result = eye(4);

result(1:3,1:3) = so3;
result(1:3,4) = tran;

end