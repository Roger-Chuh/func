function testDiffInSE3()

pose_est = [rodrigues([0.1 0.2 0.3]) [0.1 0.2 0.3]';0 0 0 1];
pose_prior = [rodrigues([0.1 0.2 -0.3]) [0.3 -0.2 0.3]';0 0 0 1];
% pose_est = [rodrigues([-0.0001 0.0002 -0.0003]) [0.0001 0.00002 0.00003]';0 0 0 1] * pose_prior;
[dr, dp] = DiffInSE3(pose_est, pose_prior);


delta_T =  Exp(dr, dp);

pose_est_correct = delta_T*pose_prior;

pose_est2 = [rodrigues(dr) dp;0 0 0 1]*pose_prior;

lie_est = [rodrigues(pose_prior(1:3,1:3)); pose_prior(1:3,4)] + [dr; dp];
pose_est3 = [rodrigues(lie_est(1:3)) [lie_est(4:6)];0 0 0 1] - pose_est


pose_est2 - pose_est

pose_est_correct - pose_est


end
function [dr, dp] = DiffInSE3(T1, T2)


R = T1(1:3,1:3) * T2(1:3,1:3)';
t = T1(1:3,4) - R * T2(1:3,4);


T1_check = [R t;0 0 0 1]*T2;

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
function drdp = LogSE3(T)
drdp = zeros(6,1);
R = T(1:3,1:3);
omega = rodrigues(R);
theta = norm(omega);
drdp(1:3) = omega;
Omega = SkewSymMat(omega);
if (theta < 1e-10)
    V_inv = eye(3) - 0.5 * Omega + (1. / 12.) * (Omega * Omega);
    drdp(4:6) = V_inv * T(1:3,4);
else
    half_theta = (0.5) * theta;
    V_inv = (eye(3) - (0.5) * Omega + ((1) - theta * cos(half_theta) / ((2) * sin(half_theta))) / (theta * theta) * (Omega * Omega));
    drdp(4:6) = V_inv * T(1:3,4);
end
end
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