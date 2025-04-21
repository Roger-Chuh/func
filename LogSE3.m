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