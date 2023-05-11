function [theta,d_func_d_theta1] = solveTheta( k1,  k2,  k3,  k4,   r_theta ,d_func_d_theta, ITER)


theta = r_theta;
for  i = ITER:-1:1
    theta2 = theta .* theta;
    
    func = k4 .* theta2;
    func = func + k3;
    func = func .* theta2;
    func = func + k2;
    func = func .* theta2;
    func = func + k1;
    func = func .* theta2;
    func = func + 1.0;
    func = func .* theta;
    
    d_func_d_theta = (9) .* k4 .* theta2;
    d_func_d_theta = d_func_d_theta + (7) .* k3;
    d_func_d_theta = d_func_d_theta .* theta2;
    d_func_d_theta = d_func_d_theta + (5) .* k2;
    d_func_d_theta = d_func_d_theta .* theta2;
    d_func_d_theta = d_func_d_theta + (3) .* k1;
    d_func_d_theta = d_func_d_theta .* theta2;
    d_func_d_theta = d_func_d_theta + (1);
    
    
    %theta = theta + (r_theta - func) / d_func_d_theta;
    theta_fix = (r_theta - func) ./ d_func_d_theta;
    theta = theta+ theta_fix;
    if 0 % (abs(theta_fix) < 1e-8)
        break;
    end
end
d_func_d_theta1 = d_func_d_theta;
end