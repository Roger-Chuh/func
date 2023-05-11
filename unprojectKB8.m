function p3d = unprojectKB8( fx,  fy,  cx,  cy,  k1,  k2,  k3,  k4, proj)


mx = (proj(:,1) - cx) ./ fx;
my = (proj(:,2) - cy) ./ fy;




theta = 0;
sin_theta = 0;
cos_theta = 1;
thetad = sqrt(mx .* mx + my .* my);

thetad = min(max(thetad,-3.141592653/2),3.141592653/2);
scaling = 1.0;
d_func_d_theta = 0;
if (thetad > 1e-8)
    [theta,d_func_d_theta1] = solveTheta(k1, k2, k3, k4, thetad, d_func_d_theta,3);
    d_func_d_theta = d_func_d_theta1;
    sin_theta = sin(theta);
    cos_theta = cos(theta);
    scaling = sin_theta ./ thetad;
end

p3d = [mx .* scaling my .* scaling cos_theta];


end