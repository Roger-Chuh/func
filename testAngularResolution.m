function testAngularResolution()
global Param
K = eye(3);
param = [236.8820281-0;
236.7437341-0;
325.0065396;
241.9557923;
0.2124214875;
-0.1719834148;
0.06211290885;
-0.009899091649];

Param = param;



ang_mat_x = zeros(480, 640);
ang_mat_y = zeros(480, 640);
for i = 1 : 639
    for j = 1 : 479
        
        [pt3d, d_pt3d_d_uv, d_pt3d_d_param] = unprojectKB8(i, j, param);
        % horizontal
        [pt3d_x, d_pt3d_d_uv, d_pt3d_d_param] = unprojectKB8(i+1, j, param);
        ang_mat_x(j,i) = CalcDegree(pt3d, pt3d_x);
        % horizontal
        [pt3d_y, d_pt3d_d_uv, d_pt3d_d_param] = unprojectKB8(i, j+1, param);
        ang_mat_y(j,i) = CalcDegree(pt3d, pt3d_y);
        
    end
    
end


end
function [pt2d, d_uv_d_pt3d, d_uv_d_param] = projectKB8(X, Y, Z, param)
fx = param(1);
fy = param(2);
cx = param(3);
cy = param(4);
k1 = param(5);
k2 = param(6);
k3 = param(7);
k4 = param(8);

x = X./Z;
y = Y./Z;
z = 1;%Z;


r2 = x .* x + y .* y;
r = sqrt(r2);

R2 = X .* X + Y .* Y;
R = sqrt(R2);


theta1 = atan2(r, z);
theta = atan(r);
theta2 = theta .* theta;
theta4 = theta2 .* theta2;
theta6 = theta4 .* theta2;
theta8 = theta6 .* theta2;

r_theta = theta .* (1 + k1 .* theta2 + k2 .* theta4 + k3 .* theta6 + k4 .* theta8);

if r > 1e-18
    norm_inv = 1./r;
else
    norm_inv = ones(length(r), 1);
end
% norm_inv = r > 1e-8 ? double(1.0) / r : 1;

mx = r_theta .* x .* norm_inv;
my = r_theta .* y .* norm_inv;

x_c = X;
y_c = Y;
z_c = Z;


NORM_inv = 1/R;
cosphi = x_c * NORM_inv;
sinphi = y_c * NORM_inv;


pt2d(:,1) = fx .* mx + cx;
pt2d(:,2) = fy .* my + cy;

d_uv_d_fxfycxcy = [ -mx     0   -1    0 ;
                      0   -my    0   -1   ]; 
d_uv_d_k1 = [-fx * theta2 * theta * x .* norm_inv
             -fy * theta2 * theta * y .* norm_inv];
d_uv_d_k2 = d_uv_d_k1.*theta2;
d_uv_d_k3 = d_uv_d_k2.*theta2;
d_uv_d_k4 = d_uv_d_k3.*theta2;

d_uv_d_param = -[d_uv_d_fxfycxcy d_uv_d_k1 d_uv_d_k2 d_uv_d_k3 d_uv_d_k4];

d_uv_d_mxmy = [fx 0;0 fy];

a = x;
b = y;
thd = r_theta;
d_thd_d_th = 1 + 3*k1.*theta2 + 5*k2.*theta4 + 7*k3.*theta6 + 9*k4.*theta8;  %% + 11*k5.*th10 + 13*k6.*th12;


if 1
    dRtheta_dtheta = d_thd_d_th;
    coff1 = (dRtheta_dtheta * cosphi * z_c) / (R * (R * R + z_c * z_c));
    coff23 = R * R * R;
    xd = r_theta * cosphi; % mx
    yd = r_theta * sinphi; % my
    dp_dxc = zeros(2,3);
    dp_dxc(1, 1) = coff1 * x_c + (r_theta * y_c * y_c) / coff23;
    dp_dxc(1, 2) = coff1 * y_c - (r_theta * x_c * y_c) / coff23;
    dp_dxc(1, 3) = -1.0 * (dRtheta_dtheta * cosphi * R) / (R * R + z_c * z_c);
    
    coff2 = (dRtheta_dtheta * sinphi * z_c) / (R * (R * R + z_c * z_c));
    dp_dxc(2, 1) = coff2 * x_c - (r_theta * x_c * y_c) / coff23;
    dp_dxc(2, 2) = coff2 * y_c + (r_theta * x_c * x_c) / coff23;
    dp_dxc(2, 3) = -1.0 * (dRtheta_dtheta * sinphi * R) / (R * R + z_c * z_c);
end


d_x_r_d_a = b.^2./r.^3.*thd + a.^2./(r.^2 + r.^4)*d_thd_d_th;
d_x_r_d_b = -a.*b./r.^3.*thd + a.*b./(r.^2 + r.^4)*d_thd_d_th;
d_y_r_d_a = -a.*b./r.^3.*thd + a.*b./(r.^2 + r.^4)*d_thd_d_th;
d_y_r_d_b = a.^2./r.^3.*thd + b.^2./(r.^2 + r.^4)*d_thd_d_th;
d_mxmy_d_xy = [d_x_r_d_a d_x_r_d_b; d_y_r_d_a d_y_r_d_b];

d_xy_d_XYZ = [ 1./Z    0     -X./Z.^2;
                 0    1./Z   -Y./Z.^2];


d_uv_d_pt3d = d_uv_d_mxmy * d_mxmy_d_xy * d_xy_d_XYZ;

end
function [theta,d_func_d_theta1] = solveTheta( k1,  k2,  k3,  k4,   r_theta ,d_func_d_theta, ITER)


theta = r_theta;
for  i = ITER:-1:0
    theta2 = theta * theta;
    
    func = k4 * theta2;
    func = func + k3;
    func = func * theta2;
    func = func + k2;
    func = func * theta2;
    func = func + k1;
    func = func * theta2;
    func = func + 1.0;
    func = func * theta;
    
    d_func_d_theta = (9) * k4 * theta2;
    d_func_d_theta = d_func_d_theta + (7) * k3;
    d_func_d_theta = d_func_d_theta * theta2;
    d_func_d_theta = d_func_d_theta + (5) * k2;
    d_func_d_theta = d_func_d_theta * theta2;
    d_func_d_theta = d_func_d_theta + (3) * k1;
    d_func_d_theta = d_func_d_theta * theta2;
    d_func_d_theta = d_func_d_theta + (1);
    
    
    %theta = theta + (r_theta - func) / d_func_d_theta;
    theta_fix = (r_theta - func) / d_func_d_theta;
    theta = theta+  theta_fix;
if 0
    if (abs(theta_fix) < 1e-8)
        break;
    end
end
end
d_func_d_theta1 = d_func_d_theta;
end
function [pt3d, d_pt3d_d_uv, d_pt3d_d_param] = unprojectKB8(u, v, param)
fx = param(1);
fy = param(2);
cx = param(3);
cy = param(4);
k1 = param(5);
k2 = param(6);
k3 = param(7);
k4 = param(8);

mx = (u - cx) / fx;
my = (v - cy) / fy;

theta = 0;
sin_theta = 0;
cos_theta = 1;
thetad = sqrt(mx * mx + my * my);

thetad = min(max(thetad,-3.141592653/2),3.141592653/2);
scaling = 1.0;
d_func_d_theta = 0;
if (thetad > 1e-18)
    [theta,d_func_d_theta1] = solveTheta(k1, k2, k3, k4, thetad, d_func_d_theta,10);
    d_func_d_theta = d_func_d_theta1;
    sin_theta = sin(theta);
    cos_theta = cos(theta);
    scaling = sin_theta / thetad;
    
    
    d_thetad_d_mx = mx / thetad;
    d_thetad_d_my = my / thetad;
    theta2 = theta * theta;
    d_scaling_d_thetad = (thetad * cos_theta / d_func_d_theta - sin_theta) / (thetad * thetad);
    d_cos_d_thetad = -sin_theta / d_func_d_theta;
    d_scaling_d_k1 = -cos_theta * theta * theta2 / (d_func_d_theta * thetad);
    d_cos_d_k1 = -d_cos_d_thetad * theta * theta2;
end

pt3d = [mx * scaling my * scaling cos_theta];
d_res0_d_mx = scaling + mx * d_scaling_d_thetad * d_thetad_d_mx;
d_res0_d_my = mx * d_scaling_d_thetad * d_thetad_d_my;

d_res1_d_mx = my * d_scaling_d_thetad * d_thetad_d_mx;
d_res1_d_my = scaling + my * d_scaling_d_thetad * d_thetad_d_my;

d_res2_d_mx = d_cos_d_thetad * d_thetad_d_mx;
d_res2_d_my = d_cos_d_thetad * d_thetad_d_my;


c0(1,:) = d_res0_d_mx / fx; 
c0(2,:) = d_res1_d_mx / fx;
c0(3,:) = d_res2_d_mx / fx;

c1(1,:) = d_res0_d_my / fy; 
c1(2,:) = d_res1_d_my / fy;
c1(3,:) = d_res2_d_my / fy;

d_pt3d_d_uv = [c0 c1];

d_pt3d_d_k1 = [mx * d_scaling_d_k1; my * d_scaling_d_k1; d_cos_d_k1];
d_pt3d_d_k2 = d_pt3d_d_k1*theta2;
d_pt3d_d_k3 = d_pt3d_d_k2*theta2;
d_pt3d_d_k4 = d_pt3d_d_k3*theta2;

d_pt3d_d_param = [-c0*mx -c1*my -c0 -c1 d_pt3d_d_k1 d_pt3d_d_k2 d_pt3d_d_k3 d_pt3d_d_k4];


end