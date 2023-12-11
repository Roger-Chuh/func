function testBestSearchLevel()

close all;
if 0
    img0 = imread('G:\matlab\data\direct\gt\D2_011\4\cam2\508373517000.png');
    img1 = imresize(img0, 0.5, 'method', 'bicubic');
    img2 = imresize(img0, 0.25, 'method', 'bicubic');
    
    
    
    figure, subplot(1,3,1);imshow(img0);
    subplot(1,3,2);imshow(img1);
    subplot(1,3,3);imshow(img2);
end


param_host = [234.3060155, 234.6350819, 319.9609246, 240.7165617, 0.2206591747,-0.1801930843, 0.07274986781, -0.01355495556];
param_target = [236.4172329, 236.8149768, 322.7260237, 238.7688593, 0.2108834396, -0.1630992223, 0.05792931434, -0.008858528184];


for (i = 1 : 3)
    host_params{i,1} = param_host;
    host_params{i,1}(1:2) = param_host(1:2)*2^-(i-1);
    host_params{i,1}(3:4) = (param_host(3:4) + 0.5)*2^-(i-1) - 0.5;
    
    target_params{i,1} = param_target;
    target_params{i,1}(1:2) = param_target(1:2)*2^-(i-1);
    target_params{i,1}(3:4) = (param_target(3:4) + 0.5)*2^-(i-1) - 0.5;
end


offset = [0, 0;
    -2 , -2;
    -1 , -1;
    -2 , 2 ;
    -1 , 1 ;
    1 , 0 ;
    2 , 0 ;
    3 , 0 ];
patches = [200 100] + offset;

[host_bearings] = unprojectKB8(patches(:,1), patches(:,2), host_params{2,1});
[target_patches_tmp] = projectKB8(host_bearings(:,1), host_bearings(:,2),host_bearings(:,3), host_params{2,1});

depth = 3000;

xyz = host_bearings.*depth;

T_th = [rodrigues([0.1 0.2 0.3]) [0 0 -1800]'; 0 0 0 1];

for i = 1 : size(offset,1)
    xyz_target(i,:) = (T_th(1:3,1:3) * xyz(i,:)' + T_th(1:3,4))';
    
end

figure,imshow(ones(480, 640));hold on;plot(patches(:,1), patches(:,2),'.r');plot(patches(1:2,1), patches(1:2,2),'or');
for i = 1 : 3
    [target_patches{i,1}] = projectKB8(xyz_target(:,1), xyz_target(:,2),xyz_target(:,3), target_params{i,1});
    plot(target_patches{i,1}(:,1), target_patches{i,1}(:,2),'.b');
    plot(target_patches{i,1}(1:2,1), target_patches{i,1}(1:2,2),'ob');
end
[target_patches{1,1} ((target_patches{2,1} + 0.5).*2^1 - 0.5) ((target_patches{3,1} + 0.5).*2^2 -0.5)]


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


if (length(X) > 1)
    d_uv_d_pt3d = [];
    d_uv_d_param = [];
    return;
end

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
    theta = theta +  theta_fix;
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

mx = (u - cx) ./ fx;
my = (v - cy) ./ fy;

theta = 0;
sin_theta = 0;
cos_theta = 1;
thetad = sqrt(mx .* mx + my .* my);

thetad0 = min(max(thetad,-3.141592653/2),3.141592653/2);
scaling = 1.0;
d_func_d_theta = 0;
if (thetad > 1e-18)
    [theta,d_func_d_theta1] = solveTheta(k1, k2, k3, k4, thetad, d_func_d_theta,10);
    d_func_d_theta = d_func_d_theta1;
    sin_theta = sin(theta);
    cos_theta = cos(theta);
    scaling = sin_theta ./ thetad;
    
    
    d_thetad_d_mx = mx ./ thetad;
    d_thetad_d_my = my ./ thetad;
    theta2 = theta .* theta;
    d_scaling_d_thetad = (thetad .* cos_theta ./ d_func_d_theta - sin_theta) ./ (thetad .* thetad);
    d_cos_d_thetad = -sin_theta ./ d_func_d_theta;
    d_scaling_d_k1 = -cos_theta .* theta .* theta2 ./ (d_func_d_theta .* thetad);
    d_cos_d_k1 = -d_cos_d_thetad .* theta .* theta2;
end

pt3d = [mx .* scaling my .* scaling cos_theta .* ones(length(u), 1)];

if (length(u) > 1)
    [~,dist] = NormalizeVector(pt3d);
    
    d_pt3d_d_uv = [];
    d_pt3d_d_param = [];
    return;
end


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