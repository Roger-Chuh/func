function testZNCC()

img0 = imread('G:\matlab\data\direct\gt\D2_011\4\cam2\508373517000.png');
img1 = imresize(img0, 0.5, 'method', 'bicubic');
img2 = imresize(img0, 0.25, 'method', 'bicubic');

imgs{1,1} = img0;
imgs{2,1} = img1;
imgs{3,1} = img2;

figure, subplot(1,3,1);imshow(img0);
subplot(1,3,2);imshow(img1);
subplot(1,3,3);imshow(img2);



patch_offset = [0, 0;
    -2 , -2;
    -1 , -1;
    -2 , 2 ;
    -1 , 1 ;
    1 , 0 ;
    2 , 0 ;
    3 , 0 ];

patch_offset = [0, 0;
    0 , -2;
    -1 , -1;
    -2 , 0 ;
    -1 , 1 ;
    0 , 2 ;
    2 , 0 ;
    1 , -1;
    1, 1];



PATCH_SIZE = size(patch_offset,1);
Mat_ZNSSD_I = eye(PATCH_SIZE);
vec_1 = ones(PATCH_SIZE,1);
J_ZNSSD_mean = Mat_ZNSSD_I - (vec_1 / PATCH_SIZE) * vec_1';
param_host = [234.3060155, 234.6350819, 319.9609246, 240.7165617, 0.2206591747,-0.1801930843, 0.07274986781, -0.01355495556];



Tcw = [rodrigues([0.01 0.02 -0.03]) [10 20 -1000]'; 0 0 0 1];
idp = 1/3000;
uv = [300 200];
uv_patch = uv + patch_offset;
[host_bearing, d_pt3d_d_uv] = unprojectKB8(uv(:,1), uv(:,2), param_host);
[host_uv, d_uv_d_pt3d] = projectKB8(host_bearing(:,1), host_bearing(:,2), host_bearing(:,3), param_host);

xyz = host_bearing./idp;
xyz_target = (Tcw(1:3,1:3) * xyz' + Tcw(1:3,4))';

[target_uv, d_uv_d_pt3d_target] = projectKB8(xyz_target(:,1), xyz_target(:,2), xyz_target(:,3), param_host);

for i = 1 : length(imgs)
    host_img = double(imgs{i,1});
    gx = zeros(size(host_img));
    gy = gx;
    
    width = size(host_img,2);
    height = size(host_img,1);
    
    gx(:,2 : width-1) = 0.5 * (host_img(:,3:width) - host_img(:, 1:width-2));
    gy(2 : height-1, :) = 0.5 * (host_img(3:height,:) - host_img(1:height-2,:));
    
    [xMat, yMat] = meshgrid(1:width, 1:height);
    pix = [xMat(:) yMat(:)];
    pix_tensor = repmat(pix, 1,1,size(patch_offset,1));
    patch_offset_tensor = repmat(patch_offset(:), width * height, 1);
    patch_offset_tensor = reshape(patch_offset_tensor,size(patch_offset, 1), 2, []);
    patch_offset_tensor = permute(patch_offset_tensor, [3 2 1]);
    
    pix_patch_tensor = pix_tensor + patch_offset_tensor;
    x_corrd = squeeze(pix_patch_tensor(:,1,:));
    y_corrd = squeeze(pix_patch_tensor(:,2,:));
    
    valid_ind = round(x_corrd(:)) >= 1 & round(y_corrd(:)) >= 1 & round(x_corrd(:)) <= width & round(y_corrd(:)) <= height;
    %     pix_patch_tensor_vec = reshape(pix__patch_tensor, )
    pix_vec = [x_corrd(:) y_corrd(:)];
    ind = sub2ind([height width], round(pix_vec(valid_ind,2)), round(pix_vec(valid_ind,1)));
    values = host_img(ind);
    value_vec = nan(length(valid_ind), 1);
    value_vec(valid_ind) = values;
    value_tensor = reshape(value_vec, height, width, []);
    if 0
        pix_check = [134 334];figure,imshow(host_img,[]);hold on;plot(pix_check(:,1) + patch_offset(:,1), pix_check(:,2) + patch_offset(:,2),'.r')
        squeeze(value_tensor(pix_check(2),pix_check(1),:))'
    end
    value_mean = sum(value_tensor,3)./size(patch_offset,1);
    value_mean_tensor = repmat(value_mean, 1, 1, size(patch_offset,1));
    value_tensor_sub_mean = value_tensor - value_mean_tensor;
    value_tensor_sub_mean_reshape = reshape(value_tensor_sub_mean, [], size(patch_offset,1));
    %     value_tensor_sub_mean_norm = norm(value_tensor_sub_mean_reshape');
    [~,dist] = NormalizeVector(value_tensor_sub_mean_reshape);
    dist_map = reshape(dist, height, width);
    norm_map_tensor = repmat(dist_map, 1, 1, size(patch_offset,1));
    value_tensor_sub_mean_normalized = value_tensor_sub_mean./norm_map_tensor;
    
    if 0
        pix_check = [334 234];
        b = squeeze(value_tensor_sub_mean_normalized(pix_check(2), pix_check(1),:))'
        norm(b)
        a = squeeze(value_tensor(pix_check(2), pix_check(1),:))'
        aa = (a - mean(a))./norm((a - mean(a)));
        err = aa - b
    end
    
    %     value_tensor_sub_mean_normalized_dot = dot(value_tensor_sub_mean_normalized);
    
    %     J_ZNSSD_J_I = (Mat_ZNSSD_I - (normalized_vals.matrix() * normalized_vals.matrix().transpose())) / sigma * J_ZNSSD_mean;
    
    eig_map = zeros(height, width);
    eig_map_pose = zeros(height, width);
    eig_map_pose_min = zeros(height, width);
    eig_map_pose_max = zeros(height, width);
    for row = 1 : height
        for col = 1 : width
            normalized_vals = squeeze(value_tensor_sub_mean_normalized(row, col,:))';
            sigma0 = dist_map(row, col);
            if isnan(sum(normalized_vals))
                continue;
            end
            sigma = dist_map(row, col);
            J_ZNSSD_J_I = sigma * (Mat_ZNSSD_I - (normalized_vals * normalized_vals')) / sigma * J_ZNSSD_mean;
            
            coords = [col row] + patch_offset;
            valid = round(coords(:,1)) >= 1 & round(coords(:,2)) >= 1 & round(coords(:,1)) <= width & round(coords(:,2)) <= height;
            if sum(valid) < PATCH_SIZE
                continue;
            end
            ind = sub2ind([height width], round(coords(:,2)), round(coords(:,1)));
            grad = [gx(ind) gy(ind)];
            J_ZNSSD_J_uv = J_ZNSSD_J_I * grad;
            
            [~,dist] = NormalizeVector(J_ZNSSD_J_uv);
            g20 = dist.^2;
            g2 = 50.0 ./ (50.0 + g20);
            gradWeightMat = diag(g2);
            gradWeightMat = diag(ones(length(g2),1));
            
            J_ZNSSD_J_dir = J_ZNSSD_J_uv * d_uv_d_pt3d;
            H_uv = J_ZNSSD_J_uv' * gradWeightMat *  J_ZNSSD_J_uv;
            H_dir = J_ZNSSD_J_dir' * gradWeightMat * J_ZNSSD_J_dir;
            
            
            n = host_bearing';
            dp_dx0 = zeros(3,6);
            dp_dx0(:,1:3) = (-SkewSymMat(n) * Tcw(1:3,1:3) + idp * SkewSymMat(Tcw(1:3,4)) * Tcw(1:3, 1:3));
            dp_dx0(:,4:6) = idp * Tcw(1:3, 1:3);
            H_x0 = dp_dx0' * H_dir * dp_dx0;
            J_x0 = J_ZNSSD_J_dir * (dp_dx0);
            
            [V,D] = eig(H_uv);
            [u,s,v] = svd(H_uv);
            eig_norm = norm([D(1,1) D(2,2)]);
            eig_map(row, col) = eig_norm;
            
            [V2,D2] = eig(H_x0);
            [u2,s2,v2] = svd(H_x0);
            eig_map_pose(row, col) = norm(diag(D2));
            eig_map_pose_min(row, col) = min(diag(D2));
            eig_map_pose_max(row, col) = max(diag(D2));
        end
    end
    mask = eig_map < 0.1 & eig_map > 0.05;
    mask_pose = eig_map_pose < 10000 & eig_map_pose > 1000;
    
    
    mask = eig_map > 50; %< 0.1 & eig_map > 0.05;
    mask_pose = eig_map_pose > 5e6; % < 10000 & eig_map_pose > 1000;
    
%     figure,subplot(1,2,1);imshowpair(host_img, mask);title('H uv');subplot(1,2,2);imshowpair(host_img, mask_pose);title('H x0');
    figure,subplot(2,2,1);imshowpair(host_img, eig_map);title('H uv');subplot(2,2,2);imshowpair(host_img, eig_map_pose);title('H x0');
    subplot(2,2,3);imshow(eig_map, []);title('H uv');subplot(2,2,4);imshow(eig_map_pose, []);title('H x0');
    Eig_Map{i,1} = eig_map;
    Eig_Map_pose{i,1} = eig_map_pose;
%     Eig_Map_pose{i,1} = eig_map_pose_max;
%     Eig_Map_pose{i,2} = eig_map_pose_min;
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