function LookupTable()
global param CENTER_FOV;


one_based = 1; 0; 1;

% close all;
param = [236.8820281-0;
    236.7437341-0;
    325.0065396;
    241.9557923;
    0.2124214875;
    -0.1719834148;
    0.06211290885;
    -0.009899091649];
param = [236.829010 237.239883 321.422058 236.294388 0.217803 -0.175568 0.064648 -0.010765]';

param = [234.3060155, 234.6350819, 319.9609246, 240.7165617, 0.2206591747, -0.1801930843, 0.07274986781, -0.01355495556]';

CENTER_FOV = 90; 170;150;  120; 90; 60; 30; 80;90;
fov = 3 * CENTER_FOV;
K = [param(1) 0 param(3); 0 param(2) param(4);0 0 1];
virtual_img_width = 64; 100; 512; 32;40; 100; 500;  250; 100; 250;500; 250 ;

width = 640;
height = 480;

[rotations, f_center,intrMatPinhole] = generateTransforms(virtual_img_width, fov);

rotations_bak = rotations;
rotations{2,1} = rotations_bak{5,1};
rotations{5,1} = rotations_bak{2,1};
if one_based
    [xMat_fish, yMat_fish] = meshgrid(1:width, 1:height);
else
    [xMat_fish, yMat_fish] = meshgrid(0:width-1, 0:height-1);
end
fish = [xMat_fish(:) yMat_fish(:)];

[pt3d, d_pt3d_d_uv, d_pt3d_d_param] = unprojectKB8(fish(:,1), fish(:,2), param);
[pt2d, d_uv_d_pt3d, d_uv_d_param] = projectKB8(pt3d(:,1), pt3d(:,2), pt3d(:,3), param);
[~,err] = NormalizeVector(pt2d-fish);
valid_fisheye_area = find(err<0.000001);
mask = zeros(height, width);
mask(valid_fisheye_area) = 1;

Valid = cell(5,1);
Pix2Pinhole = cell(5,1);
valid_stack = [];
pinhole_valid_ind = cell(5,1);
color = {'ob', '.r','xg','sm','dy'};
figure(100);subplot(2,3,1),imshow(zeros(height,width)),hold on;

% figure(200),
%% 把带畸变图映射成不带畸变图的过程为：用针孔内参unproject，施加旋转矩阵，再用畸变内参project得到的坐标映射关系
%% 我知道target图的位置是什么（都是整数像素），只是不知道它的灰度值是啥，所以要去src图上取（对应的浮点像素），进而要建立从target坐标到src坐标的映射
%% 把不带畸变图映射成带畸变图的过程为：用鱼眼内参unproject，施加旋转矩阵，再用针孔内参project得到的坐标映射关系
for i = 1 : length(rotations)
    fish2pinhole1 = (intrMatPinhole * (rotations{i}' * pt3d'));
    fish2pinhole = pflat(fish2pinhole1);
    fish2pinhole = fish2pinhole(1:2,:)';
    valid = find(fish2pinhole(:,1)>1 & fish2pinhole(:,1)<virtual_img_width & ...
        fish2pinhole(:,2)>1 & fish2pinhole(:,2)<virtual_img_width & fish2pinhole1(3,:)' > 0 & ...
        err<0.000001);
    valid = setdiff(valid, valid_stack);
    valid_stack = [valid_stack; valid];
    Valid{i,1} = valid;
    Pix2Pinhole{i,1} = fish2pinhole(valid,:);
    ind = sub2ind([virtual_img_width virtual_img_width], round(fish2pinhole(valid,2)), round(fish2pinhole(valid,1)));
    pinhole_valid_ind{i,1} = ind;
    if 1
        figure(100);subplot(2,3,1);plot(fish(valid,1),fish(valid,2), color{i});
        subplot(2,3,1+i),imshow(zeros(virtual_img_width, virtual_img_width));hold on;plot(fish2pinhole(valid,1),fish2pinhole(valid,2),color{i});
        if i ~= 1
            title(sprintf('pinhole id: %d', i));
        else
            title(sprintf('pinhole id: %d\nvirtual camera width: %d pixel\nvirtual camera fov: %f degree', i, virtual_img_width, CENTER_FOV));
        end
        
    end
    
end

assert(size(cell2mat(Pix2Pinhole),1) == size(valid_stack,1));
if one_based
    [xMat_pinhole, yMat_pinhole] = meshgrid(1:virtual_img_width, 1:virtual_img_width);
else
    [xMat_pinhole, yMat_pinhole] = meshgrid(0:virtual_img_width-1, 0:virtual_img_width-1);
end
pinhole = [xMat_pinhole(:) yMat_pinhole(:)];
pinhole(:,3) = 1;
Lut = cell(5,2);
pinhole_mask = cell(5,1);
figure(200);subplot(2,3,1),imshow(zeros(height,width)),hold on;
for i = 1 : length(rotations)
    pinhole2fish1 = (rotations{i} * inv(intrMatPinhole) * pinhole')';
    pinhole2fish = projectKB8(pinhole2fish1(:,1), pinhole2fish1(:,2), pinhole2fish1(:,3), param);
    [pt3d_check] = unprojectKB8(pinhole2fish(:,1), pinhole2fish(:,2), param);
    pinhole2fish_check = projectKB8(pt3d_check(:,1), pt3d_check(:,2), pt3d_check(:,3), param);
    [~,err_check] = NormalizeVector(pinhole2fish-pinhole2fish_check);
    
    
    valid_fish1 = find(pinhole2fish(:,1)>1 & pinhole2fish(:,1)<width & ...
        pinhole2fish(:,2)>1 & pinhole2fish(:,2)<height & pinhole2fish1(:,3) > 0 & err_check < 0.00001);
    if 0
        valid_fish = intersect(valid_fish1, pinhole_valid_ind{i,1});
    else
        valid_fish = valid_fish1;
    end
    
    xMat_pinhole_ = -1.*ones(virtual_img_width, virtual_img_width);
    yMat_pinhole_ = -1.*ones(virtual_img_width, virtual_img_width);
    
    %     pinhole_mask{i,1} = xMat_pinhole_;
    
    xMat_lut = -1.*ones(virtual_img_width, virtual_img_width);
    yMat_lut = -1.*ones(virtual_img_width, virtual_img_width);
    
    xMat_pinhole_(valid_fish) = xMat_pinhole(valid_fish);
    yMat_pinhole_(valid_fish) = yMat_pinhole(valid_fish);
    xMat_lut(valid_fish) = pinhole2fish(valid_fish,1);
    yMat_lut(valid_fish) = pinhole2fish(valid_fish,2);
    
    pinhole_mask{i,1} = xMat_pinhole_;
    Lut{i,1} = xMat_lut;
    Lut{i,2} = yMat_lut;
    if 1
        mm = zeros(virtual_img_width, virtual_img_width);
        mm(pinhole_valid_ind{i,1}) = 1;
        %         figure,imshow(mm)
        mmm = zeros(virtual_img_width, virtual_img_width);
        mmm(valid_fish1) = 1;
        %         figure,imshow(mmm);
        
        figure(200);subplot(2,3,1);plot(pinhole2fish(valid_fish,1),pinhole2fish(valid_fish,2), color{i});
        subplot(2,3,1+i),imshow(zeros(virtual_img_width, virtual_img_width));hold on;plot(pinhole(valid_fish,1),pinhole(valid_fish,2),color{i});title(sprintf('pinhole id: %d', i));
        drawnow;
        
    end
end


xyz = 1000.*(rand(10000,3)-0.5);
xyz(:,3) = 1000 + xyz(:,3);



mask_test = zeros(height,width);
mask_test(valid_fisheye_area) = 1;
SE = strel('cube',10);
mask_test_use = imerode(mask_test,SE);


valid_fisheye_area_use = find(mask_test_use == 1);

xyz = pt3d(valid_fisheye_area_use,:);
xyz = xyz + 0.01.*(rand(length(valid_fisheye_area_use),3)-0.5);


step = 2; 10; 100; 1;100;
test_id = 1 : step : size(xyz,1);
valid_fisheye_area_test = valid_fisheye_area_use(test_id);
proj_stack1 = zeros(length(test_id),2);
proj_stack2 = [];
cnt1 = 1;
tic;
for i = test_id
    proj1 = projectKB8(xyz(i,1), xyz(i,2), xyz(i, 3), param);
    
    proj_stack1(cnt1,:) = proj1;
    cnt1 = cnt1+1;
end
toc;

cnt2 = 1;
cnt22 = 1;
n_rotations = length(rotations);
found_ind = [];
tic;
for i = test_id
    found = false;
    for j = 1 : n_rotations
        %     proj1 = projectKB8(xyz(i,1), xyz(i,2), xyz(i, 3), param);
        proj_pinhole = pflat(intrMatPinhole * rotations{j}' * xyz(i,:)');
        if (proj_pinhole(1) > 1 && proj_pinhole(1) < virtual_img_width && proj_pinhole(2)>1 && proj_pinhole(2) < virtual_img_width)
            proj_pinhole_x1 = floor(proj_pinhole(1));
            proj_pinhole_x2 = ceil(proj_pinhole(1));
            proj_pinhole_y1 = floor(proj_pinhole(2));
            proj_pinhole_y2 = ceil(proj_pinhole(2));
            
            %%  【Q11 Q12】
            %%  【Q21 Q22】
            if ((pinhole_mask{j,1}(proj_pinhole_y1, proj_pinhole_x1) > 0) && ...
                (pinhole_mask{j,1}(proj_pinhole_y1, proj_pinhole_x2) > 0) && ...
                (pinhole_mask{j,1}(proj_pinhole_y2, proj_pinhole_x1) > 0) && ...
                (pinhole_mask{j,1}(proj_pinhole_y2, proj_pinhole_x2) > 0))
                
                
                proj2 = bilinearInterp(Lut, j,proj_pinhole,  proj_pinhole_x1, proj_pinhole_x2, proj_pinhole_y1, proj_pinhole_y2);
                found = true;
                break;
            end
        end
        if 0
            
           figure,imshow(pinhole_mask{j,1}, []);hold on;plot(proj_pinhole(1),proj_pinhole(2),'.r') 
            
        end
        
    end
    
    if found
        proj_stack2(cnt2,:) = proj2;
        found_ind = [found_ind; cnt22];
        cnt2 = cnt2+1;
    end
    
    cnt22 = cnt22+1;
    
end
toc;
error_check = proj_stack1(found_ind,:) - proj_stack2;

[~, error_check_norm] = NormalizeVector(error_check);

errMatX = zeros(height, width);
errMatY = zeros(height, width);
errMat = zeros(height, width);

ind_test = valid_fisheye_area_test(found_ind);

errMatX(ind_test) = error_check(:,1);
errMatY(ind_test) = error_check(:,2);
errMat(ind_test) = error_check_norm;

figure,subplot(1,2,1);quiver(errMatX,errMatY);axis equal; subplot(1,2,2);contour(errMat, 100);title('error in pixel');colorbar; axis equal;


end


function xy_final = bilinearInterp(Lut, j,proj_pinhole, proj_pinhole_x1, proj_pinhole_x2, proj_pinhole_y1, proj_pinhole_y2)

Qx11 = Lut{j,1}(proj_pinhole_y1, proj_pinhole_x1);
Qy11 = Lut{j,2}(proj_pinhole_y1, proj_pinhole_x1);

Qx12 = Lut{j,1}(proj_pinhole_y1, proj_pinhole_x2);
Qy12 = Lut{j,2}(proj_pinhole_y1, proj_pinhole_x2);

Qx21 = Lut{j,1}(proj_pinhole_y2, proj_pinhole_x1);
Qy21 = Lut{j,2}(proj_pinhole_y2, proj_pinhole_x1);

Qx22 = Lut{j,1}(proj_pinhole_y2, proj_pinhole_x2);
Qy22 = Lut{j,2}(proj_pinhole_y2, proj_pinhole_x2);


coeff1 = (proj_pinhole_x2 - proj_pinhole(1)) * (proj_pinhole_y2 - proj_pinhole(2));
coeff2 = (proj_pinhole(1) - proj_pinhole_x1) * (proj_pinhole_y2 - proj_pinhole(2));
coeff3 = (proj_pinhole_x2 - proj_pinhole(1)) * (proj_pinhole(2) - proj_pinhole_y1);
coeff4 = (proj_pinhole(1) - proj_pinhole_x1) * (proj_pinhole(2) - proj_pinhole_y1);

tempX1 = coeff1*Qx11;
tempX2 = coeff2*Qx12;
tempX3 = coeff3*Qx21;
tempX4 = coeff4*Qx22;

tempY1 = coeff1*Qy11;
tempY2 = coeff2*Qy12;
tempY3 = coeff3*Qy21;
tempY4 = coeff4*Qy22;

xy_final = [(tempX1 + tempX2 + tempX3 + tempX4) (tempY1 + tempY2 + tempY3 + tempY4)];
end
    function [t_left, f_center, intrMat] = generateTransforms(imgWidth, fov)
        
        global CENTER_FOV;
        
        centerFOV = deg2rad(CENTER_FOV);
        sideVerticalFOV = deg2rad(fov - CENTER_FOV) /2;
        rot_angle = (centerFOV + sideVerticalFOV)/2;
        
        f_center = imgWidth / 2 / tan(centerFOV / 2);
        f_side = f_center;
        sideImgHeight = round(2 * f_side * tan(sideVerticalFOV/2));
        
        cx_side = imgWidth/2;
        cy_side = sideImgHeight/2;
        intrMat = [f_center,0,(imgWidth-1) / 2;0,f_center,(imgWidth-1) / 2;0,0,1];
        
        intr_mat = intrMat;
        
        
        t_left = cell(5, 1);
        t_left{1,1} = eye(3);
        t_left{2,1} = t_left{1,1}*rotx(-rad2deg(rot_angle));
        t_left{3,1} = t_left{2,1} * rotx(-rad2deg(pi / 2 - rot_angle)) * roty(rad2deg(pi / 2) )* rotx(rad2deg(pi / 2 - rot_angle));
        t_left{4,1} = t_left{3,1} * rotx(-rad2deg(pi / 2 - rot_angle)) * roty(rad2deg(pi / 2) )* rotx(rad2deg(pi / 2 - rot_angle));
        t_left{5,1} = t_left{4,1} * rotx(-rad2deg(pi / 2 - rot_angle)) * roty(rad2deg(pi / 2) )* rotx(rad2deg(pi / 2 - rot_angle));
        
    end
    function generateAllUndistMap(imgWidth,fov, f_center, t)
        global param CENTER_FOV;
        
        fisheye2cam_pt = zeros(480, 640);
        fisheye2cam_id = ones(480, 640);
        
        
        centerFOV = deg2rad(CENTER_FOV);
        sideVerticalFOV = deg2rad(fov - CENTER_FOV) /2;
        rot_angle = (centerFOV + sideVerticalFOV)/2;
        f_center = imgWidth / 2 / tan(centerFOV / 2);
        f_side = f_center;
        sideImgHeight = 2 * f_side * tan(sideVerticalFOV/2);
        
    end

    function [map, fisheye2cam_pt, fisheye2cam_id] = genOneUndistMap(id, rotation,imgWidth,imgHeight,f_center,fisheye2cam_pt,fisheye2cam_id)
        global param;
        map = zeros(imgHeight, imgWidth,2);
        fisheye2cam_pt = zeros(imgHeight, imgWidth,2);
        fisheye2cam_id = 255.*ones(imgHeight, imgWidth);
        if one_based
            [xMat, yMat] = meshgrid(1 : imgWidth, 1 : imgHeight);
        else
            [xMat, yMat] = meshgrid(0 : imgWidth-1, 0 : imgHeight-1);
        end
        pix = [xMat(:) yMat(:)];
        
        pt3d = (rotation * [ pix(:,1) - (imgWidth-1) / 2;    pix(:,2) - (imgHeight-1) / 2;  f_center.*ones(1,size(pix,1))])';
        
        pt2d = projectKB8(pt3d(:,1), pt3d(:,2), pt3d(:,3), param);
        
        
        map(:,:,1) = reshape(pt2d(:,1), imgHeight, imgWidth);
        map(:,:,2) = reshape(pt2d(:,2), imgHeight, imgWidth);
        
        ids = find(pt2d(:,1)>0 && pt2d(:,1)<imgWidth && pt2d(:,2)>0 && pt2d(:,2)<imgHeight);
        fisheye2cam_pt(round(pt2d(ids,2)),round(pt2d(ids,1)),1) = pix(ids,1);
        fisheye2cam_pt(round(pt2d(ids,2)),round(pt2d(ids,1)),2) = pix(ids,2);
        fisheye2cam_id(round(pt2d(ids,2)),round(pt2d(ids,1))) = id;
    end

    function [pt2d, d_uv_d_pt3d, d_uv_d_param] = projectKB8(X, Y, Z, param)
        
        d_uv_d_pt3d = [];
        d_uv_d_param = [];
        
        
        
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
        
        
        NORM_inv = 1./R;
        cosphi = x_c .* NORM_inv;
        sinphi = y_c .* NORM_inv;
        
        
        pt2d(:,1) = fx .* mx + cx;
        pt2d(:,2) = fy .* my + cy;
        
        
        
        
        if 0
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
        d_pt3d_d_uv = [];
        d_pt3d_d_param = [];
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
        
        thetad = min(max(thetad,-3.141592653/2),3.141592653/2);
        scaling = ones(length(u),1);
        d_func_d_theta = zeros(length(u),1);
        if (min(thetad) > 1e-18)
            [theta,d_func_d_theta1] = solveTheta(k1, k2, k3, k4, thetad, d_func_d_theta,5);
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
        
        pt3d = [mx .* scaling my .* scaling cos_theta];
        
        
        
        if 0
            
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
        
    end