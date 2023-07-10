function testPatchMatch(varargin)
close all;

global pix_noise
if nargin < 1
   name =  'new';
elseif nargin == 1
      name = varargin{1};
else
end

inputDir = strcat('G:\matlab\data\direct\sim\',name);

ref_patch = (imread(fullfile(inputDir, 'ref_patch.png')));
cur_patch = (imread(fullfile(inputDir, 'cur_patch.png')));
ref_pt = load(fullfile(inputDir, 'host.txt'));
ref_pt = reshape(ref_pt, 2, [])' + 1;
ref_pt_x = reshape(ref_pt(:,1),10,10);
ref_pt_y = reshape(ref_pt(:,2),10,10);
ref_px = [ref_pt_x(6,6) ref_pt_y(6,6)];
ref_offset = ref_pt - ref_px;
ref_pt_orig = ref_offset + (ref_px-0.5).*2+0.5;
cur_pt = load(fullfile(inputDir, 'target.txt'));
cur_pt = reshape(cur_pt, 2, [])' + 1;
cur_pt_x = reshape(cur_pt(:,1),8,[])';
cur_pt_y = reshape(cur_pt(:,2),8,[])';
cur_px = [cur_pt_x(end-3,end-3) cur_pt_y(end-3,end-3)];
cur_offset = cur_pt - cur_px;

figure,subplot(2,2,1);imshow(ref_patch);hold on;plot(ref_pt(101:end,1), ref_pt(101:end,2), '.b');plot(ref_pt(1:100,1), ref_pt(1:100,2), '.r');plot(ref_px(1),ref_px(2), '.g');title(sprintf('ref img size: [%d %d]', size(ref_patch,2), size(ref_patch,1)));
subplot(2,2,2);imshow(cur_patch);hold on;plot(cur_pt(1:end-64,1), cur_pt(1:end-64,2), '.b');plot(cur_pt(end-63:end,1), cur_pt(end-63:end,2), '.r');plot(cur_px(1),cur_px(2), '.g');title(sprintf('cur img size: [%d %d]', size(cur_patch,2), size(cur_patch,1)));
subplot(2,2,3);imshow(imresize(ref_patch, [480 640]));hold on;plot(ref_pt_orig(:,1), ref_pt_orig(:,2), '.r');plot( (ref_px(1)-0.5).*2+0.5, (ref_px(2)-0.5).*2+0.5, '.g');title(sprintf('ref img size: [%d %d]', 640, 480));
subplot(2,2,4);imshow(imresize(cur_patch, [480 640]));


warped_gray_val = interp2(double(cur_patch), 90.1544, 223.536, 'linear', 0);

% return;

% load('G:\matlab\LK\slam_algo\pba_large_scene_2.mat');


% check_id = 39;host_id = 1;target_id = 8;load('G:\matlab\LK\slam_algo\patch_match.mat');
check_id = 36;host_id = 1;target_id = 2;load('G:\matlab\LK\slam_algo\patch_match.mat');
% check_id = 39;host_id = 1;target_id = 2;load('G:\matlab\LK\slam_algo\patch_match2.mat');
% check_id = 43;host_id = 1;target_id = 2; load('G:\matlab\LK\slam_algo\patch_match3.mat');
% load('G:\matlab\LK\slam_algo\patch_match4.mat');% 
% check_id = 54;host_id = 1;target_id = 2;
% check_id = 54;host_id = 1;target_id = 3;





Pose = [];
deltaPoseMat = [];
LR_bak = LR;
Depth_bak = Depth;
intrMat_new = intrMat;

if 1
    intrMat_new(1,1) = intrMat_new(1,1)+100;
    intrMat_new(2,2) = intrMat_new(2,2)+100;
    intrMat_new(1,3) = intrMat_new(1,3)+30;
    intrMat_new(2,3) = intrMat_new(2,3)-20;
end

[xMat, yMat] = meshgrid(1:size(LR{1,1},2), 1:size(LR{1,1},1));
pix = [xMat(:) yMat(:)];
width = size(LR{1,1},2);
height = size(LR{1,1},1);
% pixUndist_gt = Orig2Rect(pix, intrMat, intrMat_new, eye(3), zeros(5,1));
%  pixOrig_gt2_ = remapRect(pix_gt, intrMat_gt, intrMat_virtual, R, k);
%% 用老畸变图插出新畸变图
pixDistorted1 = remapRect2(pix, intrMat,zeros(5,1), intrMat,zeros(5,1), eye(3));
pixDistorted2 = remapRect2(pix, intrMat_new,zeros(5,1), intrMat,zeros(5,1), eye(3));

pixDistorted = {pixDistorted1; pixDistorted2};

IntrMats = {intrMat; intrMat_new};

IntrMatList = {};

intr_id = [1 2 2 1 2 2 2 1 2 2 1 1 1];
intr_id = [1 1 1 1 1 1 1 1 1 1 1 1 1 ];
intr_id = [intr_id intr_id intr_id intr_id intr_id intr_id intr_id];

for i = 0 : length(LR)-1
    J0(:,:,1) = interp2(xMat,yMat,double(LR_bak{i+1,1}(:,:,1)),reshape(pixDistorted{intr_id(i+1)}(:,1),height,width),reshape(pixDistorted{intr_id(i+1)}(:,2),height,width));
    J0(:,:,2) = interp2(xMat,yMat,double(LR_bak{i+1,1}(:,:,2)),reshape(pixDistorted{intr_id(i+1)}(:,1),height,width),reshape(pixDistorted{intr_id(i+1)}(:,2),height,width));
    J0(:,:,3) = interp2(xMat,yMat,double(LR_bak{i+1,1}(:,:,3)),reshape(pixDistorted{intr_id(i+1)}(:,1),height,width),reshape(pixDistorted{intr_id(i+1)}(:,2),height,width));
    J1 = interp2(xMat,yMat,double(Depth_bak{i+1,1}(:,:,1)),reshape(pixDistorted{intr_id(i+1)}(:,1),height,width),reshape(pixDistorted{intr_id(i+1)}(:,2),height,width));
    J1(isnan(J1)) = inf;
    
    if 0
        figure,imshow([uint8(J0) LR_bak{i+1,1}])
        figure,imshow([J1 Depth_bak{i+1,1}], [])
    end
    IntrMatList{i+1,1} = IntrMats{intr_id(i+1)};
    IntrMatList{i+1,2} = intr_id(i+1);
    LR{i+1,1} = uint8(J0);
    Depth{i+1,1} = J1;
    pose = inv(k2c_comp{i+1});
    Pose = [Pose; rodrigues(pose(1:3,1:3))' pose(1:3,4)'];
    deltaPoseMat = [deltaPoseMat; [reshape(pose(1:3,1:3),1,9) pose(1:3,4)']];
end




host_img = imgaussfilt(rgb2gray(LR{host_id, 1}), 0.6);
target_img = imgaussfilt(rgb2gray(LR{target_id, 1}), 0.6);
host_depth = Depth{host_id, 1};
target_depth = Depth{target_id, 1};
Twc_h = [reshape(deltaPoseMat(host_id, 1 : 9), 3, 3) deltaPoseMat(host_id, 10:12)';0 0 0 1];
Twc_t = [reshape(deltaPoseMat(target_id, 1 : 9), 3, 3) deltaPoseMat(target_id, 10:12)';0 0 0 1];
T_th = inv(Twc_t) * Twc_h;

host_pt = detectFASTFeatures(host_img,'MinQuality',0.3,'MinContrast',0.3);
host_px = double(host_pt.Location);
ind = sub2ind(size(host_img), host_px(:,2), host_px(:,1));
depths = host_depth(ind);

host_px_norm = inv(intrMat) * [pextend(host_px')];
host_xyz = repmat(depths,1,3) .* host_px_norm';
target_xyz = (T_th(1:3, 1:3) * host_xyz' + repmat(T_th(1:3,4),1,size(host_xyz,1)))';
target_proj = pflat(intrMat * target_xyz');
target_proj = target_proj(1:2,:)';

figure,subplot(1,2,1),imshow(host_img);hold on;plot(host_px(:,1),host_px(:,2), '.r'); plot(host_px(check_id,1),host_px(check_id,2), '*b');
subplot(1,2,2),imshow(target_img);hold on;plot(target_proj(:,1), target_proj(:,2), '.g'); plot(target_proj(check_id,1),target_proj(check_id,2), '*b');


gx = zeros(size(host_img));
gy = gx;

gx(:,2 : width-1) = 0.5 * (double(host_img(:,3:width)) - double(host_img(:, 1:width-2)));
gy(2 : height-1, :) = 0.5 * (double(host_img(3:height,:)) - double(host_img(1:height-2,:)));


max_levels = 5;
level_ref = 2;0; 1;2;0;1;
halfpatch_size = 4;2;10;4;
cnt = 1;
loss = {};
errVec = [];
iterNum = size(host_px,1);
noise_level = -3:3;
for j = noise_level
    pix_noise = j;
    for  check_id = 1 : iterNum
        close all;
        depth = depths(check_id);
        depth_ref = depth;
        xyz_ref = host_xyz(check_id, :);
        depth_ref_norm = norm(xyz_ref);
        % warped_gx = interp2(double(gx1), pt1s(:, 1), pt1s(:, 2), 'linear', 0);
        % warped_gy = interp2(double(gy1), pt1s(:, 1), pt1s(:, 2), 'linear', 0);
        px_ref = host_px(check_id,:);
        px_cur = pflat(intrMat * (T_th(1:3,1:3)*xyz_ref' + T_th(1:3,4)));
        if (px_cur(1) < 60 || px_cur(2) < 60 || px_cur(1) > 640 - 60 || px_cur(2) > 480-60)
           continue; 
        end
        [A_th, ~, epiDir] = GetWarpMatrixAffine(intrMat, px_ref, level_ref, halfpatch_size, T_th, xyz_ref);
        
        
        search_level = GetBestSearchLevel(A_th, max_levels - 1);
        % search_level = 0;
        
        pc_cur_check = target_proj(check_id,:);
        
        check_err = pc_cur_check - px_cur(1:2)';
        
        intrMat_scaled = intrMat./2;
        intrMat_scaled(3,3) = 1;
        intrMat_scaled(1,3) = (intrMat(1,3)-0.5) /2 +0.5;
        intrMat_scaled(2,3) = (intrMat(2,3)-0.5) /2 +0.5;
        pc_cur_scaled = pflat(intrMat_scaled * target_xyz(check_id,:)');
        pc_cur_scaled = pc_cur_scaled(1:2)';
        pc_cur_check_error = 2*(pc_cur_scaled-0.5)+0.5 - pc_cur_check;
        
        
        % is_host_square = false;
        [host_offset1, target_offset1,  px_host_scaled1, px_target_scaled1, host_image_scaled1, target_image_scaled1] = WarpAffine(A_th, host_img, px_ref, level_ref, target_img, pc_cur_check, search_level, halfpatch_size+1);
        
        % is_host_square = true;
        % WarpAffine(inv(A_th), target_img, pc_cur_check, search_level, host_img, px_ref, level_ref, halfpatch_size+1);
        [A_th_scaled, ~] = GetWarpMatrixAffineInv(intrMat, px_ref, level_ref, halfpatch_size, T_th, xyz_ref);
        search_level_scaled = GetBestSearchLevel(inv(A_th_scaled), max_levels - 1);
        
        extra_search_level = 0;1;
        
        [host_offset2, target_offset2, px_host_scaled2, px_target_scaled2, host_image_scaled2, target_image_scaled2] = WarpAffineInv(A_th, host_img, px_ref, level_ref, target_img, pc_cur_check, search_level+extra_search_level, halfpatch_size+1);
        % WarpAffineInv2(A_th_scaled, host_img, px_ref, level_ref, target_img, pc_cur_check, search_level_scaled, halfpatch_size+1);
        
        host_offset1_check = (inv(A_th) * target_offset1')'.*2^(search_level);
        target_offset2_check = ((A_th) * host_offset2')'.*2^(-(search_level+extra_search_level));
        % figure,imshow(zeros(480, 640));hold on;plot(host_offset1_check(:,1) + 100, host_offset1_check(:,2) + 100,'.r');
        
        fprintf('============================================================\n');
        [px_cur_error1, update_vec1] = align1d_zncc(eye(2), host_image_scaled1, target_image_scaled1, px_host_scaled1, px_target_scaled1, host_offset1, target_offset1, epiDir, epiDir, 2*(halfpatch_size+1));
        % fprintf('============================================================\n');
        % align1d_zncc(A_th, host_image_scaled2, target_image_scaled2, px_host_scaled2, px_target_scaled2, host_offset2, target_offset2, epiDir, epiDir, 2*(halfpatch_size+1));
        fprintf('============================================================\n');
        [px_cur_error2, update_vec2] = align1d_zncc2(A_th, search_level, host_image_scaled2, target_image_scaled2, px_host_scaled2, px_target_scaled2, host_offset2, target_offset2, epiDir, epiDir, 2*(halfpatch_size+1));
        
        errVec = [errVec; [px_cur_error1 px_cur_error2]];
        loss{cnt,1} = update_vec1;
        loss{cnt,2} = update_vec2;
        cnt = cnt+1;
    end
end

end
function [px_cur_error, update_vec] = align1d_zncc(A_th, host_image, target_image,  px_host_scaled, px_target_scaled, host_offset, target_offset, epiDir, epiDir_noise, patch_size)
global pix_noise
epiDir0 = epiDir;

use_b_only = false;


target_check = A_th * px_host_scaled';

figure,subplot(1,2,1);imshow(host_image);hold on;plot(host_offset(:,1) + px_host_scaled(:,1), host_offset(:,2) + px_host_scaled(:,2), '.r');
plot(px_host_scaled(:,1), px_host_scaled(:,2),'*g');
subplot(1,2,2);imshow(target_image);hold on;plot(target_offset(:,1) + px_target_scaled(:,1), target_offset(:,2) + px_target_scaled(:,2), '.r');
plot(px_target_scaled(:,1), px_target_scaled(:,2),'*g');

width = size(host_image,2);
height = size(host_image,1);

gx = zeros(size(host_image));
gy = gx;

gx(:,2 : width-1) = 0.5 * (double(host_image(:,3:width)) - double(host_image(:, 1:width-2)));
gy(2 : height-1, :) = 0.5 * (double(host_image(3:height,:)) - double(host_image(1:height-2,:)));



warped_host_val = interp2(double(host_image), host_offset(:, 1) + px_host_scaled(:,1), host_offset(:, 2) + px_host_scaled(:,2), 'linear', 0);
warped_host_val0 = warped_host_val;
warped_gx_val = interp2(gx, host_offset(:, 1) + px_host_scaled(:,1), host_offset(:, 2) + px_host_scaled(:,2), 'linear', 0);
warped_gy_val = interp2(gy, host_offset(:, 1) + px_host_scaled(:,1), host_offset(:, 2) + px_host_scaled(:,2), 'linear', 0);

gxxx = warped_gx_val;
gyyy = warped_gy_val;
gxxx0 = gxxx;
gyyy0 = gyyy;
if 0
    % correct gx gy
    warped_host_val_mat = reshape(warped_host_val, patch_size, patch_size);
    gxx = zeros(size(warped_host_val_mat));
    gyy = gxx;
    gxx(:,2 : size(warped_host_val_mat,2)-1) = 0.5 * (double(warped_host_val_mat(:,3:size(warped_host_val_mat,2))) - double(warped_host_val_mat(:, 1:size(warped_host_val_mat,2)-2)));
    gyy(2 : size(warped_host_val_mat,1)-1, :) = 0.5 * (double(warped_host_val_mat(3:size(warped_host_val_mat,1),:)) - double(warped_host_val_mat(1:size(warped_host_val_mat,1)-2,:)));
    
    
    
    warped_host_val_mat = warped_host_val_mat(2:end-1,2:end-1);
    warped_host_val = warped_host_val_mat(:);
        
    host_offset_x = shrinkMat(host_offset(:,1), patch_size);
    host_offset_y = shrinkMat(host_offset(:,2),patch_size);
    host_offset = [host_offset_x host_offset_y];
    target_offset_x = shrinkMat(target_offset(:,1), patch_size);
    target_offset_y = shrinkMat(target_offset(:,2),patch_size);
    target_offset = [target_offset_x target_offset_y];
    gxxx = shrinkMat(gxx(:), patch_size);
    gyyy = shrinkMat(gyy(:),patch_size);
else
    warped_host_val = shrinkMat(warped_host_val(:,1), patch_size);
     host_offset_x = shrinkMat(host_offset(:,1), patch_size);
    host_offset_y = shrinkMat(host_offset(:,2),patch_size);
    host_offset = [host_offset_x host_offset_y];
    target_offset_x = shrinkMat(target_offset(:,1), patch_size);
    target_offset_y = shrinkMat(target_offset(:,2),patch_size);
    target_offset = [target_offset_x target_offset_y];
    
    gxxx = shrinkMat(warped_gx_val, patch_size);
    gyyy = shrinkMat(warped_gy_val,patch_size);
end

warped_host_val = warped_host_val - mean(warped_host_val);
if ~use_b_only
    sigma0 = norm(warped_host_val);
    warped_host_val = warped_host_val./sigma0;
end


Mat_I = eye(size(host_offset,1));
Mat_1 = ones(size(host_offset,1), size(host_offset,1));

if use_b_only
    J_zncc = (Mat_I - Mat_1 / size(host_offset,1));
else
    J_zncc = (Mat_I - (warped_host_val * warped_host_val')) / sigma0 * (Mat_I - Mat_1 / (size(host_offset,1)));
end

epiDir = A_th * epiDir0';
epiDir = epiDir'./norm(epiDir);
JznccMat = J_zncc * [gxxx gyyy] * inv(A_th) *  epiDir';
% JznccMat = J_zncc * [gxxx gyyy] * epiDir';

Hess = JznccMat' * JznccMat;
Hess_inv = 1 / Hess;


px_target_scaled_error = px_target_scaled + pix_noise.*epiDir_noise;

u = px_target_scaled_error(1);
v = px_target_scaled_error(2);

chi2 = 0;
update = 0;



iter_num = 15;
min_update_squared = 0.001^2;


update_vec = [];
res = [];
color = {'.g','.b','.m','.c','.y','.g','.b','.m','.c','.y','.g','.b','.m','.c','.y','.g','.b','.m','.c','.y','.g','.b','.m','.c','.y','.g','.b','.m','.c','.y'};
for iter = 1 : iter_num
    
    u_r = floor(u);
    v_r = floor(v);
    
    new_chi2 = 0;
    warped_target_val = interp2(double(target_image), target_offset(:, 1) + u, target_offset(:, 2) + v, 'linear', 0);
    warped_target_val = warped_target_val - mean(warped_target_val);
    if ~use_b_only
        sigma1 = norm(warped_target_val);
        warped_target_val = warped_target_val./sigma1;
    end
    rs = warped_target_val - warped_host_val;
    new_chi2 = sum(rs.^2);
    res = [res new_chi2];
    if iter > 1 && new_chi2 > chi2
       u = u - update * epiDir(1);
       v = v - update * epiDir(2);
       fprintf(sprintf('error increased, new_chi2: %f, update: %f\n', new_chi2,update));
       break;
    end
    chi2 = new_chi2;
    update = -Hess_inv * rs' * JznccMat;
    update_vec = [update_vec  update];
    if 1
        uv_new = A_th * [u; v] + update*epiDir';
        uv_new_coord = inv(A_th) * uv_new;
        u = uv_new_coord(1);
        v = uv_new_coord(2);
    else
        u = u + update * epiDir(1);
        v = v + update * epiDir(2);
    end
    subplot(1,2,2);plot(u, v, color{iter});
    if (update * update < min_update_squared) 
      converged = true;
      break;
    end
end
update_vec
% res
px_cur_error = norm(px_target_scaled - [u v])

end
function [px_cur_error, update_vec] = align1d_zncc2(A_th_, search_level, host_image, target_image,  px_host_scaled, px_target_scaled, host_offset, target_offset, epiDir, epiDir_noise, patch_size)
global pix_noise
use_b_only = false;

wrong_code = falstrue;

target_offset_check = ((A_th_) * host_offset')'.*2^(-(search_level));

target_patch = px_target_scaled + target_offset;
host_px = px_host_scaled + [-3 -3];
target_center_check = ((A_th_) * (host_px - px_host_scaled)')'.*2^(-(search_level)) + px_target_scaled;
A_th = A_th_.*2^(-search_level);
const = -(A_th*px_host_scaled')' + px_target_scaled;
A = [A_th const';0 0 1];
target_center_check2 = A_th * host_px' + const';
target_center_check3 = A*[host_px 1]';


epiDir0 = epiDir;

if wrong_code
    epiDir = A_th * epiDir0';
    epiDir = epiDir'./norm(epiDir);
end
target_offset0 = target_offset;

figure,subplot(1,2,1);imshow(host_image);hold on;plot(host_offset(:,1) + px_host_scaled(:,1), host_offset(:,2) + px_host_scaled(:,2), '.r');
plot(px_host_scaled(:,1), px_host_scaled(:,2),'*g');
subplot(1,2,2);imshow(target_image);hold on;plot(target_offset(:,1) + px_target_scaled(:,1), target_offset(:,2) + px_target_scaled(:,2), '.r');
plot(px_target_scaled(:,1), px_target_scaled(:,2),'*g');

width = size(host_image,2);
height = size(host_image,1);

gx = zeros(size(host_image));
gy = gx;

gx(:,2 : width-1) = 0.5 * (double(host_image(:,3:width)) - double(host_image(:, 1:width-2)));
gy(2 : height-1, :) = 0.5 * (double(host_image(3:height,:)) - double(host_image(1:height-2,:)));



warped_host_val = interp2(double(host_image), host_offset(:, 1) + px_host_scaled(:,1), host_offset(:, 2) + px_host_scaled(:,2), 'linear', 0);
warped_host_val0 = warped_host_val;
warped_gx_val = interp2(gx, host_offset(:, 1) + px_host_scaled(:,1), host_offset(:, 2) + px_host_scaled(:,2), 'linear', 0);
warped_gy_val = interp2(gy, host_offset(:, 1) + px_host_scaled(:,1), host_offset(:, 2) + px_host_scaled(:,2), 'linear', 0);

gxxx = warped_gx_val;
gyyy = warped_gy_val;
gxxx0 = gxxx;
gyyy0 = gyyy;
if 0
    % correct gx gy
    warped_host_val_mat = reshape(warped_host_val, patch_size, patch_size);
    gxx = zeros(size(warped_host_val_mat));
    gyy = gxx;
    gxx(:,2 : size(warped_host_val_mat,2)-1) = 0.5 * (double(warped_host_val_mat(:,3:size(warped_host_val_mat,2))) - double(warped_host_val_mat(:, 1:size(warped_host_val_mat,2)-2)));
    gyy(2 : size(warped_host_val_mat,1)-1, :) = 0.5 * (double(warped_host_val_mat(3:size(warped_host_val_mat,1),:)) - double(warped_host_val_mat(1:size(warped_host_val_mat,1)-2,:)));
    
    
    
    warped_host_val_mat = warped_host_val_mat(2:end-1,2:end-1);
    warped_host_val = warped_host_val_mat(:);
        
    host_offset_x = shrinkMat(host_offset(:,1), patch_size);
    host_offset_y = shrinkMat(host_offset(:,2),patch_size);
    host_offset = [host_offset_x host_offset_y];
    target_offset_x = shrinkMat(target_offset(:,1), patch_size);
    target_offset_y = shrinkMat(target_offset(:,2),patch_size);
    target_offset = [target_offset_x target_offset_y];
    gxxx = shrinkMat(gxx(:), patch_size);
    gyyy = shrinkMat(gyy(:),patch_size);
else
    warped_host_val = shrinkMat(warped_host_val(:,1), patch_size);
     host_offset_x = shrinkMat(host_offset(:,1), patch_size);
    host_offset_y = shrinkMat(host_offset(:,2),patch_size);
    host_offset = [host_offset_x host_offset_y];
    target_offset_x = shrinkMat(target_offset(:,1), patch_size);
    target_offset_y = shrinkMat(target_offset(:,2),patch_size);
    target_offset = [target_offset_x target_offset_y];
    
    gxxx = shrinkMat(warped_gx_val, patch_size);
    gyyy = shrinkMat(warped_gy_val,patch_size);
end

warped_host_val = warped_host_val - mean(warped_host_val);
if ~use_b_only
    sigma0 = norm(warped_host_val);
    warped_host_val = warped_host_val./sigma0;
end


Mat_I = eye(size(host_offset,1));
Mat_1 = ones(size(host_offset,1), size(host_offset,1));

if use_b_only
    J_zncc = (Mat_I - Mat_1 / size(host_offset,1));
else
    J_zncc = (Mat_I - (warped_host_val * warped_host_val')) / sigma0 * (Mat_I - Mat_1 / (size(host_offset,1)));
end
% JznccMat = J_zncc * [gxxx gyyy] * epiDir';
if ~wrong_code
    JznccMat = J_zncc * [gxxx gyyy] * inv(A_th) * epiDir0';
else
    JznccMat = J_zncc * [gxxx gyyy] * inv(A_th) * epiDir';
end
Hess = JznccMat' * JznccMat;
Hess_inv = 1 / Hess;


px_target_scaled_error = px_target_scaled + pix_noise.*epiDir_noise;

u = px_target_scaled_error(1);
v = px_target_scaled_error(2);

chi2 = 0;
update = 0;



iter_num = 15;
min_update_squared = 0.001^2;


update_vec = [];
res = [];
color = {'.g','.b','.m','.c','.y','.g','.b','.m','.c','.y','.g','.b','.m','.c','.y','.g','.b','.m','.c','.y','.g','.b','.m','.c','.y','.g','.b','.m','.c','.y'};
for iter = 1 : iter_num
    
    u_r = floor(u);
    v_r = floor(v);
    
    new_chi2 = 0;
    warped_target_val_pad = interp2(double(target_image), target_offset0(:, 1) + u, target_offset0(:, 2) + v, 'linear', 0);
    warped_target_val = interp2(double(target_image), target_offset(:, 1) + u, target_offset(:, 2) + v, 'linear', 0);
    
    warped_target_val_pad_mat = reshape(warped_target_val_pad, patch_size, patch_size);
    gx_target = zeros(patch_size, patch_size);
    gy_target = gx_target;
    gx_target(:,2 : patch_size-1) = 0.5 * (double(warped_target_val_pad_mat(:,3:patch_size)) - double(warped_target_val_pad_mat(:, 1:patch_size-2)));
    gy_target(2 : patch_size-1, :) = 0.5 * (double(warped_target_val_pad_mat(3:patch_size,:)) - double(warped_target_val_pad_mat(1:patch_size-2,:)));
    
    gx_vec =  shrinkMat(gx_target(:), patch_size);
    gy_vec =  shrinkMat(gy_target(:), patch_size);
    
    warped_target_val = warped_target_val - mean(warped_target_val);
    if ~use_b_only
        sigma1 = norm(warped_target_val);
        warped_target_val = warped_target_val./sigma1;
    end
    rs = warped_target_val - warped_host_val;
    new_chi2 = sum(rs.^2);
    res = [res new_chi2];
    if iter > 1 && new_chi2 > chi2
       u = u - update_value(1) * epiDir0(1);
       v = v - update_value(1) * epiDir0(2);
       fprintf(sprintf('error increased, new_chi2: %f, update: %f\n', new_chi2,update));
       break;
    end
    chi2 = new_chi2;
    update = -Hess_inv * rs' * JznccMat;
%     update_vec = [update_vec  update];
    
    if wrong_code
        uv_new = A_th * [u; v] + update*epiDir';
        uv_new_coord = inv(A_th) * uv_new;
        update_orig = uv_new_coord' - [u v];
        update_value = update_orig./epiDir0;
        update_vec = [update_vec  update_value(1)];
%         u = uv_new_coord(1);
%         v = uv_new_coord(2);
        u = u + update_value(1) * epiDir0(1);
        v = v + update_value(1) * epiDir0(2);
    else
        update_value = update;
        update_vec = [update_vec  update_value(1)];
        u = u + update * epiDir0(1);
        v = v + update * epiDir0(2);
    end
    subplot(1,2,2);plot(u, v, color{iter});
    if (update * update < min_update_squared) 
      converged = true;
      break;
    end
end
update_vec
% res
px_cur_error = norm(px_target_scaled - [u v])

end
function host_offset_ = shrinkMat(host_offset, patch_size)
host_offset_x = reshape(host_offset,patch_size,patch_size);
host_offset_x = host_offset_x(2:end-1, 2:end-1);
host_offset_ = host_offset_x(:);


end
function search_level = GetBestSearchLevel(A_cur_ref, max_level)
search_level = 0;
D0 = det(A_cur_ref);
D = det(A_cur_ref);
while(D > 1.5 && search_level < max_level)
    search_level = search_level + 1;
    D = D * 0.25;
end

end
function [host_offset, target_offset, px_host_scaled, px_target_scaled, host_img0, target_img0] = WarpAffine(A_cur_ref, host_img, px_ref, level_ref, target_img, px_cur, search_level, halfpatch_size)
patch_size = halfpatch_size * 2;
patch_size_use = (halfpatch_size-1) * 2;
A_ref_cur = inv(A_cur_ref);
host_img0 = imresize(host_img, size(host_img)* 2^(-level_ref));
target_img0 = imresize(target_img, size(target_img)* 2^(-search_level));
px_ref_ = (px_ref-0.5) * 2^(-level_ref) + 0.5;
[xMat, yMat] = meshgrid(0:patch_size-1, 0:patch_size-1);
pix0 = [xMat(:) yMat(:)];
pix0 = pix0-halfpatch_size;
% pix = (pix0 -0.5) * 2^(search_level) + 0.5;
% svo 2.0 not
if 1
    pix = (pix0 ) * 2^(search_level);
else
    pix = pix0;
end
% pix = (pix0 ) * 2^(level_ref);
% pix = (pix0 -0.5) * 2^(0) + 0.5;
pix_ = (A_ref_cur * pix')' + px_ref_;       % 缩放尺度
% patch_bound = (A_ref_cur * pix')' + (px_ref_-0.5) * 2^(search_level) + 0.5; % 原尺寸
patch_bound = (A_ref_cur * pix')' + (px_ref_-0.5) * 2^(level_ref) + 0.5; % 原尺寸
warped_patchs = interp2(double(host_img0), pix_(:, 1), pix_(:, 2), 'linear', 0);
check_id = 3;
x_final = bilinearInterp(double(host_img0), pix_(check_id,:), floor(pix_(check_id,1)), ceil(pix_(check_id,1)), floor(pix_(check_id,2)), ceil(pix_(check_id,2)));
error = x_final - warped_patchs(check_id,:);

%   warped_patch = [];
%   pixs = [];
%   for  y = 2 : patch_size_use + 1
%     ref_patch_border_ptr = patch_with_border_ + y * (patch_size_ + 2) + 1;
%     for x = 1 : patch_size_use
%       ref_patch_ptr[x] = ref_patch_border_ptr[x];
%     end
%   end

px_cur_scaled = (px_cur -0.5) * 2^(-search_level) + 0.5;
px_curs = px_cur + pix0;
px_curs_scaled = px_cur_scaled + pix0;
figure,subplot(2,2,1);imshow(host_img0);hold on;plot(pix_(:,1), pix_(:,2),'.r');plot(px_ref_(:,1), px_ref_(:,2),'.g');title(sprintf('small image, covers more area, patch src, level: %d',level_ref));
subplot(2,2,2);imshow(target_img0);hold on;plot(px_curs_scaled(:,1), px_curs_scaled(:,2),'.r');plot(px_cur_scaled(:,1), px_cur_scaled(:,2),'.g');title(sprintf('target patch at search level: %d', search_level));
subplot(2,2,3);imshow(host_img);hold on;plot(patch_bound(:,1), patch_bound(:,2),'.r');plot(px_ref(:,1), px_ref(:,2),'.g');title('orig host image, covers less area, level: 0');
subplot(2,2,4);imshow(target_img);hold on;plot(px_curs(:,1), px_curs(:,2),'.r');plot(px_cur(:,1), px_cur(:,2),'.g');title('orig target image, covers less area, level: 0');



 host_offset = pix_ - px_ref_;
 target_offset = pix0;
px_host_scaled = px_ref_;
px_target_scaled = px_cur_scaled;

end
function [host_offset, target_offset, px_host_scaled, px_target_scaled, host_img0, target_img0] = WarpAffineInv(A_cur_ref, host_img, px_ref, level_ref, target_img, px_cur, search_level, halfpatch_size)
patch_size = halfpatch_size * 2;
patch_size_use = (halfpatch_size-1) * 2;
% A_cur_ref = inv(A_ref_cur);
host_img0 = imresize(host_img, size(host_img)* 2^(-level_ref));
target_img0 = imresize(target_img, size(target_img)* 2^(-search_level));
px_cur_ = (px_cur-0.5) * 2^(-search_level) + 0.5;
[xMat, yMat] = meshgrid(0:patch_size-1, 0:patch_size-1);
pix0 = [xMat(:) yMat(:)];
pix0 = pix0-halfpatch_size;
% pix = (pix0 -0.5) * 2^(search_level) + 0.5;
% svo 2.0 not
if 1
    pix = (pix0 ) * 2^(-search_level);
else
    pix = pix0;
end
% pix = (pix0 ) * 2^(level_ref);
% pix = (pix0 -0.5) * 2^(0) + 0.5;
pix_ = (A_cur_ref * pix')' + px_cur_;       % 缩放尺度
patch_bound = (A_cur_ref * pix')' + (px_cur_-0.5) * 2^(search_level) + 0.5; % 原尺寸
warped_patchs = interp2(double(target_img0), pix_(:, 1), pix_(:, 2), 'linear', 0);
check_id = 3;
x_final = bilinearInterp(double(target_img0), pix_(check_id,:), floor(pix_(check_id,1)), ceil(pix_(check_id,1)), floor(pix_(check_id,2)), ceil(pix_(check_id,2)));
error = x_final - warped_patchs(check_id,:);

%   warped_patch = [];
%   pixs = [];
%   for  y = 2 : patch_size_use + 1
%     ref_patch_border_ptr = patch_with_border_ + y * (patch_size_ + 2) + 1;
%     for x = 1 : patch_size_use
%       ref_patch_ptr[x] = ref_patch_border_ptr[x];
%     end
%   end

px_host_scaled = (px_ref -0.5) * 2^(-level_ref) + 0.5;
px_hosts = px_ref + pix0;
px_hosts_scaled = px_host_scaled + pix0;
figure,subplot(2,2,1);imshow(host_img0);hold on;plot(px_hosts_scaled(:,1), px_hosts_scaled(:,2),'.r');plot(px_host_scaled(:,1), px_host_scaled(:,2),'.g');title(sprintf('small image, covers more area, patch src, level: %d',level_ref));
subplot(2,2,2);imshow(target_img0);hold on;plot(pix_(:,1), pix_(:,2),'.r');plot(px_cur_(:,1), px_cur_(:,2),'.g');title(sprintf('target patch at search level: %d', search_level));
subplot(2,2,3);imshow(host_img);hold on;plot(px_hosts(:,1), px_hosts(:,2),'.r');plot(px_ref(:,1), px_ref(:,2),'.g');title('orig host image, covers less area, level: 0');
subplot(2,2,4);imshow(target_img);hold on;plot(patch_bound(:,1), patch_bound(:,2),'.r');plot(px_cur(:,1), px_cur(:,2),'.g');title('orig target image, covers less area, level: 0');

host_offset = pix0;
target_offset = pix_ - px_cur_;
px_target_scaled = px_cur_;
end
function WarpAffineInv2(A_cur_ref, host_img, px_ref, level_ref, target_img, px_cur, search_level, halfpatch_size)
patch_size = halfpatch_size * 2;
patch_size_use = (halfpatch_size-1) * 2;
% A_cur_ref = inv(A_ref_cur);
host_img0 = imresize(host_img, size(host_img)* 2^(-level_ref));
target_img0 = imresize(target_img, size(target_img)* 2^(-search_level));
px_cur_ = (px_cur-0.5) * 2^(-search_level) + 0.5;
[xMat, yMat] = meshgrid(0:patch_size-1, 0:patch_size-1);
pix0 = [xMat(:) yMat(:)];
pix0 = pix0-halfpatch_size;
% pix = (pix0 -0.5) * 2^(search_level) + 0.5;
% svo 2.0 not
if 1
    pix = (pix0 ) * 2^(search_level);
else
    pix = pix0;
end
% pix = (pix0 ) * 2^(level_ref);
% pix = (pix0 -0.5) * 2^(0) + 0.5;
pix_ = (A_cur_ref * pix')' + px_cur_;       % 缩放尺度
patch_bound = (A_cur_ref * pix')' + (px_cur_-0.5) * 2^(search_level) + 0.5; % 原尺寸
warped_patchs = interp2(double(target_img0), pix_(:, 1), pix_(:, 2), 'linear', 0);
check_id = 3;
x_final = bilinearInterp(double(target_img0), pix_(check_id,:), floor(pix_(check_id,1)), ceil(pix_(check_id,1)), floor(pix_(check_id,2)), ceil(pix_(check_id,2)));
error = x_final - warped_patchs(check_id,:);

%   warped_patch = [];
%   pixs = [];
%   for  y = 2 : patch_size_use + 1
%     ref_patch_border_ptr = patch_with_border_ + y * (patch_size_ + 2) + 1;
%     for x = 1 : patch_size_use
%       ref_patch_ptr[x] = ref_patch_border_ptr[x];
%     end
%   end

px_host_scaled = (px_ref -0.5) * 2^(-level_ref) + 0.5;
px_hosts = px_ref + pix0;
px_hosts_scaled = px_host_scaled + pix0;
figure,subplot(2,2,1);imshow(host_img0);hold on;plot(px_hosts_scaled(:,1), px_hosts_scaled(:,2),'.r');plot(px_host_scaled(:,1), px_host_scaled(:,2),'.g');title(sprintf('small image, covers more area, patch src, level: %d',level_ref));
subplot(2,2,2);imshow(target_img0);hold on;plot(pix_(:,1), pix_(:,2),'.r');plot(px_cur_(:,1), px_cur_(:,2),'.g');title(sprintf('target patch at search level: %d', search_level));
subplot(2,2,3);imshow(host_img);hold on;plot(px_hosts(:,1), px_hosts(:,2),'.r');plot(px_ref(:,1), px_ref(:,2),'.g');title('orig host image, covers less area, level: 0');
subplot(2,2,4);imshow(target_img);hold on;plot(patch_bound(:,1), patch_bound(:,2),'.r');plot(px_cur(:,1), px_cur(:,2),'.g');title('orig target image, covers less area, level: 0');



end
function x_final = bilinearInterp(Lut, proj_pinhole, proj_pinhole_x1, proj_pinhole_x2, proj_pinhole_y1, proj_pinhole_y2)

Qx11 = Lut(proj_pinhole_y1, proj_pinhole_x1);
% Qy11 = Lut{j,2}(proj_pinhole_y1, proj_pinhole_x1);

Qx12 = Lut(proj_pinhole_y1, proj_pinhole_x2);
% Qy12 = Lut{j,2}(proj_pinhole_y1, proj_pinhole_x2);

Qx21 = Lut(proj_pinhole_y2, proj_pinhole_x1);
% Qy21 = Lut{j,2}(proj_pinhole_y2, proj_pinhole_x1);

Qx22 = Lut(proj_pinhole_y2, proj_pinhole_x2);
% Qy22 = Lut{j,2}(proj_pinhole_y2, proj_pinhole_x2);


coeff1 = (proj_pinhole_x2 - proj_pinhole(1)) * (proj_pinhole_y2 - proj_pinhole(2));
coeff2 = (proj_pinhole(1) - proj_pinhole_x1) * (proj_pinhole_y2 - proj_pinhole(2));
coeff3 = (proj_pinhole_x2 - proj_pinhole(1)) * (proj_pinhole(2) - proj_pinhole_y1);
coeff4 = (proj_pinhole(1) - proj_pinhole_x1) * (proj_pinhole(2) - proj_pinhole_y1);

tempX1 = coeff1*Qx11;
tempX2 = coeff2*Qx12;
tempX3 = coeff3*Qx21;
tempX4 = coeff4*Qx22;

% tempY1 = coeff1*Qy11;
% tempY2 = coeff2*Qy12;
% tempY3 = coeff3*Qy21;
% tempY4 = coeff4*Qy22;

% xy_final = [(tempX1 + tempX2 + tempX3 + tempX4) (tempY1 + tempY2 + tempY3 + tempY4)];
x_final = [(tempX1 + tempX2 + tempX3 + tempX4)];
end
function [A_th, px_cur, epidir] = GetWarpMatrixAffine(intrMat, px_ref, level_ref, halfpatch_size, T_th, xyz_ref)

% xyz_ref(3) = xyz_ref(3) + 100;

is_fisheye = 0;1;

xyz_du_ref = inv(intrMat) * [px_ref + 2^level_ref * [halfpatch_size+1 0] 1]';
xyz_dv_ref = inv(intrMat) * [px_ref + 2^level_ref * [0 halfpatch_size+1] 1]';


if ~is_fisheye
    % 假设patch内的z分量一样
    xyz_du_ref1 = xyz_du_ref * xyz_ref(3) / xyz_du_ref(3);
    xyz_dv_ref1 = xyz_dv_ref * xyz_ref(3) / xyz_dv_ref(3);
    xyz_du_ref = xyz_du_ref1;
    xyz_dv_ref = xyz_dv_ref1;
else
    % 假设patch内的xyz的模长一样
    xyz_du_ref0 = xyz_du_ref./norm(xyz_du_ref);
    xyz_dv_ref0 = xyz_dv_ref./norm(xyz_dv_ref);
    xyz_du_ref2 = xyz_du_ref0 * depth_ref_norm;
    xyz_dv_ref2 = xyz_dv_ref0 * depth_ref_norm;
    xyz_du_ref = xyz_du_ref2;
    xyz_dv_ref = xyz_dv_ref2;
end



xyz_ref2 = 1.01.*xyz_ref;
px_cur2 = pflat(intrMat * (T_th(1:3,1:3)*xyz_ref2' + T_th(1:3,4)));

px_cur = pflat(intrMat * (T_th(1:3,1:3)*xyz_ref' + T_th(1:3,4)));

epidir = px_cur2(1:2) - px_cur(1:2);
epidir = epidir'./norm(epidir);

px_du = pflat(intrMat * (T_th(1:3,1:3)*xyz_du_ref + T_th(1:3,4)));
px_dv = pflat(intrMat * (T_th(1:3,1:3)*xyz_dv_ref + T_th(1:3,4)));
A_th = zeros(2,2);
A_th(:,1) = (px_du(1:2) - px_cur(1:2))/(halfpatch_size+1);
A_th(:,2) = (px_dv(1:2) - px_cur(1:2))/(halfpatch_size+1);
end
function [A_th, px_cur] = GetWarpMatrixAffineInv(intrMat, px_ref, level_ref, halfpatch_size, T_th, xyz_ref)

is_fisheye = 0;1;

if 0
    xyz_du_ref = inv(intrMat) * [px_ref + 2^level_ref * [halfpatch_size+1 0] 1]';
    xyz_dv_ref = inv(intrMat) * [px_ref + 2^level_ref * [0 halfpatch_size+1] 1]';
else
    xyz_du_ref = inv(intrMat) * [px_ref + 2^(-level_ref) * [halfpatch_size+1 0] 1]';
    xyz_dv_ref = inv(intrMat) * [px_ref + 2^(-level_ref) * [0 halfpatch_size+1] 1]';
end


if ~is_fisheye
    % 假设patch内的z分量一样
    xyz_du_ref1 = xyz_du_ref * xyz_ref(3) / xyz_du_ref(3);
    xyz_dv_ref1 = xyz_dv_ref * xyz_ref(3) / xyz_dv_ref(3);
    xyz_du_ref = xyz_du_ref1;
    xyz_dv_ref = xyz_dv_ref1;
else
    % 假设patch内的xyz的模长一样
    xyz_du_ref0 = xyz_du_ref./norm(xyz_du_ref);
    xyz_dv_ref0 = xyz_dv_ref./norm(xyz_dv_ref);
    xyz_du_ref2 = xyz_du_ref0 * depth_ref_norm;
    xyz_dv_ref2 = xyz_dv_ref0 * depth_ref_norm;
    xyz_du_ref = xyz_du_ref2;
    xyz_dv_ref = xyz_dv_ref2;
end

px_cur = pflat(intrMat * (T_th(1:3,1:3)*xyz_ref' + T_th(1:3,4)));
px_du = pflat(intrMat * (T_th(1:3,1:3)*xyz_du_ref + T_th(1:3,4)));
px_dv = pflat(intrMat * (T_th(1:3,1:3)*xyz_dv_ref + T_th(1:3,4)));
A_th = zeros(2,2);
A_th(:,1) = (px_du(1:2) - px_cur(1:2))/(halfpatch_size+1);
A_th(:,2) = (px_dv(1:2) - px_cur(1:2))/(halfpatch_size+1);
end