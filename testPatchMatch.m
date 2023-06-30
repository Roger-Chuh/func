function testPatchMatch()
close all;

inputDir = 'G:\matlab\data\direct\sim';

ref_patch = imread(fullfile(inputDir, 'ref_patch.png'));
cur_patch = imread(fullfile(inputDir, 'cur_patch.png'));
ref_pt = load(fullfile(inputDir, 'host.txt'));
ref_pt = reshape(ref_pt, 2, [])' + 1;
ref_pt_x = reshape(ref_pt(:,1),10,10);
ref_pt_y = reshape(ref_pt(:,2),10,10);
ref_px = [ref_pt_x(6,6) ref_pt_y(6,6)];
ref_offset = ref_pt - ref_px;
ref_pt_orig = ref_offset + (ref_px-0.5).*2+0.5;
cur_pt = load(fullfile(inputDir, 'target.txt'));
cur_pt = reshape(cur_pt, 2, [])' + 1;
cur_pt_x = reshape(cur_pt(:,1),8,8);
cur_pt_y = reshape(cur_pt(:,2),8,8);
cur_px = [cur_pt_x(5,5) cur_pt_y(5,5)];
cur_offset = cur_pt - cur_px;

figure,subplot(2,2,1);imshow(ref_patch);hold on;plot(ref_pt(:,1), ref_pt(:,2), '.r');plot(ref_px(1),ref_px(2), '.g');title(sprintf('ref img size: [%d %d]', size(ref_patch,2), size(ref_patch,1)));
subplot(2,2,2);imshow(cur_patch);hold on;plot(cur_pt(:,1), cur_pt(:,2), '.r');plot(cur_px(1),cur_px(2), '.g');title(sprintf('cur img size: [%d %d]', size(cur_patch,2), size(cur_patch,1)));
subplot(2,2,3);imshow(imresize(ref_patch, [480 640]));hold on;plot(ref_pt_orig(:,1), ref_pt_orig(:,2), '.r');plot( (ref_px(1)-0.5).*2+0.5, (ref_px(2)-0.5).*2+0.5, '.g');title(sprintf('ref img size: [%d %d]', 640, 480));
subplot(2,2,4);imshow(imresize(cur_patch, [480 640]));




% check_id = 39;host_id = 1;target_id = 8;load('G:\matlab\LK\slam_algo\patch_match.mat');
check_id = 39;host_id = 1;target_id = 2;load('G:\matlab\LK\slam_algo\patch_match2.mat');
% check_id = 43;host_id = 1;target_id = 2; load('G:\matlab\LK\slam_algo\patch_match3.mat');
% load('G:\matlab\LK\slam_algo\pba_large_scene_2.mat');




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
host_px = host_pt.Location;
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
level_ref = 2; 1;2;0;1;
halfpatch_size = 4;
depth = depths(check_id);
depth_ref = depth;
xyz_ref = host_xyz(check_id, :);
depth_ref_norm = norm(xyz_ref);
% warped_gx = interp2(double(gx1), pt1s(:, 1), pt1s(:, 2), 'linear', 0);
% warped_gy = interp2(double(gy1), pt1s(:, 1), pt1s(:, 2), 'linear', 0);
px_ref = host_px(check_id,:);
[A_th, px_cur] = GetWarpMatrixAffine(intrMat, px_ref, level_ref, halfpatch_size, T_th, xyz_ref);


search_level = GetBestSearchLevel(A_th, max_levels - 1);


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
WarpAffine(A_th, host_img, px_ref, level_ref, target_img, pc_cur_check, search_level, halfpatch_size+1);

% is_host_square = true;
% WarpAffine(inv(A_th), target_img, pc_cur_check, search_level, host_img, px_ref, level_ref, halfpatch_size+1);
[A_th_scaled, ~] = GetWarpMatrixAffineInv(intrMat, px_ref, level_ref, halfpatch_size, T_th, xyz_ref);
search_level_scaled = GetBestSearchLevel(inv(A_th_scaled), max_levels - 1);
WarpAffineInv(A_th, host_img, px_ref, level_ref, target_img, pc_cur_check, search_level, halfpatch_size+1);


end
function search_level = GetBestSearchLevel(A_cur_ref, max_level)
search_level = 0;
D0 = det(A_cur_ref);
D = det(A_cur_ref);
while(D > 2.0 && search_level < max_level)
    search_level = search_level + 1;
    D = D * 0.25;
end

end
function WarpAffine(A_cur_ref, host_img, px_ref, level_ref, target_img, px_cur, search_level, halfpatch_size)
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


end
function WarpAffineInv(A_cur_ref, host_img, px_ref, level_ref, target_img, px_cur, search_level, halfpatch_size)
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
function [A_th, px_cur] = GetWarpMatrixAffine(intrMat, px_ref, level_ref, halfpatch_size, T_th, xyz_ref)

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

px_cur = pflat(intrMat * (T_th(1:3,1:3)*xyz_ref' + T_th(1:3,4)));
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