function testH22()

check_id = 36;host_id = 1;target_id = 2;load('G:\matlab\LK\slam_algo\patch_match.mat');


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

host_pix = host_px(check_id,:);
target_pix = target_proj(check_id,:);

host_patch = host_pix + patch_offset;

host_ind = sub2ind(size(host_img), round(host_patch(:,2)), round(host_patch(:,1)));
host_value = double(host_img(host_ind));

end