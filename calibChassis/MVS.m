function MVS()

im=double(imread('D:\Work\project1_IM_2018-2019-master\project1_IM_2018-2019-master\src\21\test0000.jpg'));
[H,W,~]=size(im);

cameras=fopen('D:\Work\project1_IM_2018-2019-master\project1_IM_2018-2019-master\src\21\cameras2.txt','r');

camera_mat=fscanf(cameras,'%f %f %f',[3,Inf]);

dMin = 10;
dMax = 100;

inputDir = 'D:\Work\project1_IM_2018-2019-master\project1_IM_2018-2019-master\src\21';

dirInfo = dir(fullfile(inputDir, '*.jpg'));
imageNum = 5;
step = 1; 2; 1; 3;
base = 200;0;200; 303;
cnt = 1;

fastCalc = 1;

for i = 1 : step : 1+step*(imageNum-1)
    if ~fastCalc
        image = imresize(rgb2gray(imread(fullfile(inputDir, dirInfo(base+i).name))),1);
    else
        image = imresize(rgb2gray(imread(fullfile(inputDir, dirInfo(base+i).name))),0.5);
    end
    img0.img = image;
    img1.img = imresize(image, 0.5);
    img2.img = imresize(image, 0.25);
    img0.dMin = dMin;  img0.dMax = dMax;
    img1.dMin = dMin;  img1.dMax = dMax;
    img2.dMin = dMin;  img2.dMax = dMax;
    img0.ID = i;
    img1.ID = i;
    img2.ID = i;
    %         data{i,1} = img;
    seq = (base+i-1)*7;
    K_cen = camera_mat(:,1+seq:3+seq)';
    
    R_cen = camera_mat(:,4+seq:6+seq)';
    T_cen = camera_mat(:,7+seq);
    if ~fastCalc
        intrMat0 = K_cen;
    else
        intrMat0 = [K_cen(1,1)/2 0 (1+K_cen(1,3))/2; 0 K_cen(2,2)/2 (1+K_cen(2,3))/2;0 0 1];
    end
    intrMat1 = [intrMat0(1,1)/2 0 (1+intrMat0(1,3))/2; 0 intrMat0(2,2)/2 (1+intrMat0(2,3))/2;0 0 1];
    intrMat2 = [intrMat0(1,1)/4 0 (1+intrMat0(1,3))/4; 0 intrMat0(2,2)/4 (1+intrMat0(2,3))/4;0 0 1];
    
    img0.intrMat = intrMat0;
    img1.intrMat = intrMat1;
    img2.intrMat = intrMat2;
    
    
    pose = [R_cen T_cen; 0 0 0 1];
    pose_inv = inv(pose);
    
    img0.P = intrMat0*pose(1:3,:);
    img1.P = intrMat1*pose(1:3,:);
    img2.P = intrMat2*pose(1:3,:);
    
    img0.R = pose(1:3,1:3);
    img1.R = pose(1:3,1:3);
    img2.R = pose(1:3,1:3);
    
    img0.T = pose(1:3,4);
    img1.T = pose(1:3,4);
    img2.T = pose(1:3,4);
    
    img0.C = pose_inv(1:3,4);
    img1.C = pose_inv(1:3,4);
    img2.C = pose_inv(1:3,4);
    
    %         img.l2g_l2g = [R_cen' T_cen; 0 0 0 1];
    %         g2l_g2l = inv([R_cen T_cen; 0 0 0 1]);
    
    data{i,1} = img0;
    data{i,2} = img1;
    data{i,3} = img2;
    data{i,4} = intrMat0;
    data{i,5} = intrMat1;
    data{i,6} = intrMat2;
    
    data{i,7} = pose;
    %         data{i,5} = [R_cen' T_cen; 0 0 0 1];%rot L2G, trans: L2G
    %         data{i,6} = inv([R_cen' T_cen; 0 0 0 1]); %rotG2L, trans: G2L
    
    
    if 0
        imageT =  image';
        calibRomFid = fopen(fullfile(inputDir, sprintf('img_%04d.bin', cnt)),'w');
        fwrite(calibRomFid, im2double(imageT(:)),'double');
        fclose(calibRomFid);
        
        poseT = inv(pose)';
        poseT(4,1:3) = poseT(4,1:3)./100;
        calibRomFid = fopen(fullfile(inputDir, sprintf('pose_%04d.bin', cnt)),'w');
        fwrite(calibRomFid, poseT(:),'single');
        fclose(calibRomFid);
        
        intrMatT = intrMat0';
        calibRomFid = fopen(fullfile(inputDir, sprintf('intrMat_%04d.bin', 1)),'w');
        fwrite(calibRomFid, intrMatT(:),'double');
        fclose(calibRomFid);
        cnt = cnt+1;
    end
end

%     return;

%% main

downSample1 = data(:,2);
downSample2 = data(:,3);
downSample2 = Init(downSample2);

geo_consistency = true;
for i = 1 : imageNum
    weightMapArr{i,1} = zeros(size(downSample2{1, 1}.img));
end
for i = 1 :imageNum
    for j = 1 : imageNum
        weightMap{j,1} = zeros(size(downSample2{1, 1}.img));
    end
    [downSample2, weightMap] = ACMH(downSample2, i, weightMap, 0);
    
end
fprintf('finish ACMH at coarsest scale\n');
fprintf('start ACMH & geometry at coarsest scale\n');
for i = 1 :imageNum
    for j = 1 : imageNum
        weightMap{j,1} = zeros(size(downSample2{1, 1}.img));
    end
    [downSample2, weightMap] = ACMH(downSample2, i, weightMap, 1);
    weightMapArr{i} = weightMap;
end
fprintf('finish ACMH & geometry at coarsest scale\n');


for i = 1 : imageNum
    high = downSample1{i};
    low = downSample2{i};
    %         high.depth = zeros(size(high.img));
    high.depth = JointBilateralUpsample2(high.img, low.img, low.depth);
    high.normalMap = JointBilateralUpsample2(high.img, low.img, low.normalMap);
    downSample1{i}.depth =  high.depth;
    downSample1{i}.normalMap =  high.normalMap;
end
% downSample1 =  high;
fprintf('finish up sample at coarsest scale\n');

for i = 1 : imageNum
    %         DetailRestore(downSample1, i, weightMapArr{i});
    [initCostMap, refImage] = DetailRestore(downSample1, i, weightMapArr{i});
    downSample1{i} = refImage;
end
fprintf('finish detail restoring at coarsest scale\n');








%% 中间尺度的深度计算
for i = 1 : imageNum
    weightMapArr{i,1} = zeros(size(downSample1{1, 1}.img));
end

for i = 1: imageNum
    
    %          WeightMap weightMap(imageNum);
    for j = 1: imageNum
        weightMap{j,1} = zeros(size(downSample1{1}.img));
    end
    [downSample1, weightMap] = ACMH(downSample1, i, weightMap, true);
    weightMapArr{i} = weightMap;
    
end

% 提升一个尺度
downSample0 = data(:,1);
for i = 1 : imageNum
    high = downSample0{i};
    low = downSample1{i};
    high.depth = JointBilateralUpsample2(high.img, low.img, low.depth);
    high.normalMap = JointBilateralUpsample2(high.img, low.img, low.normalMap);
    downSample0{i}.depth =  high.depth;
    downSample0{i}.normalMap =  high.normalMap;
end
% downSample0 =  high;


for i = 1 : imageNum
    %         DetailRestore(downSample1, i, weightMapArr{i});
    [initCostMap, refImage] = DetailRestore(downSample0, i, weightMapArr{i});
    downSample0{i} = refImage;
end



%% 原始分辨率的深度计算
for i = 1 : imageNum
    weightMapArr{i,1} = zeros(size(downSample1{1, 1}.img));
end

for i = 1 : imageNum
    
    %             WeightMap weightMap(imageNum);
    for j = 1 : imageNum
        weightMap{j,1} = zeros(size(downSample0));
    end
    [downSample0, weightMap] = ACMH(downSample0, i, weightMap, true);
    weightMapArr{i} = weightMap;
    
    
end

fprintf('finish ACMH & geometry at raw scale\n');



end
function data1 = Init(data)
imageNum = length(data);
image = data{1,1};
% % [x,y,z] = sphere(size(image.img,1) * size(image.img,2));

[xGrid, yGrid] = meshgrid(1:size(image.img,2), 1:size(image.img,1));
pix = [xGrid(:) yGrid(:)];
pt3d = inv(image.intrMat)*pextend(pix');
for i = 1 : imageNum
    image = data{i,1};
    image.depth = image.dMin + (image.dMax-image.dMin).*rand(size(image.img));
    x= -1 + (1-(-1)).*rand(size(image.img,1)*size(image.img,2), 1);
    y= -1 + (1-(-1)).*rand(size(image.img,1)*size(image.img,2), 1);
    z= -1 + (1-(-1)).*rand(size(image.img,1)*size(image.img,2), 1);
    [vec,~] = NormalizeVector([x y z]);
    sign = dot(pt3d, vec');
    id = find(sign>0);
    vec(id,:) = -vec(id,:);
    image.normalMap = zeros([size(image.img) 3]);
    image.normalMap(:,:,1) = reshape(vec(:,1), size(image.img));
    image.normalMap(:,:,2) = reshape(vec(:,2), size(image.img));
    image.normalMap(:,:,3) = reshape(vec(:,3), size(image.img));
    data1{i,1} = image;
end





end
function [images, weightMapArr] = ACMH(images, refId, weightMapArr, geo_consistency)

imageNum = length(images);
w = size(images{1}.img, 2);
h = size(images{1}.img, 1);


maxIteration = 5;
fprintf(sprintf('Iteration times: %d, reference Id: %d\n',maxIteration,refId));
refImage = images{refId,1};

refImage.confMap = 2.*ones(h, w);

fprintf(sprintf('Current scale is %d x %d\n',size(refImage.img,2), size(refImage.img,1)));
lastImportantView = inf(h, w);

for it = 1 : maxIteration
    
    initCost = InitMatchingCost(images, refId);
    fprintf( 'Finish init matching cost\n');
    
    for a = 0 : 1
        
        for i = 0 : h-1
            
            for j = mod(mod(i,2)+a,2):2:w-1
                pt = [j+1, i+1];
                
                good = CheckerboardSampling(refImage, pt, initCost);
                good_depth_normal = [];
                for ii = 1 : size(good,1)
                    point = good(ii,:);
                    d = refImage.depth(point(:,2), point(:,1));
                    n = [refImage.normalMap(point(:,2), point(:,1),1); refImage.normalMap(point(:,2), point(:,1),2); refImage.normalMap(point(:,2), point(:,1),3)];
                    depth_normal = [d; n];
                    good_depth_normal = [good_depth_normal; depth_normal'];
                end
                bilateralCost = ComputeBilateralNCC(images, good_depth_normal, pt, refId);
                
                if geo_consistency
                    projError = ComputeReprojectionError(images, good_depth_normal, refId, pt);
                    [viewWeight, lastId]= ComputeViewWeight(bilateralCost, it, refId, lastImportantView(pt(:,2), pt(:,1)), false);
                    if 0 %因为updateViewWeight的flag被置成false，所以不用更新lastImportantView
                        lastImportantView(pt(:,2), pt(:,1)) = lastId;
                    end
                    idx_cost = SelectHypotheses(bilateralCost, viewWeight, geo_consistency, projError);
                else
                    [viewWeight, lastId]= ComputeViewWeight(bilateralCost, it, refId, lastImportantView(pt(:,2), pt(:,1)), false);
                    if 0 %因为updateViewWeight的flag被置成false，所以不用更新lastImportantView
                        lastImportantView(pt(:,2), pt(:,1)) = lastId;
                    end
                    idx_cost = SelectHypotheses(bilateralCost, viewWeight, geo_consistency, []);
                end
                
                if idx_cost(2) < refImage.confMap(pt(:,2), pt(1))
                    refImage.confMap(pt(:,2), pt(1)) = idx_cost(2);
                    refImage.depth(pt(:,2), pt(1)) = good_depth_normal(idx_cost(1),1);
                    refImage.normalMap(pt(:,2), pt(1),1) = good_depth_normal(idx_cost(1),2);
                    refImage.normalMap(pt(:,2), pt(1),2) = good_depth_normal(idx_cost(1),3);
                    refImage.normalMap(pt(:,2), pt(1),3) = good_depth_normal(idx_cost(1),4);
                end
            end
        end
    end
    
    
    for i = 1 : h
        for j = 1 : w
            pt = [j ,i];
            depth = refImage.depth(pt(:,2), pt(:,1));
            normal = [refImage.normalMap(pt(2), pt(1), 1) refImage.normalMap(pt(2), pt(1), 2) refImage.normalMap(pt(2), pt(1), 3)]';
            refine_depth_normal = GenerateDepthNormal(depth, normal, it, refImage, pt);
            bilateralCost = ComputeBilateralNCC(images, refine_depth_normal, pt, refId);
            
            
            if(geo_consistency)
                
                projError = ComputeReprojectionError(images, refine_depth_normal, refId, pt);
                [viewWeight, lastId] = ComputeViewWeight(bilateralCost, it, refId, lastImportantView(pt(2),pt(1)), true);
                lastImportantView(pt(:,2), pt(:,1)) = lastId;
                idx_cost = SelectHypotheses(bilateralCost, viewWeight, geo_consistency, projError);
                
            else
                
                [viewWeight, lastId] = ComputeViewWeight(bilateralCost, it, refId, lastImportantView(pt(2),pt(1)), true);
                lastImportantView(pt(:,2), pt(:,1)) = lastId;
                idx_cost = SelectHypotheses(bilateralCost, viewWeight, geo_consistency, []);
            end
            
            
            assert(length(weightMapArr) == length(viewWeight));
            for k = 1 : imageNum
                weightMapArr{k}(pt(2),pt(1)) = viewWeight(k);
            end
            
            if idx_cost(2) < refImage.confMap(pt(:,2), pt(1))
                refImage.confMap(pt(:,2), pt(1)) = idx_cost(2);
                refImage.depth(pt(:,2), pt(1)) = refine_depth_normal(idx_cost(1),1);
                refImage.normalMap(pt(:,2), pt(1),1) = refine_depth_normal(idx_cost(1),2);
                refImage.normalMap(pt(:,2), pt(1),2) = refine_depth_normal(idx_cost(1),3);
                refImage.normalMap(pt(:,2), pt(1),3) = refine_depth_normal(idx_cost(1),4);
            end
            
        end
    end
    
end


depthFilt = refImage.depth;
depthFilt_ = medfilt2(depthFilt, [5 5]);
refImage.depth = depthFilt_;
images{refId} = refImage;
end

function initCost = InitMatchingCost(images, refId)

w = size(images{1}.img, 2);
h = size(images{1}.img, 1);


imageNum = length(images);
refImage = images{refId,1};
initCost = zeros(size(refImage.img));
for i = 1 : imageNum
    img = images{i};
    neibor.ID = img.ID;
    neibor.Hl = img.intrMat * img.R * refImage.R';
    neibor.Hm = img.intrMat * img.R * (refImage.C - img.C);
    neibor.Hr = inv(refImage.intrMat);
    imageInfos{i,1} = neibor;
end

for i = 1 : w
    
    for j = 1 : h
        
        pt = [i j];
        normal = [refImage.normalMap(pt(2),pt(1),1) refImage.normalMap(pt(2),pt(1),2) refImage.normalMap(pt(2),pt(1),3)]';
        depth = refImage.depth(pt(2),pt(1));
        
        costs = [];
        for m = 1 : imageNum
            srcImage = images{m};
            
            
            if (m == refId)
                
                costs = [costs; 2];
                continue;
            end
            
            
            
            cost = BilateralNCC(refImage, srcImage, imageInfos{m}, pt, depth, normal);
            costs = [costs; cost];
        end
        K = 5;
        [~, minIndex] = sort(costs,'ascend');
        cost = mean(costs(minIndex(1:K)));
        initCost(pt(2), pt(1)) = cost;
    end
    
    
end


end
function cost = BilateralNCC(refImage, srcImage, imageInfo, pix, depth, normal)
sigma_g = 0.2;
sigma_x = 3;
nSizeHalfWindow = 3;
w = size(refImage.img, 2);
h = size(refImage.img, 1);
% [xGrid, yGrid] = meshgrid(1:size(refImage.img,2), 1:size(refImage.img,1));
nSizeWindow = 2 * nSizeHalfWindow + 1;
nTexels = nSizeWindow * nSizeWindow;
lt0 = [pix(1) - nSizeHalfWindow pix(2) - nSizeHalfWindow];
rb0 = [pix(1) + nSizeHalfWindow pix(2) + nSizeHalfWindow];
if ~isInside(lt0, w, h) || ~isInside(rb0, w, h)
    cost = 2;
    return;
end

center = double(refImage.img(pix(2), pix(1)));

k=1;
for i = 1 : nSizeWindow
    for j = 1 : nSizeWindow
        pix_ = [lt0(1) + j lt0(2) + i];
        current = double(refImage.img(pix_(2), pix_(1)));
        delta_g = (current - center) / 255;
        ref_texels(k) = current;
        ref_texels_square(k) = current * current;
        distance_square(k) = (i - nSizeHalfWindow) * (i - nSizeHalfWindow) + (j - nSizeHalfWindow) * (j - nSizeHalfWindow);
        bilateral_weight(k) = exp(- delta_g * delta_g / (2 * sigma_g * sigma_g) - distance_square(k) / (2 * sigma_x * sigma_x)); %双边权重
        k = k+1;
    end
    
end

X0 = inv(refImage.intrMat)*[pix 1]';
H = (imageInfo.Hl + (imageInfo.Hm) * normal' * (1 / dot(normal ,  X0 * depth))) * imageInfo.Hr;
k = 1;
for i = 1 : nSizeWindow
    for j = 1 : nSizeWindow
        
        pt = ProjectH(H, [lt0(1)+j  lt0(2)+i]);
        if ~isInside(pt, w, h)
            cost = 2;
            return;
        end
        
        src_texels(k) = Sample(srcImage.img, pt); % 通过双线性插值计算灰度值
        src_texels_square(k) = src_texels(k) * src_texels(k);
        k=k+1;
        
    end
end



bilateral_weight_sum = sum(bilateral_weight);
src_ref_avg = sum(ref_texels .* src_texels .* bilateral_weight) ./ bilateral_weight_sum;
src_avg = sum(src_texels .* bilateral_weight) ./ bilateral_weight_sum;
ref_avg = sum(ref_texels .* bilateral_weight) ./ bilateral_weight_sum;
src_avg_square = sum(src_texels_square .* bilateral_weight) ./ bilateral_weight_sum;
ref_avg_square = sum(ref_texels_square .* bilateral_weight) ./ bilateral_weight_sum;
src_src_cov = src_avg_square - src_avg .* src_avg;
ref_ref_cov = ref_avg_square - ref_avg .* ref_avg;

kMinVar = 1e-5;
if (src_src_cov < kMinVar || ref_ref_cov < kMinVar)
    cost = 2;
    return;
end

src_ref_cov = src_ref_avg - src_avg .* ref_avg;
ncc = src_ref_cov ./ sqrt(src_src_cov .* ref_ref_cov);
cost = max(0, min(2, 1 - ncc));
end
function flag = isInside(pix, w, h)
flag = pix(1) >= 1 && pix(2) >= 1 && pix(1) < w && pix(2) + 1 < h;
end
function proj = ProjectH(H, X)
HH = H';
invZ = 1 / (HH(7) * X(1) + HH(8) * X(2) + HH(9));
proj = [(HH(1) * X(1) + HH(2) * X(2) + HH(3))*invZ, (HH(4) * X(1) + HH(5) * X(2) + HH(6))*invZ];
end

function interp = Sample(imageGray, pix)

imageGray = double(imageGray);

lx = floor(pix(1));
ly = floor(pix(2));
x = pix(1) - lx;
y = pix(2) - ly;
x1 = 1 - x;
y1 = 1 - y;
pt1 = [lx, ly];
pt2 = [lx + 1, ly];
pt3 = [lx, ly + 1];
pt4 = [lx + 1, ly + 1];
interp =  (imageGray(pt1(2),pt1(1)) * x1 + imageGray(pt2(:,2),pt2(:,1)) * x) * y1 + (imageGray(pt3(:,2), pt3(:,1)) * x1 + imageGray(pt4(:,2),pt4(:,1)) * x) * y;
end
function Candidate = CheckerboardSampling(refImage, pt, costMap)
w = size(refImage.img, 2);
h = size(refImage.img, 1);
candidate = [pt(1)-1 pt(2); pt(1)+1 pt(2); pt(1) pt(2)-1;pt(1) pt(2)+1];


for i = 2:4
    candidate = [candidate; [pt(1)-i pt(2)+i-1]];
    candidate = [candidate; [pt(1)-i pt(2)-i+1]];
    
    candidate = [candidate; [pt(1)+i pt(2)+i-1]];
    candidate = [candidate; [pt(1)+i pt(2)-i+1]];
    
    candidate = [candidate; [pt(1)+i-1 pt(2)-i]];
    candidate = [candidate; [pt(1)-i+1 pt(2)-i]];
    
    candidate = [candidate; [pt(1)+i-1 pt(2)+i]];
    candidate = [candidate; [pt(1)-i+1 pt(2)+i]];
    
end



for i = 3:2:24
    
    candidate = [candidate; pt(1), pt(2) - i]; % 上
    candidate = [candidate; pt(1), pt(2) + i]; % 下
    candidate = [candidate; pt(1) - i, pt(2)]; % 左
    candidate = [candidate; pt(1) + i, pt(2)]; % 右
end


hypothesis = {};
Score = [];
for i = 1 : size(candidate, 1)
    
    point = candidate(i,:);
    if ~isInside(point,w,h)
        continue;
    end
    ps.pt = point;
    ps.score = costMap(point(2), costMap(1));
    hypothesis = [hypothesis; {ps}];
    Score = [Score; costMap(point(2), costMap(1))];
end

[~, ind] = sort(Score,'ascend');
hypothesis = hypothesis(ind,:);
hypothesis = hypothesis(1:min(8,length(hypothesis)),:);
Candidate = [];

for i = 1 : length(hypothesis)
    Candidate = [Candidate; hypothesis{i}.pt];
end
end

function bilateralCost = ComputeBilateralNCC(images, depth_normal, pt, refID)
imageNum = length(images);
bilateralCost = zeros(size(depth_normal,1), imageNum);
refImage = images{refID};

imageInfos = {};
for i = 1 : imageNum
    image = images{i};
    neibor.ID = image.ID;
    neibor.Hl = image.intrMat * image.R * refImage.R';
    neibor.Hm = image.intrMat * image.R * (refImage.C - image.C);
    neibor.Hr = inv(refImage.intrMat);
    imageInfos = [imageInfos; {neibor}];
end

for i = 1 : size(bilateralCost, 1)
    
    depth =  depth_normal(i,1);
    normal = depth_normal(i,2:4)';
    for j = 1 : size(bilateralCost, 2)
        if refID == j
            
            bilateralCost(i, j) = 2;
            continue;
        end
        srcImage = images{j};
        imageInfo = imageInfos{j};
        bilateralCost(i,j) =  BilateralNCC(refImage, srcImage, imageInfo, pt, depth, normal);
    end
    
end

end
function  projError = ComputeReprojectionError(images, depth_normal, refID, pt)
refImage = images{refID};
w = size(refImage.img,2);
h = size(refImage.img,1);
projError = zeros(size(depth_normal,1), length(images));

for m = 1 : size(projError, 1)
    for n = 1 : size(projError, 2)
        if n == refID
            projError(m,n) = 3;
            continue;
        end
        d = depth_normal(m,1);
        srcImage = images{n};
        X = TransformPointI2W(pt(1), pt(2), d, refImage.intrMat);
        srcX = TransformPointC2Ii(TransformPointW2C(X, srcImage.R, srcImage.C ), srcImage.intrMat);
        
        if ~isInside(srcX, w, h)
            projError(m, n) = 3;
            continue;
        end
        
        X = TransformPointI2W(srcX(1), srcX(2), srcImage.depth(round(srcX(:,2), round(srcX(:,1)))), srcImage.intrMat);
        refX = TransformPointC2Ii(TransformPointW2C(X, refImage.R, refImage.C ), refImage.intrMat);
        if ~isInside(refX, w, h)
            projError(m, n) = 3;
            continue;
        end
        error = sqrt((pt(1)-refX(1)).^2 + (pt(2)-refX(2)).^2);
        projError(m, n) = min([error,3]);
    end
end

end
function xyz = TransformPointI2W(x, y, d, K)

xyz = [(x - K(1, 3)) * d / K(1, 1)	(y - K(2, 3)) * d / K(2, 2) 	d]';

end
function pix = TransformPointC2Ii(xyz, K)
x = xyz(1);
y = xyz(2);
z = xyz(3);
pix = [round(K(1, 3) + K(1, 1) * x / z),		round(K(2, 3) + K(2, 2) * y / z)];

end
function pt3d = TransformPointW2C(X,R,C)
pt3d = R * (X - C);
end
function [viewWeight, last] = ComputeViewWeight(cost, iteration, refID, last, varargin)


if (nargin == 4)
    update = 0;
elseif (nargin == 5)
    update = varargin{1};
else
    
end




init_good_threshold = 0.8;
alpha = 90;
beta = 0.3;
good_threshold = init_good_threshold * exp(-iteration * iteration / alpha);
bad_threshold = 1.2;
n1 = 2;
n2 = 3;
viewWeight = zeros(size(cost,2),1);
S_t = [];

for i = 1 : size(cost, 2)
    if i == refID
        viewWeight(i) = 0;
    end
    
    S_good = [];
    S_bad = [];
    
    for j = 1 : size(cost, 1)
        c = cost(j, i);
        if c < good_threshold
            S_good = [S_good; c];
        elseif c > bad_threshold
            S_bad = [S_bad; c];
        end
        
    end
    
    if length(S_good) > n1 && length(S_bad) < n2
        S_t = [S_t; i];
    else
        viewWeight(i) = 0.2 * (i == last);
        continue;
    end
    confidence = zeros(length(S_good), 1);
    
    
    for ii = 1: length(confidence)
        c = S_good(ii);
        conf = exp(-c * c / (2 * beta * beta));
        confidence(ii) = conf;
    end
    weight = sum(confidence) / length(S_good);
    weight = ((i == last) + 1) * weight;
    viewWeight(i) = weight;
end

if update
    
    [~,id] = max(viewWeight);
    last = id;
    
end

end
function idx_cost = SelectHypotheses(cost, viewWeight, varargin)




if (nargin == 2)
    geo_consistency = 0;
    projError = [];
else
    geo_consistency = varargin{1};
    projError = varargin{2};
end






viewCost = zeros(size(cost, 2), 1);
photometric = [];


for i = 1 : size(cost, 1)
    
    for j = 1 : size(cost, 2)
        if geo_consistency
            viewCost(j) = cost(i,j) + 0.2*projError(i,j);
        else
            viewCost(j) = cost(i,j);
        end
    end
    
    photometric = [photometric; sum(viewCost .* viewWeight)./sum(viewWeight)];
    
end
[c, idx] = min(photometric);
idx_cost = [idx c];
end
function depth_normal = GenerateDepthNormal(depth, normal, iteration, refImage, pt)

% normal is 3 by 1 vector

pt3d = inv(refImage.intrMat)*pextend(pt');
rand_depth = refImage.dMin + (refImage.dMax-refImage.dMin).*rand(1,1);
x= -1 + (1-(-1)).*rand(1, 1);
y= -1 + (1-(-1)).*rand(1, 1);
z= -1 + (1-(-1)).*rand(1, 1);
[vec,~] = NormalizeVector([x y z]);
sign = dot(pt3d, vec');
id = find(sign>0);
vec(id,:) = -vec(id,:);
rand_normal = vec';
perturbation = 1./(2.^iteration);
max_depth = (1 + perturbation) .* depth;
min_depth = (1 - perturbation) .* depth;
prt_depth = (0 + (1-0).*rand(1,1)).*(max_depth - min_depth) + min_depth;
prt_normal = PerturbNormal(normal, perturbation * pi);
view_ray = [(pt(1) - refImage.intrMat(1, 3)) / refImage.intrMat(1, 1);(pt(2) - refImage.intrMat(2, 3)) / refImage.intrMat(2, 2); 1];

m = 0;
while dot(prt_normal, view_ray) > 0
    if m == 3
        prt_normal = normal;
        break;
    end
    prt_normal = PerturbNormal(normal, 0.5 * perturbation * pi);
    perturbation  =  0.5 * perturbation;
    m = m+1;
end
prt_normal = prt_normal./norm(prt_normal);

dn1 = [depth normal'];
dn2 = [depth, prt_normal'];
dn3 = [depth, rand_normal'];
dn4 = [rand_depth, normal'];
dn5 = [rand_depth, prt_normal'];
dn6 = [rand_depth, rand_normal'];
dn7 = [prt_depth, normal'];
dn8 = [prt_depth, rand_normal'];
dn9 = [prt_depth, prt_normal'];
depth_normal = [dn1; dn2; dn3; dn4; dn5; dn6; dn7; dn8; dn9];
end


function  normal = GenerateRandomNormal(refImage, pt)
pt3d = inv(refImage.intrMat)*pextend(pt');
x= -1 + (1-(-1)).*rand(1, 1);
y= -1 + (1-(-1)).*rand(1, 1);
z= -1 + (1-(-1)).*rand(1, 1);
[vec,~] = NormalizeVector([x y z]);
sign = dot(pt3d, vec');
id = find(sign>0);
vec(id,:) = -vec(id,:);
normal = vec';


end
function  prt_normal = PerturbNormal( normal,  perturbation)

a1 = (rand(1,1)-0.5)*perturbation;
a2 = (rand(1,1)-0.5)*perturbation;
a3 = (rand(1,1)-0.5)*perturbation;
sin_a1 = sin(a1);
sin_a2 = sin(a2);
sin_a3 = sin(a3);
cos_a1 = cos(a1);
cos_a2 = cos(a2);
cos_a3 = cos(a3);

R(1,1) = cos_a2 * cos_a3;
R(1,2) = -cos_a2 * sin_a3;
R(1,3) = sin_a2;
R(2,1) = cos_a1 * sin_a3 + cos_a3 * sin_a1 * sin_a2;
R(2,2) = cos_a1 * cos_a3 - sin_a1 * sin_a2 * sin_a3;
R(2,3) = -cos_a2 * sin_a1;
R(3,1) = sin_a1 * sin_a3 - cos_a1 * cos_a3 * sin_a2;
R(3,2) = cos_a3 * sin_a1 + cos_a1 * sin_a2 * sin_a3;
R(3,3) = cos_a1 * cos_a2;

RR = R';
prt_normal = [RR(1) * normal(1) + RR(2) * normal(2) + RR(3) * normal(3);...
    RR(4) * normal(1) + RR(5) * normal(2) + RR(6) * normal(3);...
    RR(7) * normal(1) + RR(8) * normal(2) + RR(9) * normal(3)];
end
function upSampled = JointBilateralUpsample(high, low, coarseMap, varargin)


if ndims(high) < 3
    high = cat(3, high, high, high);
    low = cat(3, low, low, low);
end

width = size(high, 2);
height = size(high, 1);
factor = 2;  % 高像素与低像素图像尺度之比为2

if (nargin == 3)
    halfWindow = 5;
    sigma_d = 0.5;
    sigma_r = 0.1;
elseif (nargin == 4)
    halfWindow = varargin{1};
    sigma_d = 0.5;
    sigma_r = 0.1;
elseif (nargin == 5)
    halfWindow = varargin{1};
    sigma_d =  varargin{2};
    sigma_r = 0.1;
elseif (nargin == 6)
    halfWindow = varargin{1};
    sigma_d = varargin{2};
    sigma_r =  varargin{3};
else
    error('Too many input arguments');
end

type = 1;
if ndims(coarseMap) < 3
    upSampled = zeros(size(high));
else
    type = 0;
    upSampled = zeros([size(high) 3]);
end

for i = 1:height
    for j = 1 : width
        p = high(i,j,:);
        low_i = i / factor;
        low_j = j / factor;
        iMax = floor(min(size(low,1) - 0, low_i + halfWindow));
        iMin = ceil(max(1, low_i - halfWindow));
        jMax = floor(min(size(low,2) - 0, low_j + halfWindow));
        jMin = ceil(max(1, low_j - halfWindow));
        
        lowWindow = coarseMap(iMin:iMax,jMin:jMax,:);
        if ndims(coarseMap) < 3
            assert(min(lowWindow(:)) > 0);
        end
        
        iw = (iMin:iMax) - low_i;
        jw = (jMin:jMax) - low_j;
        [mx,my] = meshgrid(jw,iw);
        spatial = exp( -(mx.^2 + my.^2) ./ (2*sigma_d.^2) );
        
        %         highWindow是像素P在高分辨率下的支撑窗口，如果低分辨率时窗口为5*5
        %         高分辨率下应该就是10*10，但公式中要求二者一致，因此高分辨率下窗口也是5*5
        %         这就需要对窗口的数据降采样，隔一行采样一行，隔一列采样一列
        highWindow = high(iMin * factor:factor:iMax * factor,jMin * factor:factor:jMax * factor,:);
        
        dR = highWindow(:,:,1) - p(1);
        dG = highWindow(:,:,2) - p(2);
        dB = highWindow(:,:,3) - p(3);
        range = exp( -(dR.^2 + dG.^2 + dB.^2) ./ (2*sigma_r.^2));
        
        
        spatial_range = spatial .* range;
        Kp = sum(spatial_range(:));
        
        if type
            depth = lowWindow.*spatial_range;
            depth = sum(depth(:))./Kp;
            upSampled(i,j) = depth;
            
        else
            
            Normal = lowWindow.*cat(3, spatial_range, spatial_range, spatial_range);
            normal = [sum(sum(Normal(:,:,1))); sum(sum(Normal(:,:,2))); sum(sum(Normal(:,:,3)))]./Kp;
            upSampled(i,j,:) = normal;
            
        end
        
    end
end

end
function upSampled = JointBilateralUpsample2(high, low, coarseMap, varargin)


if ndims(high) < 3
    high = double(cat(3, high, high, high));
    low = double(cat(3, low, low, low));
else
    high = double(high);
    low = double(low);
end

width = size(high, 2);
height = size(high, 1);
factor = 2;  % 高像素与低像素图像尺度之比为2

if (nargin == 3)
    halfWindow = 5;
    sigma_d = 0.5;
    sigma_r = 0.1;
elseif (nargin == 4)
    halfWindow = varargin{1};
    sigma_d = 0.5;
    sigma_r = 0.1;
elseif (nargin == 5)
    halfWindow = varargin{1};
    sigma_d =  varargin{2};
    sigma_r = 0.1;
elseif (nargin == 6)
    halfWindow = varargin{1};
    sigma_d = varargin{2};
    sigma_r =  varargin{3};
else
    error('Too many input arguments');
end

type = 1;
if ndims(coarseMap) < 3
    upSampled = zeros(size(high(:,:,1)));
else
    type = 0;
    upSampled = zeros([size(high)]);
end

for i = 1:height
    for j = 1 : width
        p = high(i,j,:);
        low_i = i / factor;
        low_j = j / factor;
        iMax = floor(min(size(low,1) - 0, low_i + halfWindow));
        iMin = ceil(max(1, low_i - halfWindow));
        jMax = floor(min(size(low,2) - 0, low_j + halfWindow));
        jMin = ceil(max(1, low_j - halfWindow));
        
        lowWindow = coarseMap(iMin:iMax,jMin:jMax,:);
        if ndims(coarseMap) < 3
%             assert(min(lowWindow(:)) > 0);
        end
        
        iw = (iMin:iMax) - low_i;
        jw = (jMin:jMax) - low_j;
        [mx,my] = meshgrid(jw,iw);
        spatial = exp( -(mx.^2 + my.^2) ./ (2*sigma_d.^2) );
        
        %         highWindow是像素P在高分辨率下的支撑窗口，如果低分辨率时窗口为5*5
        %         高分辨率下应该就是10*10，但公式中要求二者一致，因此高分辨率下窗口也是5*5
        %         这就需要对窗口的数据降采样，隔一行采样一行，隔一列采样一列
        highWindow = high(iMin * factor:factor:iMax * factor,jMin * factor:factor:jMax * factor,:);
        
        dR = highWindow(:,:,1) - p(1);
        dG = highWindow(:,:,2) - p(2);
        dB = highWindow(:,:,3) - p(3);
        range = exp( -(dR.^2 + dG.^2 + dB.^2) ./ (2*sigma_r.^2));
        
        
        spatial_range = spatial .* range;
        Kp = sum(spatial_range(:));
        
        if type
            depth = lowWindow.*spatial_range;
            depth = sum(depth(:))./Kp;
            upSampled(i,j) = depth;
            
        else
            
            Normal = lowWindow.*cat(3, spatial_range, spatial_range, spatial_range);
            normal = [sum(sum(Normal(:,:,1))); sum(sum(Normal(:,:,2))); sum(sum(Normal(:,:,3)))]./Kp;
            upSampled(i,j,:) = normal;
            
        end
        
    end
end

end

function [initCostMap, refImage] = DetailRestore(images, refID,  weightMap)
assert(sum(abs(round(size(weightMap{1})) - round(size(images{1}.depth)./2))) == 0);
refImage = images{refID};
initCostMap = zeros(size(refImage.depth));

height = size(refImage.img, 1);
width = size(refImage.img, 2);
for i = 1 : height
    for j = 1 : width
        pt = [j i];
        depth = refImage.depth(pt(:,2), pt(:,1));
        normal = squeeze(refImage.normalMap(pt(:,2),pt(:,1),:));
        good_depth_normal = [depth normal'];
        bilateralNCC = ComputeBilateralNCC(images, good_depth_normal, pt, refID);
        assert(size(bilateralNCC,1) == 1);
        
        
        %           weightSum = 0; a = 0;
        ptDown = ceil([j/2 i/2]);
        
        weightMat = [];
        for kk = 1 : length(weightMap)
            weightMat = [weightMat weightMap{kk}(ptDown(:,2), ptDown(:,1))];
        end
        a = sum(sum(bilateralNCC.*weightMat));
        weightSum = sum(weightMat);
        initPhotoCost = a / weightSum;
        initCostMap(pt(2), pt(1)) = initPhotoCost;
    end
end
depthMapTmp = refImage.depth;
normalMapTmp = refImage.normalMap;
confMapTmp = zeros(size(refImage.depth));

initCost = InitMatchingCost(images, refID);
for a = 0 : 1
    
    for i = 0 : height-1
        
        for j = mod(mod(i, 2) + a, 2) : 2 : width-1
            pt = [j+1, i+1];
            good = CheckerboardSampling(refImage, pt, initCost);
            
            good_depth_normal = [];
            for ii = 1 : size(good,1)
                point = good(ii,:);
                d = refImage.depth(round(point(:,2)), round(point(:,1)));
                n = squeeze(refImage.normalMap(round(point(:,2)), round(point(:,1)),:));
                depth_normal = [d, n'];
                good_depth_normal = [good_depth_normal; depth_normal];
            end
            bilateralCost = ComputeBilateralNCC(images, good_depth_normal, pt, refID);
            idx_cost = [];
            tmp = 10000000;
            viewWeight = ComputeViewWeight(bilateralCost, 0, refID, tmp, 0);
            if 0
                % 因为updateViewWeight的flag被置成false，所以不用更新lastImportantView
                lastImportantView(pt(:,2), pt(:,1)) = lastId;
            end
            idx_cost = SelectHypotheses(bilateralCost, viewWeight);
            depthMapTmp(pt(2), pt(1)) = good_depth_normal(idx_cost(1), 1);
            normalMapTmp(pt(:,2), pt(1),:) = good_depth_normal(idx_cost(1),2:4);
            confMapTmp(pt(:,2), pt(1)) = idx_cost(2);
        end
    end
end

for i = 1 : height
    for j = 1 : width
        pt = [j i];
        depth = refImage.depth(pt(:,2), pt(:,1));
        normal = squeeze(refImage.normalMap(pt(:,2), pt(:,1), :));
        refine_depth_normal = GenerateDepthNormal(depth, normal, 0, refImage, pt);
        refine_bilateralCost = ComputeBilateralNCC(images, refine_depth_normal, pt, refID);
        bilateralCost = ComputeBilateralNCC(images, refine_depth_normal, pt, refID);
        tmp = 1000000;
        viewWeight = ComputeViewWeight(bilateralCost, 0, refID, tmp, 0);
        idx_cost = SelectHypotheses(bilateralCost, viewWeight);
        if idx_cost(2) < confMapTmp(pt(:,2), pt(1))
            
            depthMapTmp(pt(:,2), pt(1)) = refine_depth_normal(idx_cost(1),1);
            normalMapTmp(pt(:,2), pt(1),:) = refine_depth_normal(idx_cost(1),2:4);
            confMapTmp(pt(:,2), pt(1)) = idx_cost(2);
        end
    end
end

for i = 1 : height
    for j = 1 : width
        pt = [j,i];
        if initCostMap(pt(:,2), pt(1)) - confMapTmp(pt(:,2), pt(1)) > 0.1            
            refImage.depth(pt(:,2), pt(1)) = depthMapTmp(pt(:,2), pt(1));
            refImage.normalMap(pt(:,2), pt(1),:) = normalMapTmp(pt(:,2), pt(1),:);
        end
    end
end

%   images{refID} = refImage;

end