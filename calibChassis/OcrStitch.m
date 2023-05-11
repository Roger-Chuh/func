function OcrStitch()
global mergeDist ovlpRatio stopThr

minOvlp = 0.3;0.5;
minEdgeSize = 1600;400;
doOpt2 = 0;
stopThr = 1.5; 0.2;1e-4;
maxErr =  2; 1.5; 1.5;  2; 1.5;2.00500008;
result = [];
inlierRati=0.80007000050060007008;
ovlpRatio = 0.45; 0.75; 0.7; -0.8;

mergeDist = 5; 7; 5;7; 5; 3; 5; 2; 3; 2; 1.1; 2; 3; 10; 2; 3; 1.3; 1; 1.3; 2; 3; 2.5; 3; 6; 2;1.1; 2; 1.5; 2; 3; 10; 20; 200;
inputDir = 'C:\Users\rongjiezhu\Documents\WXWork\1688852647983342\Cache\File\2021-03\Youdao-OCR-拼图Sample_0309\Youdao-OCR-拼图Sample_0309\4';
inputDir = 'C:\Users\rongjiezhu\Documents\WXWork\1688852647983342\Cache\File\2021-03\Youdao-OCR-拼图Sample_0309\Youdao-OCR-拼图Sample_0309\5';
% inputDir = 'C:\Users\rongjiezhu\Documents\WXWork\1688852647983342\Cache\File\2021-03\Youdao-OCR-拼图Sample_0309\Youdao-OCR-拼图Sample_0309\6';

new =1;0;
if ~new
    mergeDist = 20;
end
ovwt = 0; 1;
dirInfo = dir(fullfile(inputDir,'*.jpg'));
transRatio = 1;
imgRng = 50:10000;  45:100000;
imgRng = [imgRng(1) : min(imgRng(end), length(dirInfo))];

thetaList = [-9:3:9];  % degree
alphaList = [-90:30:90];  % degree
scaleList = [-45:1:45]; % pixel
scaleList = [-5:1:5]; % pixel

[thetaMat, alphaMat] = meshgrid(thetaList, alphaList);

thetaAlphaList = [thetaMat(:) alphaMat(:)];
thetaAlphaScaleList = repmat(thetaAlphaList, length(scaleList), 1);
scaleListMat = repmat(scaleList, size(thetaAlphaList,1),1);
thetaAlphaScaleList = [thetaAlphaScaleList scaleListMat(:)];
if new
    maskStep = 2; 1;2; 3; 4; 2; 4; 2; 1; 2;
else
    maskStep =1;
end
basePose = eye(4);
pad_up = 0;
for i = 1 : 1: length(imgRng)
%     imgCur = imgaussfilt(imresize(imread(fullfile(inputDir, dirInfo(imgRng(i)).name)),2),1.2);
    imgCur = imresize(imread(fullfile(inputDir, dirInfo(imgRng(i)).name)),0.25);
    
    if 1 %~exist('mask', 'var')
        mask1 = zeros(size(imgCur(:,:,1)));
        mask1(1:maskStep:end,:) = 1;
        mask2 = zeros(size(imgCur(:,:,1)));
        mask2(:,1:maskStep:end) = 1;
        mask = (mask1 + mask2) == 1; 2;1; 2; 1; 2;   % 1
        if sum(mask(:)) == 0
            mask = ones(size(mask));
        end
    end
    
    if 0
        [pixCur, edgeCur] = detectEdge(imgCur,imgCur(:,:,1),[0.1 0.13]);
    else
        [pixCur, edgeCur] = detectEdge2(imgCur,imgCur(:,:,1),[0.1 0.13], mask);
        if 0 %size(pixCur,1) > minEdgeSize
            maskStep = 2;
            mask1 = zeros(size(imgCur(:,:,1)));
            mask1(1:maskStep:end,:) = 1;
            mask2 = zeros(size(imgCur(:,:,1)));
            mask2(:,1:maskStep:end) = 1;
            mask = (mask1 + mask2) == 1; 2;   % 1
            if sum(mask(:)) == 0
                mask = ones(size(mask));
            end
            maskStep = 1;
            [pixCur, edgeCur] = detectEdge2(imgCur,imgCur(:,:,1),[0.1 0.13], mask);
        end
    end
    pixCur = [pixCur ones(size(pixCur,1),1)];
%     pixCur = [pixCur 0.1.*(rand(size(pixCur,1),1)-0.5)];
    if i > 1
        if ~exist('ttt','var')
            
            initT = eye(4);
        else
            initT = ttt;
            initT(1:3,4) = transRatio.*ttt(1:3,4);
        end
        if 0
            initT(1:3,1:3) = eye(3);
            initT(2,4) = 0;
        end
        pixPrv0 = (initT(1:3,1:3)*pixPrv' + repmat(initT(1:3,4),1,size(pixPrv,1)))';
         [indexes] = nn_mutual(pixPrv0(:,1:2), pixCur(:,1:2), mergeDist, inlierRati);
         [~,err] = NormalizeVector(pixPrv0(indexes(:,1),:) - pixCur(indexes(:,2),:));
         if 1
             figure(5),clf;showMatchedFeatures(zeros(size(imgCur(:,:,1))), zeros(size(imgCur(:,:,1))),pixPrv0(indexes(:,1),1:2),pixCur(indexes(:,2),1:2));   %(pixPrv0(indexes(:,1),1), pixPrv0(indexes(:,1),2), 'o');hold on;plot(pixCur(indexes(:,2),1),pixCur(indexes(:,2),2),'.');axis equal;drawnow;
%              figure(6),clf;showMatchedFeatures(zeros(size(imgCur(:,:,1))), zeros(size(imgCur(:,:,1))),pixCur00__(:,1:2),pixCur00_(:,1:2));
         end
        
        [pixPrv00, pixCur00, ovlpRatio_] = FindMutualPart(pixPrv0 ,pixCur);
        
        
        
        
        if ovwt
            pixPrv00 = pixPrv0(indexes(:,1),:);
            pixCur00 = pixCur(indexes(:,2),:);
        else
            pixPrv00_ = pixPrv0(indexes(:,1),:);
            pixCur00_ = pixCur(indexes(:,2),:);
        end
        
        
        if new
            [T ,result]= alignEdge(pixPrv0, pixCur, thetaAlphaScaleList, thetaAlphaList, scaleList, mergeDist, thetaMat, alphaMat, imgRng(i), ovlpRatio, result, inlierRati, doOpt2);
            cnt=1;
            while(result(end,4) > maxErr)
                if inlierRati-0.1*cnt < minOvlp
                    break;
                end
                [T ,result]= alignEdge(pixPrv0, pixCur, thetaAlphaScaleList, thetaAlphaList, scaleList, mergeDist, thetaMat, alphaMat, imgRng(i), ovlpRatio, result, inlierRati-0.1*cnt, doOpt2);
                cnt=cnt+1;
            end
        end
        
        if 1 %~ovwt
            try
                [tform,~,rmse] = pcregrigid(pointCloud(pixPrv00), pointCloud(pixCur00),'InlierRatio',0.9,'MaxIterations',1000); %,'Extrapolate',true); %,'Metric','pointToPlane');
            catch
                sgkhj = 1;
            end
            if 0
                [Tx,Ty,w] = rel_pose(pixPrv00_(:,1:2)',pixCur00_(:,1:2)');
                [ data_g, data_p, err, data_pp, R ] = icp_process( pixCur00_, pixPrv00_  );
                figure(6),clf;plot(data_g(:,1), data_g(:,2),'o');hold on;plot(data_p(:,1), data_p(:,2),'.');axis equal;drawnow;
            else
                [rt, pixCur00__] = AlignCoord(pixCur00_,pixPrv00_);
%                 tform.T = [rt' [0 0 0 1]'];
            end
%              figure(6),clf;showMatchedFeatures(zeros(size(imgCur(:,:,1))), zeros(size(imgCur(:,:,1))),pixCur00__(:,1:2),pixCur00_(:,1:2));
        else
            [tform,~,rmse] = pcregrigid(pointCloud(pixPrv00), pointCloud(pixCur00),'InlierRatio',0.9,'MaxIterations',1000); %,'Extrapolate',true); %,'Metric','pointToPlane');
            [Tx,Ty,w] = rel_pose(pixPrv00(:,1:2)',pixCur00(:,1:2)');
        end
        
        if new
            tform.T = [T' ];
        end
        ttt = tform.T'*initT;
        pixCurNew = (ttt(1:3,1:3)*pixPrv' + repmat(ttt(1:3,4),1,size(pixPrv,1)))';
        
        basePose = basePose*inv(ttt);
        basePose_ = basePose;
        if 0
            basePose_(2,4) = basePose_(2,4) + pad_up;
        else
            basePose_inv = inv(basePose_);
            basePose_inv(2,4) = basePose_inv(2,4) - pad_up;
            basePose_ = inv(basePose_inv);
        end
        [img1, pad_up0] = stitch(double(imgCur), double(img0), basePose_);
        pad_up = pad_up + pad_up0;
        figure(1),subplot(1,2,1);cla;plot(pixCurNew(:,1), pixCurNew(:,2),'.r');hold on;plot(pixCur00(:,1), pixCur00(:,2),'.b');
        axis equal;title(num2str(imgRng(i)));drawnow;
        subplot(1,2,2);imshow(img1);
        img0 = img1;
        basePose_inv = inv(basePose);
         pixPrvNew2 = (basePose(1:3,1:3)*pixCur00' + repmat(basePose(1:3,4),1,size(pixCur00,1)))';
         
         figure(3),clf;plot(pix0(:,1), pix0(:,2),'.r');hold on;plot(pixPrvNew2(:,1), pixPrvNew2(:,2),'.b');
          axis equal;title(num2str(imgRng(i)));drawnow;
        pix0 = [pix0; pixPrvNew2];
        imgPrv = imgCur;
        pixPrv = pixCur;
        edgePrv = edgeCur;
        sgu = 1;
    else
        
        
        
        img0 = imgCur;
        pix0 = pixCur; 
        imgPrv = imgCur;
        pixPrv = pixCur;
        edgePrv = edgeCur;
    end
    
end

theta = atand(pad_up/size(img1,2));
img11 = imrotate(img1, -theta);
figure,subplot(1,2,1);imshow(img11);subplot(1,2,2);plot(1./(result(:,5)./result(:,6)));
end

function [pt00, pt11, ovlpRatio] = FindMutualPart(pt0 ,pt1)
global mergeDist
% mergeDist = 1;
scal = 1* mergeDist;
nearestIdx1 = knnsearch(pt0, pt1, 'NSMethod', 'kdtree');
    try
        samePtFlag1 = VecNorm(pt1 - pt0(nearestIdx1, :), 2) <= 1.*scal; 0.2.*scal; 0.4; % 0.2; %obj.featPtManager.configParam.radius_thresh_for_combine; %
    catch
        sghkbj = 1;
    end
    nearestIdx2 = knnsearch(pt1, pt0, 'NSMethod', 'kdtree');
    try
        samePtFlag2 = VecNorm(pt0 - pt1(nearestIdx2, :), 2) <= 1.*scal; 0.2.*scal; 0.4;  % 0.2;% obj.featPtManager.configParam.radius_thresh_for_combine; %
    catch
        dfth = 1;
    end
   pt00 = pt0(samePtFlag2,:);
   pt11 = pt1(samePtFlag1,:);
   if 1
%        figure(2),clf;pcshow(pt0(samePtFlag2,:),repmat([1 0 0],sum(samePtFlag2),1),'MarkerSize',100);hold on;pcshow(pt1(samePtFlag1,:),repmat([0 0 1],sum(samePtFlag1),1),'MarkerSize',100);drawnow;
       figure(2),clf;plot(pt0(samePtFlag2,1),pt0(samePtFlag2,2),'.r');hold on;plot(pt1(samePtFlag1,1),pt1(samePtFlag1,2),'.b');drawnow;
   end
   ovlpRatio = max(sum(samePtFlag2)/length(samePtFlag2), sum(samePtFlag1)/length(samePtFlag1));
end
function n = VecNorm(v, dim)

% VecNorm computes the 2-norms of a list of vectors. 
% v: M by N matrix, a list of vectors. In which dimension it's a vector is
% specified by the second parameter dim.
% dim: integer scalar specifying in which dimension the data in v are seen
% as vectors. Default is 1.
% Returned n is a vector of size M (if dim == 2) or N (if dim == 1) which
% are the 2-norms of given vectors.
%
% By Ji Zhou

if ~exist('dim', 'var')
    dim = 1;
end

n = sqrt(dot(v,v,dim));

end

function [indexes] = nn_mutual(scanA, scanB, thres, inlierRatio)
warning off all;
%if nargin<3, thres = 0.5; end

NA = size(scanA, 1);% NA size [N, 2]; N=361;
NB = size(scanB, 1);% NB size [N, 2];
XA2 = sum(scanA.^2,2);% XA2 size [N, 1]; x^2+y^2;
XB2 = sum(scanB.^2,2);% XB2 size [N, 1];
XB2T = XB2';% 1 by N;
distance = XA2(:,ones(1,NB))+XB2T(ones(1,NA),:)-2*scanA(:,1)*scanB(:,1)'-2*scanA(:,2)*scanB(:,2)';

[mindist,indexes] = min(distance,[],2);
% mindist is N*1, the minimum distance between the nth point in A and all of the points in B.
% indexes is N*1, the index of the points in B which is close to the nth point in A

j=1;
if 0
    for i=1:length(indexes) % i is the index of point in A
        if 1 % length(find(indexes==indexes(i)))==1 % if only the point in B which is close to nth point in A, is only close to nth point in A
            indexes2(j,:)=[i, indexes(i)]; % [idx of point in A, idx of point in B which is close to this point in A]
            j=j+1;
            
        else
            
            sdjgfkh = 1;
            
        end
    end
elseif 1
    indexes2 = [[1:length(indexes)]' indexes];
else
    upperBound = max(1, round(inlierRatio(1)*size(scanA, 1)));    

    keepInlierA = false(size(scanA, 1), 1); 
    [~, idx] = sort(dists);
    keepInlierA(idx(1:upperBound)) = true;
    inlierIndicesA = find(keepInlierA);
    inlierIndicesB = indices(keepInlierA);
    
end
%indexes = [[1:NA]' indexes'];
ind = (mindist(indexes2(:,1))<thres^2); % check whether the minimum distance is less than the gate
if 1
    upperBound = max(1, round(inlierRatio(1)*min([size(scanA, 1) size(scanB, 1)])));
    
    keepInlierA = false(size(scanA, 1), 1);
    [~, idx] = sort(mindist);
    keepInlierA(idx(1:upperBound)) = true;
    inlierIndicesA = find(keepInlierA);
    inlierIndicesB = indexes2(keepInlierA,2);
    indexes = [inlierIndicesA inlierIndicesB];
else
    
    % the result is false or true. NOTE: not 0 or 1!!!
    indexes = indexes2(ind,:);
end
end
function [normVec, sqrtsumvec2] = NormalizeVector(vec)
if 1
    vec2 = vec.^2;
    sumvec2 = sum(vec2,2);
    sqrtsumvec2 = sqrt(sumvec2);
    if 0
        for j = 1:size(npixels,1)
            ndescrppp(j,:) = npixels(j,:)/norm(npixels(j,:));
        end
    else
        for k = 1 : size(vec,2)
            normVec(:,k) = vec(:,k)./sqrtsumvec2;
        end
    end
else
    [a] = normc(vec')';
    b=a;
    a(a == 0) = 1;
    scal = (vec./a);
    scal
    id = ~isnan(scal);
    idd = sum(id,1) == size(scal,1);   % idd = find(sum(id) == size(scal,1));
    scall = mean(scal(:,idd),2);
    
    sqrtsumvec2 = scall;
    normVec = b;
end
end
%function [Tx,Ty,w] = rel_pose(scanA,scanB)
%Compute the relative position of two scans based on the least-square
%solution.
%
%INPUTS:
% scanA = N x 2 matrix where N is the number of points
% scanB = N x 2 matrix where N is the number of points
% scanA and scanB need to have the same size
%
%OUTPUTS:
% w relative rotation
% Tx relative translation in X
% Ty relative translation in y
%
% Implemented by Fabio Tozeto Ramos 17/08/05
function [Tx,Ty,w] = rel_pose(scanA,scanB)
N = length(scanA);
M = length(scanB);

%x_dashA = mean(scanA(:,1));
%x_dashB = mean(scanB(:,1));
%y_dashA = mean(scanA(:,2));
%y_dashB = mean(scanB(:,2));
%SxAxB = sum((scanA(:,1)-x_dashA*ones(N,1)).*(scanB(:,1)-x_dashB*ones(N,1)),1);
%SyAyB = sum((scanA(:,2)-y_dashA*ones(N,1)).*(scanB(:,2)-y_dashB*ones(N,1)),1);
%SxAyB = sum((scanA(:,1)-x_dashA*ones(N,1)).*(scanB(:,2)-y_dashB*ones(N,1)),1);
%SyAxB = sum((scanA(:,2)-y_dashA*ones(N,1)).*(scanB(:,1)-x_dashB*ones(N,1)),1);
xA = sum(scanA(:,1)); % sum of the X coords in scan A
xB = sum(scanB(:,1));
yA = sum(scanA(:,2));
yB = sum(scanB(:,2));

SxAxB = sum(scanA(:,1).*scanB(:,1)) - xA*xB/N; % sum(XAn*XBm)- sum(XA)*sum(XB)/N
SyAyB = sum(scanA(:,2).*scanB(:,2)) - yA*yB/N;
SxAyB = sum(scanA(:,1).*scanB(:,2)) - xA*yB/N;
SyAxB = sum(scanA(:,2).*scanB(:,1)) - yA*xB/N;

x_dashA = xA/N;
x_dashB = xB/N;
y_dashA = yA/N;
y_dashB = yB/N;


w = atan2((SxAyB-SyAxB),(SxAxB+SyAyB)); 
Tx = x_dashB - (x_dashA*cos(w) - y_dashA*sin(w));
Ty = y_dashB - (x_dashA*sin(w) + y_dashA*cos(w));
end

function [ data_g, data_p, err, data_pp, R ] = icp_process( data_g, data_p )
% 对两个点集data_g和data_p应用ICP算法
% data_g\data_P:样本点集
% 返回旋转后的两个点集，以及误差均值以及data_p对应data_g的对应点集data_pp

    [k1, n] = size(data_g);
    [k2, m] = size(data_p);
    
    data_p1 = zeros(k2, 3);     % 中间点集
    data_pp = zeros(k1, 3);     % 对应点集
    distance = zeros(k1, 1);    % 点集之间各点的距离
    error = zeros(k1, 1);       % 对应点之间的误差
    
    % 对两个点集作去中心化
    data_g = normal_gravity(data_g);
    data_p = normal_gravity(data_p);
    
    if 1
        % 遍历整个data_g点集，寻找每个点对应data_p点集中距离最小的点，作为对应点
        for i = 1:k1
            data_p1(:, 1) = data_p(:, 1) - data_g(i, 1);    % 两个点集中的点x坐标之差
            data_p1(:, 2) = data_p(:, 2) - data_g(i, 2);    % 两个点集中的点y坐标之差
            data_p1(:, 3) = data_p(:, 3) - data_g(i, 3);    % 两个点集中的点z坐标之差
            distance = data_p1(:, 1).^2 + data_p1(:, 2).^2 + data_p1(:, 3).^2;  % 欧氏距离
            [min_dis, min_index] = min(distance);   % 找到距离最小的那个点
            data_pp(i, :) = data_p(min_index, :);   % 将那个点保存为对应点
            error(i) = min_dis;     % 保存距离差值
        end
    else
        data_pp = data_p;
    end
    % 求出协方差矩阵
    V = (data_g' * data_pp) ./ k1;
    
    % 构建正定矩阵Q（这部分不是很理解，直接套公式了）
    matrix_Q = [V(1,1)+V(2,2)+V(3,3),V(2,3)-V(3,2),V(3,1)-V(1,3),V(1,2)-V(2,1);  
                V(2,3)-V(3,2),V(1,1)-V(2,2)-V(3,3),V(1,2)+V(2,1),V(1,3)+V(3,1);  
                V(3,1)-V(1,3),V(1,2)+V(2,1),V(2,2)-V(1,1)-V(3,3),V(2,3)+V(3,2);  
                V(1,2)-V(2,1),V(1,3)+V(3,1),V(2,3)+V(3,2),V(3,3)-V(1,1)-V(2,2)];
    
    [V2, D2] = eig(matrix_Q);       % 对矩阵Q作特征值分解
    lambdas = [D2(1, 1), D2(2, 2), D2(3, 3), D2(4, 4)]; % 取出特征值
    [lambda, ind] = max(lambdas);   % 求出最大的那个特征值
    Q = V2(:, ind); % 取出那个最大的特征值所对应的特征向量
    
    % 构建旋转矩阵（四元数）
    R=[Q(1,1)^2+Q(2,1)^2-Q(3,1)^2-Q(4,1)^2,     2*(Q(2,1)*Q(3,1)-Q(1,1)*Q(4,1)),        2*(Q(2,1)*Q(4,1)+Q(1,1)*Q(3,1));  
       2*(Q(2,1)*Q(3,1)+Q(1,1)*Q(4,1)),         Q(1,1)^2-Q(2,1)^2+Q(3,1)^2-Q(4,1)^2,    2*(Q(3,1)*Q(4,1)-Q(1,1)*Q(2,1));  
       2*(Q(2,1)*Q(4,1)-Q(1,1)*Q(3,1)),         2*(Q(3,1)*Q(4,1)+Q(1,1)*Q(2,1)),        Q(1,1)^2-Q(2,1)^2-Q(3,1)^2-Q(4,1)^2;  
    ];
    
    % 对data_p点集所有的点都做R的旋转变化，然后再作中心平移
    data_p = data_p * R;
    data_pp = data_pp * R;
    data_p = normal_gravity(data_p);
    data_pp = normal_gravity(data_pp);
    err = mean(error);
    
end
function ret = normal_gravity( data )
% 将点集data去中心化
% 其中行为不同样本，列为一个样本的特征

    [m, n] = size(data);
    data_mean = mean(data);
    ret = data - ones(m, 1) * data_mean;

end

function [rt,mapped] = AlignCoord(matchPtBody,matchPtVR)
bodyCtrMass = mean(matchPtBody);
vrCtrMass = mean(matchPtVR);
bodyPt = matchPtBody - repmat(bodyCtrMass,size(matchPtBody,1),1);
vrPt = matchPtVR - repmat(vrCtrMass,size(matchPtVR,1),1);
[U,S,V] = svd(vrPt'*bodyPt);
H = U*S*V';
X = V*U';
% %     [U,S,V] = svd(bodyPt'*vrPt);
% %     H = U*S*V';
% %     X = U'*V;
if abs(det(X) - (-1)) < 0.0001
    %     if X(2,2) < 0
    V(:,3) = -V(:,3);
    R = V*U';
else
    R = X;
    
end
% %     R(3,3) = -R(3,3);
t = bodyCtrMass' - R*vrCtrMass';
rt = [R t];


mapped = (R*matchPtVR' + repmat(t,1,size(matchPtVR,1)))';
if 0
    figure,plot3(matchPtBody(:,1),matchPtBody(:,2),matchPtBody(:,3),'or');hold on;plot3(mapped(:,1),mapped(:,2),mapped(:,3),'.g');axis equal;
    [~,err] = NormalizeVector(matchPtBody - mapped);
end
end

function [posOpt , result]= alignEdge(pixPrv, pixCur, thetaAlphaScaleList, thetaAlphaList, scaleList0, mergeDist0, thetaMat, alphaMat, imgInd, ovlpRatio, result, inlierRatio, doOpt2)
% global result
mergeDist = mergeDist0;
scaleList = scaleList0;
range = [min(scaleList) max(scaleList)];


if 0
    
    for i = 1 : size(thetaAlphaScaleList, 1)
        theta = thetaAlphaScaleList(i,1);
        alpha = thetaAlphaScaleList(i,2);
        scale = thetaAlphaScaleList(i,3);
        
        r = rotz(theta);
        t = [-cosd(alpha); sind(alpha); 0];
        
        pixCurMapped = (r*pixPrv' + repmat(scale.*t, 1, size(pixPrv, 1)))';
        pixPrvX(:,i) = pixCurMapped(:,1);
        pixPrvY(:,i) = pixCurMapped(:,2);
        pixPrvZ(:,i) = pixCurMapped(:,3);
        
        
    end
    pixCurX = repmat(pixCur(:,1), 1, size(thetaAlphaScaleList, 1));
    pixCurY = repmat(pixCur(:,2), 1, size(thetaAlphaScaleList, 1));
    pixCurZ = repmat(pixCur(:,3), 1, size(thetaAlphaScaleList, 1));
    
    prvList = [pixPrvX(:) pixPrvY(:)];
    curList = [pixCurX(:) pixCurY(:)];
    if 1
        D_temp = pdist2(prvList, curList,'euclidean');
    else
        [indexes11] = nn_mutual2(prvList(:,1:2), curList(:,1:2), mergeDist);
    end
end

tStart = tic;
for i = 1 : size(thetaAlphaList, 1)
   theta = thetaAlphaList(i,1); 
   alpha = thetaAlphaList(i,2); 
   r = rotz(theta);
   t = [-cosd(alpha); sind(alpha); 0];
   
   if 0
       [optimX(i,1), ~, exytFlag] = fminbnd(@ObjectiveFunction, range(1), range(end), optimset('TolX', 0.01, 'Display', 'off'));
   else
       for j = 1 : length(scaleList)
           
           [er(j,1)] = ObjectiveFunction(scaleList(j));
           
       end
       [~, minInd] = min(er);
       optimX(i,1) = scaleList(minInd);
   end
   [err(i,1), Err, pixPrv00_, pixCur00_, ovlp(i,1)] = ObjectiveFunction2(optimX(i,1));
   Var(i,1) = var(Err);
end

tElapsed = toc(tStart);

if 1
    ind = find(ovlp >=ovlpRatio);
    while(isempty(ind))
        ovlpRatio = ovlpRatio/2;
        ind = find(ovlp >=ovlpRatio);
    end
elseif 0
    ind = find(ovlp >=ovlpRatio & Var <= mean(Var));
    while(isempty(ind))
        ovlpRatio = ovlpRatio/2;
        ind = find(ovlp >=ovlpRatio & Var <= mean(Var));
    end
elseif 0
    ind = find(Var <= mean(Var));
else
    ind = [1:length(err)]';
end
[~,minId_] = min(err(ind));

minId = ind(minId_);
ovlp_ = ovlp(minId);
thetaOpt = thetaAlphaList(minId, 1); 
alphaOpt = thetaAlphaList(minId, 2); 
scaleOpt = optimX(minId);
errOpt = err(minId);
if 1
    if 0
        result = [result; [thetaOpt alphaOpt scaleOpt ovlp_]];
    else
        result = [result; [thetaOpt alphaOpt scaleOpt errOpt tElapsed size(pixCur,1)]];
    end
end
T = [rotz(thetaOpt) scaleOpt.*[-cosd(alphaOpt); sind(alphaOpt); 0]; 0 0 0 1];

 r = rotz(thetaOpt);
 t = [-cosd(alphaOpt); sind(alphaOpt); 0];
 

[errOpt, ErrOpt, pixPrv00_opt, pixCur00_opt] = ObjectiveFunction2(scaleOpt);

 pixPrv_1 = (r*pixPrv' + repmat(scaleOpt.*t, 1, size(pixPrv, 1)))';

 
 errMat = reshape(err, size(thetaMat));
 if 0
     if 1
         figure(20),subplot(2,2,1);cla;surf(thetaMat, alphaMat, errMat);hold on; plot3(thetaOpt, alphaOpt, errOpt, '*r');title(sprintf('imgInd: %d\ntheta: %0.2f, alpha: %0.2f, scale: %0.2f\novlpRatio: %0.2f,   err: %0.5f',imgInd, thetaOpt, alphaOpt, scaleOpt, ovlp_, errOpt));subplot(2,2,2); showMatchedFeatures(zeros(240,320), zeros(240, 320),pixPrv00_opt(:,1:2),pixCur00_opt(:,1:2));subplot(2,2,4);cla; plot(pixPrv_1(:,1),pixPrv_1(:,2),'o'); hold on;plot(pixCur(:,1),pixCur(:,2), '.');axis equal;subplot(2,2,3),cla;plot(result);
     else
         figure(20),subplot(1,2,1);cla;surf(thetaMat, alphaMat, errMat);hold on; plot3(thetaOpt, alphaOpt, errOpt, '*r');title(sprintf('imgInd: %d\ntheta: %0.2f, alpha: %0.2f, scale: %0.2f\novlpRatio: %0.2f',imgInd, thetaOpt, alphaOpt, scaleOpt, ovlp_));subplot(1,2,2);cla; plot(pixPrv_1(:,1),pixPrv_1(:,2),'o'); hold on;plot(pixCur(:,1),pixCur(:,2), '.');axis equal;
     end
 end
 posOpt = icpIter(pixPrv, pixCur,  T, mergeDist, inlierRatio);
 
 
 
  r = posOpt(1:3,1:3);
 t = posOpt(1:3,4)./norm(posOpt(1:3,4));
 

[errOpt11, ErrOpt11, pixPrv00_opt11, pixCur00_opt11] = ObjectiveFunction2(norm(posOpt(1:3,4)));

 pixPrv_11 = (posOpt(1:3,1:3)*pixPrv' + repmat(posOpt(1:3,4), 1, size(pixPrv, 1)))';
 
  figure(20),subplot(2,3,1);cla;surf(thetaMat, alphaMat, errMat);hold on; plot3(thetaOpt, alphaOpt, errOpt, '*r');title(sprintf('imgInd: %d\ntheta: %0.2f, alpha: %0.2f, scale: %0.2f\novlpRatio: %0.2f,   err: %0.5f',imgInd, thetaOpt, alphaOpt, scaleOpt, ovlp_, errOpt));
  subplot(2,3,2); showMatchedFeatures(zeros(240,320), zeros(240, 320),pixPrv00_opt(:,1:2),pixCur00_opt(:,1:2));
  subplot(2,3,3); showMatchedFeatures(zeros(240,320), zeros(240, 320),pixPrv00_opt11(:,1:2),pixCur00_opt11(:,1:2));
  subplot(2,3,5);cla; plot(pixPrv_1(:,1),pixPrv_1(:,2),'o'); hold on;plot(pixCur(:,1),pixCur(:,2), '.');axis equal;
  subplot(2,3,6);cla; plot(pixPrv_11(:,1),pixPrv_11(:,2),'o'); hold on;plot(pixCur(:,1),pixCur(:,2), '.');axis equal;
  subplot(2,3,4),cla;plot(result);
  
  
  
  if ~doOpt2
      posOpt = T;
  else
      result(end,4) = errOpt11;
  end
 hadhf=1;
 
 sdgkhj = 1;
 
    function err = ObjectiveFunction(X)
        
        pixPrv_ = (r*pixPrv' + repmat(X.*t, 1, size(pixPrv, 1)))';
        [indexes] = nn_mutual(pixPrv_(:,1:2), pixCur(:,1:2), mergeDist, inlierRatio);
        
        pixPrv00_ = pixPrv_(indexes(:,1),:);
        pixCur00_ = pixCur(indexes(:,2),:);
        Error = sqrt((pixPrv00_(:,1) - pixCur00_(:,1)).^2 + (pixPrv00_(:,2) - pixCur00_(:,2)).^2);
        err = mean(Error);
    end

    function [err1, Err1, pixPrv00_, pixCur00_, overlap] = ObjectiveFunction2(X)
        
        pixPrv_ = (r*pixPrv' + repmat(X.*t, 1, size(pixPrv, 1)))';
        [indexes] = nn_mutual(pixPrv_(:,1:2), pixCur(:,1:2), mergeDist, inlierRatio);
        %         [indexes] = nn_mutual(pixCur(:,1:2), pixPrv_(:,1:2), mergeDist);
        if 0
            overlap = length(unique(indexes(:,2)))./size(indexes,1);
        else
            overlap = size(indexes,1)/min([size(pixPrv_,1) size(pixCur,1)]);
        end
        pixPrv00_ = pixPrv_(indexes(:,1),:);
        pixCur00_ = pixCur(indexes(:,2),:);
        Err1 = sqrt((pixPrv00_(:,1) - pixCur00_(:,1)).^2 + (pixPrv00_(:,2) - pixCur00_(:,2)).^2);
        err1 = mean(Err1);
        if 0
            figure(15),clf;showMatchedFeatures(zeros(240,320), zeros(240, 320),pixPrv00_(:,1:2),pixCur00_(:,1:2));
        end
    end

end
function [indexes] = nn_mutual2(scanA, scanB, thres)
warning off all;
%if nargin<3, thres = 0.5; end

NA = size(scanA, 1);% NA size [N, 2]; N=361;
NB = size(scanB, 1);% NB size [N, 2];
XA2 = sum(scanA.^2,2);% XA2 size [N, 1]; x^2+y^2;
XB2 = sum(scanB.^2,2);% XB2 size [N, 1];
XB2T = XB2';% 1 by N;
distance = XA2(:,ones(1,NB))+XB2T(ones(1,NA),:)-2*scanA(:,1)*scanB(:,1)'-2*scanA(:,2)*scanB(:,2)';

[mindist,indexes] = min(distance,[],2);
% mindist is N*1, the minimum distance between the nth point in A and all of the points in B.
% indexes is N*1, the index of the points in B which is close to the nth point in A

j=1;
if 0
    for i=1:length(indexes) % i is the index of point in A
        if 1 % length(find(indexes==indexes(i)))==1 % if only the point in B which is close to nth point in A, is only close to nth point in A
            indexes2(j,:)=[i, indexes(i)]; % [idx of point in A, idx of point in B which is close to this point in A]
            j=j+1;
            
        else
            
            sdjgfkh = 1;
            
        end
    end
else
    indexes2 = [[1:length(indexes)]' indexes];
end
%indexes = [[1:NA]' indexes'];
ind = (mindist(indexes2(:,1))<thres^2); % check whether the minimum distance is less than the gate
% the result is false or true. NOTE: not 0 or 1!!!
indexes = indexes2(ind,:);
end
function posOpt = icpIter(pointsA, pointsB,  relpos_ba, thres_trans, inlierRatio)

global stopThr

% pointsA = get_laser_points(rawscanA, 5.8); % [x; y]; (2*361)
% pointsB = get_laser_points(rawscanB, 5.8);
rotVec = (rodrigues(relpos_ba(1:3,1:3)));
pos = [relpos_ba(1,4) relpos_ba(2,4) rotVec(3)];


nit = 200;
relpos_ba = pos;

prevpos=relpos_ba;
for i=1:nit

    pointsA_rel_B = (rotz(rad2deg(pos(3))) * pointsA' + repmat([pos(1:2) 0]', 1, size(pointsA,1)))';
    indexes=nn_mutual(pointsA_rel_B,pointsB,thres_trans, inlierRatio);
    
    
    
    
    if 0
        figure,subplot(1,2,1);showMatchedFeatures(zeros(240,320),zeros(240,320),pointsA_rel_B(indexes(:,1),1:2),pointsB(indexes(:,2),1:2));        subplot(1,2,2);plot(pointsA_rel_B(indexes(:,1),1),pointsA_rel_B(indexes(:,1),2),'o');hold on;plot(pointsB(indexes(:,2),1),pointsB(indexes(:,2),2),'.');axis equal
    end
    if 0
        [post(1),post(2),post(3)]=rel_pose(pointsA_rel_B(indexes(:,1),1:2),pointsB(indexes(:,2),1:2));
    else
        [post] = rigidTransform3D(pointsA_rel_B(indexes(:,1),:),pointsB(indexes(:,2),:));
    end
    pos=pos+post;
    posdist=sqrt((pos(1)-prevpos(1))^2+(pos(2)-prevpos(2))^2);
   
    if posdist<stopThr
        break;
    else
        prevpos=pos;
    end

end

posOpt = [rotz(rad2deg(pos(3))) [pos(1:2) 0]'; 0 0 0 1];

end
function [pos] = rigidTransform3D(p, q)

n = cast(size(p, 1), 'like', p);
m = cast(size(q, 1), 'like', q);

% Find data centroid and deviations from centroid
pmean = sum(p,1)/n;
p2 = bsxfun(@minus, p, pmean);

qmean = sum(q,1)/m;
q2 = bsxfun(@minus, q, qmean);

% Covariance matrix
C = p2'*q2;

[U,~,V] = svd(C);

% Handle the reflection case
R = V*diag([1 1 sign(det(U*V'))])*U';

% Compute the translation
t = qmean' - R*pmean';
relpos_ba = [R t;0 0 0 1];
rotVec = rodrigues(relpos_ba(1:3,1:3));
pos = [relpos_ba(1,4) relpos_ba(2,4) rotVec(3)];

end
function	[out,dout]=rodrigues(in)

% RODRIGUES	Transform rotation matrix into rotation vector and viceversa.
%		
%		Sintax:  [OUT]=RODRIGUES(IN)
% 		If IN is a 3x3 rotation matrix then OUT is the
%		corresponding 3x1 rotation vector
% 		if IN is a rotation 3-vector then OUT is the 
%		corresponding 3x3 rotation matrix
%

%%
%%		Copyright (c) March 1993 -- Pietro Perona
%%		California Institute of Technology
%%

%% ALL CHECKED BY JEAN-YVES BOUGUET, October 1995.
%% FOR ALL JACOBIAN MATRICES !!! LOOK AT THE TEST AT THE END !!

%% BUG when norm(om)=pi fixed -- April 6th, 1997;
%% Jean-Yves Bouguet

%% Add projection of the 3x3 matrix onto the set of special ortogonal matrices SO(3) by SVD -- February 7th, 2003;
%% Jean-Yves Bouguet

% BUG FOR THE CASE norm(om)=pi fixed by Mike Burl on Feb 27, 2007


[m,n] = size(in);
%bigeps = 10e+4*eps;
bigeps = 10e+20*eps;

if ((m==1) & (n==3)) | ((m==3) & (n==1)) %% it is a rotation vector
    theta = norm(in);
    if theta < eps
        R = eye(3);

        %if nargout > 1,

        dRdin = [0 0 0;
            0 0 1;
            0 -1 0;
            0 0 -1;
            0 0 0;
            1 0 0;
            0 1 0;
            -1 0 0;
            0 0 0];

        %end;

    else
        if n==length(in)  in=in'; end; 	%% make it a column vec. if necess.

        %m3 = [in,theta]

        dm3din = [eye(3);in'/theta];

        omega = in/theta;

        %m2 = [omega;theta]

        dm2dm3 = [eye(3)/theta -in/theta^2; zeros(1,3) 1];

        alpha = cos(theta);
        beta = sin(theta);
        gamma = 1-cos(theta);
        omegav=[[0 -omega(3) omega(2)];[omega(3) 0 -omega(1)];[-omega(2) omega(1) 0 ]];
        A = omega*omega';

        %m1 = [alpha;beta;gamma;omegav;A];

        dm1dm2 = zeros(21,4);
        dm1dm2(1,4) = -sin(theta);
        dm1dm2(2,4) = cos(theta);
        dm1dm2(3,4) = sin(theta);
        dm1dm2(4:12,1:3) = [0 0 0 0 0 1 0 -1 0;
            0 0 -1 0 0 0 1 0 0;
            0 1 0 -1 0 0 0 0 0]';

        w1 = omega(1);
        w2 = omega(2);
        w3 = omega(3);

        dm1dm2(13:21,1) = [2*w1;w2;w3;w2;0;0;w3;0;0];
        dm1dm2(13: 21,2) = [0;w1;0;w1;2*w2;w3;0;w3;0];
        dm1dm2(13:21,3) = [0;0;w1;0;0;w2;w1;w2;2*w3];

        R = eye(3)*alpha + omegav*beta + A*gamma;

        dRdm1 = zeros(9,21);

        dRdm1([1 5 9],1) = ones(3,1);
        dRdm1(:,2) = omegav(:);
        dRdm1(:,4:12) = beta*eye(9);
        dRdm1(:,3) = A(:);
        dRdm1(:,13:21) = gamma*eye(9);

        dRdin = dRdm1 * dm1dm2 * dm2dm3 * dm3din;


    end;
    out = R;
    dout = dRdin;

    %% it is prob. a rot matr.
elseif ((m==n) & (m==3) & (norm(in' * in - eye(3)) < bigeps)...
        & (abs(det(in)-1) < bigeps))
    R = in;

    % project the rotation matrix to SO(3);
    [U,S,V] = svd(R);
    R = U*V';

    tr = (trace(R)-1)/2;
    dtrdR = [1 0 0 0 1 0 0 0 1]/2;
    theta = real(acos(tr));


    if sin(theta) >= 1e-5,     %1e-4, i changed here a little bit, 2017/08/31, for imu intgration

        dthetadtr = -1/sqrt(1-tr^2);

        dthetadR = dthetadtr * dtrdR;
        % var1 = [vth;theta];
        vth = 1/(2*sin(theta));
        dvthdtheta = -vth*cos(theta)/sin(theta);
        dvar1dtheta = [dvthdtheta;1];

        dvar1dR =  dvar1dtheta * dthetadR;


        om1 = [R(3,2)-R(2,3), R(1,3)-R(3,1), R(2,1)-R(1,2)]';

        dom1dR = [0 0 0 0 0 1 0 -1 0;
            0 0 -1 0 0 0 1 0 0;
            0 1 0 -1 0 0 0 0 0];

        % var = [om1;vth;theta];
        dvardR = [dom1dR;dvar1dR];

        % var2 = [om;theta];
        om = vth*om1;
        domdvar = [vth*eye(3) om1 zeros(3,1)];
        dthetadvar = [0 0 0 0 1];
        dvar2dvar = [domdvar;dthetadvar];


        out = om*theta;
        domegadvar2 = [theta*eye(3) om];

        dout = domegadvar2 * dvar2dvar * dvardR;


    else
        if tr > 0; 			% case norm(om)=0;

            out = [0 0 0]';

            dout = [0 0 0 0 0 1/2 0 -1/2 0;
                0 0 -1/2 0 0 0 1/2 0 0;
                0 1/2 0 -1/2 0 0 0 0 0];
        else

            % case norm(om)=pi;
            if(0)

                %% fixed April 6th by Bouguet -- not working in all cases!
                out = theta * (sqrt((diag(R)+1)/2).*[1;2*(R(1,2:3)>=0)'-1]);
                %keyboard;

            else

                % Solution by Mike Burl on Feb 27, 2007
                % This is a better way to determine the signs of the
                % entries of the rotation vector using a hash table on all
                % the combinations of signs of a pairs of products (in the
                % rotation matrix)

                % Define hashvec and Smat
                hashvec = [0; -1; -3; -9; 9; 3; 1; 13; 5; -7; -11];
                Smat = [1,1,1; 1,0,-1; 0,1,-1; 1,-1,0; 1,1,0; 0,1,1; 1,0,1; 1,1,1; 1,1,-1;
                    1,-1,-1; 1,-1,1];

                M = (R+eye(3,3))/2;
                uabs = sqrt(M(1,1));
                vabs = sqrt(M(2,2));
                wabs = sqrt(M(3,3));

                mvec = [M(1,2), M(2,3), M(1,3)];
                syn  = ((mvec > 1e-4) - (mvec < -1e-4)); % robust sign() function
                hash = syn * [9; 3; 1];
                idx = find(hash == hashvec);
                svec = Smat(idx,:)';

                out = theta * [uabs; vabs; wabs] .* svec;

            end;

            if nargout > 1,
                fprintf(1,'WARNING!!!! Jacobian domdR undefined!!!\n');
                dout = NaN*ones(3,9);
            end;
        end;
    end;

else
    error('Neither a rotation matrix nor a rotation vector were provided');
end;

return;

%% test of the Jacobians:

%%%% TEST OF dRdom:
om = randn(3,1);
dom = randn(3,1)/1000000;

[R1,dR1] = rodrigues(om);
R2 = rodrigues(om+dom);

R2a = R1 + reshape(dR1 * dom,3,3);

gain = norm(R2 - R1)/norm(R2 - R2a)

%%% TEST OF dOmdR:
om = randn(3,1);
R = rodrigues(om);
dom = randn(3,1)/10000;
dR = rodrigues(om+dom) - R;

[omc,domdR] = rodrigues(R);
[om2] = rodrigues(R+dR);

om_app = omc + domdR*dR(:);

gain = norm(om2 - omc)/norm(om2 - om_app)


%%% OTHER BUG: (FIXED NOW!!!)

omu = randn(3,1);   
omu = omu/norm(omu)
om = pi*omu;        
[R,dR]= rodrigues(om);
[om2] = rodrigues(R);
[om om2]

%%% NORMAL OPERATION

om = randn(3,1);         
[R,dR]= rodrigues(om);
[om2] = rodrigues(R);
[om om2]

return

% Test: norm(om) = pi

u = randn(3,1);
u = u / sqrt(sum(u.^2));
om = pi*u;
R = rodrigues(om);

R2 = rodrigues(rodrigues(R));

norm(R - R2)
end