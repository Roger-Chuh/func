function [pose_rel, score, warped_image, valid_mask, weightMap, warped_z, validMask] = estimateVisualOdometry(img_curr, img_prev, dep_curr, dep_prev, K, num_levels, pose_init, isVerbose,initial_sigma,default_dof,para)
% estimate visual odometry from two adjacent frames
%
% INPUT:
%   img_curr: current RGB frame
%   img_prev: previous RGB frame
%   dep_prev: previous depth frame
%   K: intrinsic parameters of the camera
%   num_levels: pyramid levels for the estimation
%   isVerbose: whether to show intermediate results
%
% OUTPUT:
%   pose_rel: relative pose
%   score: fitting score

% construct image pyramids
img_curr_pyr = constructPyramid(img_curr, num_levels);
img_prev_pyr = constructPyramid(img_prev, num_levels);

% construct depth pyramid for the previous frame
dep_prev_pyr = constructPyramid(dep_prev, num_levels);
dep_curr_pyr = constructPyramid(dep_curr, num_levels);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimate the relative pose from coarse to fine scales

% initialize the relative pose and its increment
pose_rel = pose_init; % eye(4);
increment = zeros(6, 1);

% convert the vectorize motion increment to update the relative pose
increment = twistexp(increment);
pose_rel = increment * pose_rel;

% modify the camera parameters to fit each pyramid
K_pyr = K;
K_pyr(1:2, :) = K_pyr(1:2, :) / (2^(num_levels-1));
% K_pyr(1:2, :) = K_pyr(1:2, :) / ((1/0.8)^(num_levels-1));

[height_end, width_end] = size(dep_prev_pyr{end});
if 0
    K_pyr(1, 3) = (1 + width_end) / 2;
    K_pyr(2, 3) = (1 + height_end) / 2;
end

% weightMap1 = zeros(size(img_curr(:,:,1)));

if 1
    minErr = 1e-10;
    minStep = 1e-10;
    maxIter = 100; 15; 100; 50; 500; 100; 50;
else
    minErr = 1e-7;
    minStep = 1e-9;
    maxIter = 500; 100; 50; 100; 50;
end

minIterNum = 1; 3; 10;  3;

% prev_sqError = inf;
for n = num_levels:-1:1
    
    % image size
    if 0
        [height, width] = size(dep_prev_pyr{n});
        
        % get valid point clouds in the previous frame
        [pointclouds, valid_mask] = reprojectDepthImage(dep_prev_pyr{n}, K_pyr);
        
        % warp pointclouds and prune invalid points
        
        warped_pointclouds = warpPointCloud(pointclouds, pose_rel);
        [warped_img_coordinates, valid_points, ind] = projectPointCloud(warped_pointclouds, K_pyr, height, width);
        
        warped_pointclouds = warped_pointclouds(valid_points, :);
        pointclouds = pointclouds(valid_points, :);
        
        % spatial gradient in the current frame
        [Gx_curr, Gy_curr] = imgradientxy(img_curr_pyr{n}, 'CentralDifference');
        Gx_curr = interp2(Gx_curr, warped_img_coordinates(:, 1), warped_img_coordinates(:, 2), 'linear', 0);
        Gy_curr = interp2(Gy_curr, warped_img_coordinates(:, 1), warped_img_coordinates(:, 2), 'linear', 0);
        Gs_curr = cat(2, Gx_curr, Gy_curr);
        
        % temporal visual difference
        Gt_prev = img_prev_pyr{n};
        Gt_prev = Gt_prev(valid_mask);
        Gt_prev = Gt_prev(valid_points);
        Gt_curr = img_curr_pyr{n};
        Gt_curr = interp2(Gt_curr, warped_img_coordinates(:, 1), warped_img_coordinates(:, 2), 'linear', 0);
        Gt = Gt_curr - Gt_prev;
        
        % calculate the warping Jacobian
        warp_jacobian = calculateWarpingJacobian(warped_pointclouds, pointclouds, pose_rel, K_pyr);
        
        % calculate the compositive jacobian
        comp_jacobian = squeeze(sum(bsxfun(@times, Gs_curr, warp_jacobian), 2));
        if 1
            % adjust the importance of each residual elements, when required
            % necessary pruning of bad correspondences
            % compute the weight
            variance = computeResidualVariance(Gt, para);
            weights = sqrt((para + 1) ./ (para + Gt.^2./variance));
            comp_jacobian = bsxfun(@times, comp_jacobian, weights);
            Gt = Gt.*weights;
            
            weightMap = zeros(height, width);
            weightMap(ind) = weights;
        end
        
        % calculate the increment motion
        increment = -(comp_jacobian'*comp_jacobian)\(comp_jacobian'*Gt);
        
        % get the current relative pose
        increment = twistexp(increment);
        pose_rel = increment * pose_rel;
    else
        InParas = [K_pyr(1,1) K_pyr(2,2) K_pyr(1,3) K_pyr(2,3)];
        if 1
            prev_sqError = inf;
            for kk = 1 : maxIter
                [pose_rel2, error, weightMap] = iterPose(pose_rel, K_pyr, img_curr_pyr, img_prev_pyr, dep_curr_pyr, dep_prev_pyr, n,InParas,initial_sigma,default_dof,para);
                if length(prev_sqError) > minIterNum
                    if (error > prev_sqError(end))
%                     if (error > prev_sqError(end) || error < minErr || prev_sqError(end) - error < minStep)
                        %                 pose_rel = pose_rel2;
                        break;
                    end
                     if ( error < minErr || prev_sqError(end) - error < minStep)
                          pose_rel = pose_rel2;
                         break;
                     end
                end
                prev_sqError = [prev_sqError; error];
                pose_rel = pose_rel2;
            end
        elseif 1
            %             InParas = [K_pyr(1,1) K_pyr(2,2) K_pyr(1,3) K_pyr(2,3)];
            [pose_rel2, e1, r1, W1, weightMap ]= EstimateCameraMotion(img_prev_pyr{n},dep_prev_pyr{n},img_curr_pyr{n},dep_curr_pyr{n},pose_rel,...
                InParas, minStep,maxIter,initial_sigma,default_dof,1);
            
            pose_rel = pose_rel2;
        elseif 1
            pose_rel_ = pose_rel;
            pose_rel_(1:3,4) = pose_rel_(1:3,4)./1000;
            pose_rel2 = estimateOdometry(K_pyr, img_prev_pyr{n}, img_curr_pyr{n}, dep_prev_pyr{n}./1000, dep_curr_pyr{n}./1000, pose_rel_);
            pose_rel = pose_rel2;
            pose_rel(1:3,4) = pose_rel(1:3,4)*1000;
            weightMap = zeros(size(img_prev_pyr{1}));
        else
            [ro,co] = size(img_prev_pyr{n});
            [xGridAll, yGridAll] = meshgrid(1:co, 1:ro);
            PixAll = [xGridAll(:) yGridAll(:)];
            
            rgbMatPrv = [reshape(img_prev_pyr{n}(:,:,1),[],1) reshape(img_prev_pyr{n}(:,:,1),[],1) reshape(img_prev_pyr{n}(:,:,1),[],1)];
            rgbMatCur = [reshape(img_curr_pyr{n}(:,:,1),[],1) reshape(img_curr_pyr{n}(:,:,1),[],1) reshape(img_curr_pyr{n}(:,:,1),[],1)];
            
            [XYZ_prv] = GetXYZFromDepth(K_pyr, PixAll,dep_prev_pyr{n}(:));
            [XYZ_cur] = GetXYZFromDepth(K_pyr, PixAll,dep_curr_pyr{n}(:));
            
            xyzMat_prv(:,:,1) = reshape(XYZ_prv(:,1), ro, co);xyzMat_prv(:,:,2) = reshape(XYZ_prv(:,2), ro, co);xyzMat_prv(:,:,3) = reshape(XYZ_prv(:,3), ro, co);
            xyzMat_cur(:,:,1) = reshape(XYZ_cur(:,1), ro, co);xyzMat_cur(:,:,2) = reshape(XYZ_cur(:,2), ro, co);xyzMat_cur(:,:,3) = reshape(XYZ_cur(:,3), ro, co);
            if 1
                ptCloudPrv = pointCloud(XYZ_prv(~isnan(XYZ_prv(:,3)),:),'Color',rgbMatPrv(~isnan(XYZ_prv(:,3)),:));
                ptCloudCur = pointCloud(XYZ_cur(~isnan(XYZ_cur(:,3)),:),'Color',rgbMatCur(~isnan(XYZ_cur(:,3)),:));
            else
                ptCloudPrv = pointCloud(xyzMat_prv,'Color',cat(3,img_prev_pyr{n},img_prev_pyr{n},img_prev_pyr{n}));
                ptCloudCur = pointCloud(xyzMat_cur,'Color',cat(3,img_curr_pyr{n},img_curr_pyr{n},img_curr_pyr{n}));
            end
            
            fixed1 = [];
            fixed1.ptcloud = ptCloudPrv; %removeInvalidPoints(pcdownsample(ptCloud2, 'random', percen));
            fixed1.image = img_prev_pyr{n}; %ptCloudPrv.Color(:,:,1);  % rgb2gray(ptCloudPrv.Color);
            moving1 = ptCloudCur; %removeInvalidPoints(pcdownsample(ptCloud1, 'random', percen));
            % make rkhs registration object
            dvo = rgbd_dvo((pose_rel), K_pyr);
            dvo.set_ptclouds(fixed1, moving1);
%             dvo.tform = tform12;
            dvo.align();
            pose_rel_ = dvo.tform.T';
            %
            weightMap = zeros(size(img_prev_pyr{1}));
            
            
        end
    end
    % intermediate results
    if isVerbose
        
        [warped_image1, valid_mask1] = warpImage(img_curr, dep_prev, pose_init, K, K);
        error1 = mean((warped_image1(valid_mask1) - img_prev(valid_mask1)).^2);
        [warped_image, valid_mask] = warpImage(img_curr, dep_prev, pose_rel, K, K);
        error = mean((warped_image(valid_mask) - img_prev(valid_mask)).^2);
        disp(['visual consistency score in level ' num2str(n) ' is ' num2str(error)]);
        
        figure(1);
        imshow(abs(warped_image - img_prev).*valid_mask, []);
        figure,subplot(1,2,1);imshowpair(warped_image1 ,  img_prev);title('before'); subplot(1,2,2),imshowpair(warped_image ,  img_prev);title('after');
    end
    
    % increse the focal length
    if 0
        K_pyr(1:2, :) = K_pyr(1:2, :) * 2;
    else
        if n ~= 1
%             K_pyr(1, 1) = K_pyr(1, 1) * (1/0.8);
%             K_pyr(2, 2) = K_pyr(2, 2) * (1/0.8);
            
            K_pyr(1, 1) = K_pyr(1, 1) * (2);
            K_pyr(2, 2) = K_pyr(2, 2) * (2);
            [height_next, width_next] = size(dep_prev_pyr{n-1});
            if 0
                K_pyr(1,3) = (1 + width_next) / 2;
                K_pyr(2,3) = (1 + height_next) / 2;
            else
%                 K_pyr(1, 3) = K_pyr(1, 3) * (1/0.8);
%                 K_pyr(2, 3) = K_pyr(2, 3) * (1/0.8);
                
                K_pyr(1, 3) = K_pyr(1, 3) * (2);
                K_pyr(2, 3) = K_pyr(2, 3) * (2);
            end
        end
    end
end


[UU,SS,VV] = svd(pose_rel(1:3,1:3));
pose_rel(1:3,1:3) = UU*VV';

% get the final score
[warped_image, valid_mask, warped_z, validMask] = warpImage(img_curr, dep_prev, pose_rel, K, K);
score = mean((warped_image(valid_mask) - img_prev(valid_mask)).^2);

if isVerbose
    disp(['The fitting score is ' num2str(error)]);
end

end
function [pose_rel2, error, weightMap] = iterPose(pose_rel, K_pyr, img_curr_pyr, img_prev_pyr, dep_curr_pyr, dep_prev_pyr, n,InParas,initial_sigma,default_dof, para)

[height, width] = size(dep_prev_pyr{n});

% get valid point clouds in the previous frame
[pointclouds, valid_mask] = reprojectDepthImage(dep_prev_pyr{n}, K_pyr);

% warp pointclouds and prune invalid points
warped_pointclouds = warpPointCloud(pointclouds, pose_rel);
[warped_img_coordinates, valid_points, ind] = projectPointCloud(warped_pointclouds, K_pyr, height, width);

warped_pointclouds = warped_pointclouds(valid_points, :);
pointclouds = pointclouds(valid_points, :);

% spatial gradient in the current frame
[Gx_curr, Gy_curr] = imgradientxy(img_curr_pyr{n}, 'CentralDifference');
Gx_curr = interp2(Gx_curr, warped_img_coordinates(:, 1), warped_img_coordinates(:, 2), 'linear', 0);
Gy_curr = interp2(Gy_curr, warped_img_coordinates(:, 1), warped_img_coordinates(:, 2), 'linear', 0);
Gs_curr = cat(2, Gx_curr, Gy_curr);

% temporal visual difference
Gt_prev = img_prev_pyr{n};
Gt_prev = Gt_prev(valid_mask);
Gt_prev = Gt_prev(valid_points);
Gt_curr = img_curr_pyr{n};
Gt_curr = interp2(Gt_curr, warped_img_coordinates(:, 1), warped_img_coordinates(:, 2), 'linear', 0);
Gt = Gt_curr - Gt_prev;

% calculate the warping Jacobian
warp_jacobian = calculateWarpingJacobian(warped_pointclouds, pointclouds, pose_rel, K_pyr);

% calculate the compositive jacobian
comp_jacobian = squeeze(sum(bsxfun(@times, Gs_curr, warp_jacobian), 2));
weightMap = ones(height, width);
if 1
    % adjust the importance of each residual elements, when required
    % necessary pruning of bad correspondences
    % compute the weight
    variance = computeResidualVariance(Gt, para);
    if 0
        
        I22 = WarpImage(img_prev_pyr{n},dep_curr_pyr{n},inv(pose_rel),InParas);
        r = abs(img_curr_pyr{n} - I22);
        r = immultiply(r, I22 > 0);
        [delta iteration]= TDistributionScaleEstimator(initial_sigma,default_dof,r);
        %Use t-distribution
        W = TDistributionInfluenceFunction(r,delta,default_dof);
        weightMap1 =reshape(W, size(img_prev_pyr{n},2), size(img_prev_pyr{n},1))';
        weights = weightMap1(ind);
        weights0 = sqrt((para + 1) ./ (para + Gt.^2./variance));
    else
        
        
        weights = sqrt((para + 1) ./ (para + Gt.^2./variance));
    end
    comp_jacobian = bsxfun(@times, comp_jacobian, weights);
    Gt = Gt.*weights;
    
    weightMap = zeros(height, width);
    weightMap(ind) = weights;
end

% calculate the increment motion
increment = -(comp_jacobian'*comp_jacobian)\(comp_jacobian'*Gt);

% get the current relative pose
increment = twistexp(increment);
pose_rel2 = increment * pose_rel;
[warped_image1, valid_mask1] = warpImage(img_curr_pyr{n}, dep_prev_pyr{n}, pose_rel, K_pyr, K_pyr);
[warped_image, valid_mask] = warpImage(img_curr_pyr{n}, dep_prev_pyr{n}, pose_rel2, K_pyr, K_pyr);

error0 = mean((warped_image1(valid_mask1) - img_prev_pyr{n}(valid_mask1)).^2);
error1 = mean((warped_image(valid_mask) - img_prev_pyr{n}(valid_mask)).^2);
if 0
    figure,subplot(1,2,1);imshowpair(warped_image1 ,  img_prev_pyr{n});title('before'); subplot(1,2,2),imshowpair(warped_image ,  img_prev_pyr{n});title('after');
    
    [warped_image2, valid_mask2] = warpImage(img_prev_pyr{n}, dep_curr_pyr{n}, inv(pose_rel2), K_pyr, K_pyr);
%     error1 = mean((warped_image(valid_mask) - img_prev_pyr{n}(valid_mask)).^2);
    error2 = mean((warped_image2(valid_mask2) - img_curr_pyr{n}(valid_mask2)).^2);
    error = error1  + error2;
else
    error = error1;
end
end
function [XYZ] = GetXYZFromDepth(intrMat, Pix,depthList)
metricPrevPtCcsGT = intrMat\HomoCoord(Pix',1);
metricPrevPtCcsGT = normc(metricPrevPtCcsGT);
if 1
    if 0
        if 0
            scaleAllGT = depthListGT(inlierId)./metricPrevPtCcsGT(3,:)';
        else
            scaleAllGT = depthListGT(:)./metricPrevPtCcsGT(3,:)';
        end
    else
        scaleAllGT = depthList./metricPrevPtCcsGT(3,:)';
    end
else
    dispGTComp = dispList(inlierId) - disparityErrorRound;
    depthListGTComp = intrMat(1,1).*norm(obj.camModel.transVec1To2)./(dispGTComp + (princpPtR(1) - princpPtL(1)));
    scaleAllGT = depthListGTComp./metricPrevPtCcsGT(3,:)';
end

XYZ = [repmat(scaleAllGT',3,1).*metricPrevPtCcsGT]';


end