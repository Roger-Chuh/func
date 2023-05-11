function [imgDiffMean, maskPrv, maskWarp, valid_mask_in_warp_prv, inValidMask_in_cur] = ImageAlignmentError(imgPrv, imgCur, prv_disp, cur_disp, intrMatL, intrMatRGB, baseline, optPose, varargin)


dispDiffThr = 1000; 3; 2; 1.5; 3; 8;

if 0
    imgPrv = histeq(imgPrv);
    
    imgCur = histeq(imgCur);
    
end

if (nargin == 8)
    draw = 0;
elseif (nargin == 9)
    draw = varargin{1};
else
    error('Too many input arguments');
end

cannyParam = [0.1 0.15];
if 0
    cannyParam = [0.05 0.1];
    % cannyParam = [0.2 0.25];
    cannyParam = [];
end
if 0
    se = strel('square', 7);
    se_ = strel('square', 3);
elseif 1
    se = strel('square', 3);
    se_ = strel('square', 1);
else
    se = strel('square', 7);
    se_ = strel('square', 3);
end
% 把rgb的内容渲染到depth坐标系
% 把6目的内容渲染到0坐标系
% imgCur->6目 / rgb
% imgPrv->0目 / depth
% optPose-> depth2rgb / 026 
    [warped_image1, valid_mask11, warped_z1, warped_z_mask1] = warpImage(double(imgCur), intrMatL(1)*baseline./prv_disp, optPose, intrMatL, intrMatRGB, double(imgPrv));
    
    dispDiff = abs(intrMatL(1)*baseline./warped_z1 - cur_disp);
    if 0
        inValidMask_in_cur1 = isnan(dispDiff) | (dispDiff > dispDiffThr & dispDiff < 20000);
    else
        inValidMask_in_cur1 = (dispDiff > dispDiffThr & dispDiff < 20000);
    end
    inValidMask_in_cur = inValidMask_in_cur1 | ~warped_z_mask1;
    if 0
        
        figure,subplot(1,3,1);imshowpair(isnan(dispDiff), ~warped_z_mask1);subplot(1,3,2);imshowpair(~warped_z_mask1, isnan(cur_disp));subplot(1,3,3);imshowpair(isnan(prv_disp), isnan(cur_disp));
        
        figure,subplot(1,2,1);imshow([warped_image1;double(imgPrv)], []);subplot(1,2,2);imshow([warped_z1; intrMatL(1)*baseline./cur_disp], [])
        
        figure,imshow([intrMatL(1)*baseline./prv_disp;warped_z1;intrMatL(1)*baseline./cur_disp], [])
        figure,imshow(abs(intrMatL(1)*baseline./warped_z1 - cur_disp), []);
        
    end
    
    valid_mask1 = valid_mask11 > 0 & warped_image1>0;
    valid_mask_in_warp_prv = valid_mask1;
    maskWarp = immultiply(double(warped_image1), valid_mask1);
    maskPrv = immultiply(double(imgPrv), valid_mask1);
    score1 = (( maskWarp(:)- maskPrv(:)).^2);
    
    valid_mask = imerode(valid_mask1, se);
    if 0
        [maskWarpEdgePix, maskWarpEdgeImg] = detectEdge(maskPrv, maskPrv, cannyParam);
        
    else
        if 1
            [maskWarpEdgePix_, maskWarpEdgeImg_] = detectEdge(maskPrv, maskPrv, cannyParam);
        else
            [maskWarpEdgePix_, maskWarpEdgeImg_] = detectEdge(maskWarp, maskWarp, cannyParam);
        end
        maskWarpEdgeImg = imdilate(maskWarpEdgeImg_, se_);
        [y, x] = ind2sub(size(maskWarpEdgeImg), find(maskWarpEdgeImg(:) > 0));
        maskWarpEdgePix = [x y];
    end
    
    ind = sub2ind(size(maskWarpEdgeImg), round(maskWarpEdgePix(:,2)), round(maskWarpEdgePix(:,1)));
    if 0
        validPix = maskWarpEdgePix(valid_mask(ind)>0,:);
    elseif 0
        validPix = maskWarpEdgePix(valid_mask(ind)>0 & maskPrv(ind)>0,:);
    else
        validPix = maskWarpEdgePix(valid_mask(ind)>0 & maskWarp(ind)>0,:);
    end
    
    indValid = sub2ind(size(maskWarpEdgeImg), round(validPix(:,2)), round(validPix(:,1)));
    imgDiff = abs(maskWarp(indValid) - maskPrv(indValid));
    
    if isempty(imgDiff)
        imgDiff = 2000000;
    end
    
    imgDiffMean = mean(imgDiff);
 
if  draw  % 1  % draw
    if draw == 1
        figure,imshowpair(maskPrv, maskWarp);hold on;plot(validPix(:,1), validPix(:,2), '.r'); title(sprintf('gray value diff: %0.3f',imgDiffMean));
    elseif draw == 2
        figure,imagesc(abs(maskPrv - maskWarp)); %hold on;plot(validPix(:,1), validPix(:,2), '.r');
        colorbar;
        title(sprintf('gray value diff: %0.3f  |  image diff: %0.3f | trans: %0.3fm',imgDiffMean, mean(abs(maskPrv(:) - maskWarp(:))), norm(optPose(1:3,4))));
    else
        figure,subplot(2,1,1),imshowpair(maskPrv, maskWarp);%hold on;plot(validPix(:,1), validPix(:,2), '.r');
        subplot(2,1,2),imagesc(abs(maskPrv - maskWarp));colorbar;
%         colorbar;
        title(sprintf('gray value diff: %0.3f  |  image diff: %0.3f | trans: %0.3fm',imgDiffMean, mean(abs(maskPrv(:) - maskWarp(:))), norm(optPose(1:3,4))));
    end
end     


end