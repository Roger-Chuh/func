function [stereoParamRGB] = CheckL2RGBReprojection(imgListL, imgListRGB, camParamRGB, cbcXYRGB, cbGridL, cbGridRGB, markerIdList, xyzCellGlobal, configL, configRGB)
global cfg inputDir0

height = cfg.img_height;
height_rgb = cfg.img_height_rgb;
width = cfg.img_width;

imgL = imread(imgListL{1});

% [nr, nc, nn] = size(imgL);


camParamL.foc = [cfg.KK_newL(1,1); cfg.KK_newL(2,2)];
camParamL.cen = [cfg.KK_newL(1,3); cfg.KK_newL(2,3)];
camParamL.alpha = 0;
camParamL.kc = zeros(5,1);


       

gIdLeft = cfg.Corner{1,1};
cbcXYL = cell(length(gIdLeft),1);
% 构造左相机在标定板坐标系中优化后的pose
for k = 1 : length(gIdLeft)
    pt3d_temp = [cbGridL{k}'];
    pt3d_temp(:,3) = 0;
    [poseVec_temp, outlierIdtemp] = posest(cfg.Corner{k,2}, pt3d_temp,  0.95, cfg.KK_newL, 'repr_err');
    camParamL.rotVec(:,k) = poseVec_temp(1:3);
    camParamL.tranVec(:,k) = poseVec_temp(4:6);
    cbcXYL{k,1} = cfg.Corner{k,2}';
end



goodIdL = gIdLeft;

goodId = goodIdL;

% 以fix内参和畸变的方式优化左目和RGB之间的外参
stereoParam = CbCalibStereo2(camParamL, camParamRGB, cbcXYL, cbcXYRGB, cbGridL, configL, configRGB);



stereoParamRGB = stereoParam;


% load预先保存的remap矩阵

load(fullfile(inputDir0,'rectMap.mat'));


if cfg.isRotate
    % 把横置双目的外参转为纵置双目的外参
    t_ = rotz(90)*stereoParam.transVecRef;
    r_ = rotz(90)*rodrigues(stereoParam.rotVecRef);
    r__ = [-r_(:,2) r_(:,1) r_(:,3)];
else
    
    r__ = rodrigues(stereoParam.rotVecRef);
    t_ = stereoParam.transVecRef;
end

cfg.L2RGB = [r__ t_;0 0 0 1];

KL = [stereoParam.focLeft(1) 0 stereoParam.cenLeft(1);0 stereoParam.focLeft(2) stereoParam.cenLeft(2);0 0 1];
KR = [stereoParam.focRight(1) 0 stereoParam.cenRight(1);0 stereoParam.focRight(2) stereoParam.cenRight(2);0 0 1];

if cfg.isRotate
    % 把横置图像的内参等效转成纵置方向的内参
    KL_= [KL(2,2) 0 1+height-KL(2,3); 0 KL(1,1) KL(1,3); 0 0 1];
    KR_= [KR(2,2) 0 1+height_rgb-KR(2,3); 0 KR(1,1) KR(1,3); 0 0 1];
else
    KL_ =KL;
    KR_ =KR;
end
cbcXYLL = cbcXYL;
cbcXYRR = cbcXYRGB;



imgLst11 = imgListRGB;


kr1_ = stereoParam.kcRight(1);
kr2_ = stereoParam.kcRight(2);
kr3_ = stereoParam.kcRight(5);
pr1_ = stereoParam.kcRight(3);
pr2_ = stereoParam.kcRight(4);





for i = 1 : length(goodId) 
    
    
    [xr] = normalize_pixel(cbcXYRR{(i)},stereoParam.focRight,stereoParam.cenRight,stereoParam.kcRight,stereoParam.alphaRight);
    [xl] = normalize_pixel(cbcXYLL{(i)},stereoParam.focLeft,stereoParam.cenLeft,stereoParam.kcLeft,stereoParam.alphaLeft);
    
    
    xll = pflat((KL)*pextend(xl));
    xrr = pflat((KR)*pextend(xr));
    if cfg.isRotate
        % 把横置图像的坐标等效转成纵置方向的坐标
        xrr_ = [1 + height_rgb - xrr(2,:); xrr(1,:)];
        xll_ = [1 + height - xll(2,:); xll(1,:)];
    else
        xrr_ = xrr(1:2,:);
        xll_ = xll(1:2,:);
    end
    
    
    
    
    pixRectR = xrr(1:2,:)';
    
    pixRectL = xll(1:2,:)';
    
    
    
    if 0
        
        imgL__1 = imread(cfg.Corner{i,3});
        if cfg.isRotate

            rectImgLtemp = uint8(interp2(matX, matY, double(imgL__1(:,:,1)),reMapRX,reMapRY));
        else

            rectImgLtemp = uint8(interp2(matX, matY, double(imgL__1(:,:,1)),reMapLX,reMapLY));
        end
        rectImgRtemp = undistortimage(imread(imgLst11{i}), stereoParam.focRight, KR(1,3), KR(2,3), kr1_, kr2_, kr3_, pr1_, pr2_);
        figure,subplot(1,2,1);imshow(rectImgRtemp);hold on;plot(pixRectR(:,1), pixRectR(:,2),'.r'); subplot(1,2,2),imshow(imread(imgLst11{i}));hold on;plot(cbcXYRR{(i)}(1,:), cbcXYRR{(i)}(2,:),'.r')
        figure,subplot(1,2,1);imshow(rectImgLtemp);hold on;plot(pixRectL(:,1), pixRectL(:,2),'.r'); subplot(1,2,2),imshow(imgL__1);hold on;plot(cbcXYLL{(i)}(1,:), cbcXYLL{(i)}(2,:),'.r')
    end
    
    
    rt = [rodrigues(stereoParam.optPose(1:3,(i))) stereoParam.optPose(4:6,(i));0 0 0 1];
    
    
    
    
    
    pt_3d_aruto = cbGridL{i}';
    pt_3d_aruto(:,3) = 0;
    
    
    
    cbPt = rt(1:3,1:3)*pt_3d_aruto' + repmat(rt(1:3,4),1,size(pt_3d_aruto,1));
    
    ptIcs_L = pflat(KL*cbPt);  ptIcs_L = ptIcs_L(1:2,:)';
    
    if 0
        imgL__1 = imread(cfg.Corner{i,3});
        if cfg.isRotate

            rectImgLtemp = uint8(interp2(matX, matY, double(imgL__1(:,:,1)),reMapRX,reMapRY));
        else

            rectImgLtemp = uint8(interp2(matX, matY, double(imgL__1(:,:,1)),reMapLX,reMapLY));
        end
        figure,imshow(rectImgLtemp);hold on;plot(pixRectL(:,1), pixRectL(:,2),'or');plot(ptIcs_L(:,1),ptIcs_L(:,2),'xg')
    end
    if cfg.isRotate
        cbPt_ = [-cbPt(2,:); cbPt(1,:); cbPt(3,:)];
    else
        cbPt_ = cbPt;
    end
    ptIcs = TransformAndProject(cbPt', KR, rodrigues(stereoParam.rotVecRef), stereoParam.transVecRef);
    
    
    ptIcs_ = TransformAndProject(cbPt_', KR_, r__, t_);
    
    
    
    if 0
        
        ptIcs_L = TransformAndProject(cbPt_', KL_, eye(3), [0 0 0]');
        ptIcs_L0 = TransformAndProject(pt_3d_aruto, KL, rt(1:3,1:3), rt(1:3,4));
        reprojErr_L = xll_' - ptIcs_L;
        if cfg.isRotate
            ptIcs_L00 = [1 + height - ptIcs_L0(:,2) ptIcs_L0(:,1)];
        else
            ptIcs_L00 = ptIcs_L0;
        end
        reprojErr_L0 = ptIcs_L00 - ptIcs_L;
        reprojErr_L1 = xll(1:2,:)' - ptIcs_L0;
        
        [~, reproj1] = NormalizeVector(reprojErr_L);
        [~, reproj2] = NormalizeVector(reprojErr_L1);
        reprojDiff = reproj1 - reproj2;
        
        
        rectImgRtemp = undistortimage(imread(imgLst11{i}), stereoParam.focRight, KR(1,3), KR(2,3), kr1_, kr2_, kr3_, pr1_, pr2_);
        
        figure,imshow(rectImgRtemp);hold on;plot(xrr(1,:),xrr(2,:),'or');plot(ptIcs(:,1),ptIcs(:,2),'.g')
        
        
        imgL__1 = imread(cfg.Corner{i,3});
        if cfg.isRotate
            rectImgLtemp = uint8(interp2(matX, matY, double(imgL__1(:,:,1)),reMapRX,reMapRY));
        else
            rectImgLtemp = uint8(interp2(matX, matY, double(imgL__1(:,:,1)),reMapLX,reMapLY));
        end
        rectImgRtemp = undistortimage(imread(imgLst11{i}), stereoParam.focRight, KR(1,3), KR(2,3), kr1_, kr2_, kr3_, pr1_, pr2_);
        if cfg.isRotate
            figure,subplot(1,2,1);imshow(imrotate(rectImgLtemp, -90));hold on;plot(xll_(1,:),xll_(2,:),'or');plot(ptIcs_L(:,1),ptIcs_L(:,2),'.g')
                   subplot(1,2,2);imshow(imrotate(rectImgRtemp, -90));hold on;plot(xrr_(1,:),xrr_(2,:),'or');plot(ptIcs_(:,1),ptIcs_(:,2),'.g')
        else
            figure,subplot(1,2,1);imshow(imrotate(rectImgLtemp, 0));hold on;plot(xll_(1,:),xll_(2,:),'or');plot(ptIcs_L(:,1),ptIcs_L(:,2),'.g')
                   subplot(1,2,2),imshow(imrotate(rectImgRtemp, 0));hold on;plot(xrr_(1,:),xrr_(2,:),'or');plot(ptIcs_(:,1),ptIcs_(:,2),'.g')
        end
    end
    err(i,1) = norm(mean(abs(xrr(1:2,:)'-ptIcs)));
    err_rotate(i,1) = norm(mean(abs(xrr_(1:2,:)'-ptIcs_)));
end

fprintf(sprintf('\n\n\n### average L to RGB error: %0.4f pixel ###\n\n\n',mean(err)));



cfg.RGBIntrMat = [KR];
cfg.RGBIntrMat(1,1) = KR(1,1) + cfg.ext_focal;
cfg.RGBIntrMat(2,2) = KR(2,2) + cfg.ext_focal;

% rectImgRtemp = undistortimage(imread(imgLst11{i}), stereoParam.focRight, KR(1,3), KR(2,3), kr1_, kr2_, kr3_, pr1_, pr2_);

end