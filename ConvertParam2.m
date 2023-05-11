function ConvertParam2(inputDir0, paraDir, stereoParam, imgSize, isStereo)

global cfg
% load('calib_zed.mat');

% if cfg.img_width >= 1920
if cfg.img_width > 1280
    imgSize0 = imgSize;
    cheat = 1;
    imgSize = [720 1280];
    scale = imgSize0(2)/1280;
else
    imgSize0 = imgSize;
    cheat = 0;
    scale = 1;
end


if imgSize(2) <= 640
    load(fullfile(paraDir, '640.mat'));
    sampledX1 = newXSample0;
    sampledY1 = newYSample0;
    
    sampledX2 = newXSample;
    sampledY2 = newYSample;
    
elseif imgSize(2) <= 1280 && imgSize(2) > 640
    load(fullfile(paraDir, '1280.mat'));
    sampledX1 = newXSample0;
    sampledY1 = newYSample0;
    
    sampledX2 = newXSample;
    sampledY2 = newYSample;
elseif imgSize(2) <= 1920 && imgSize(2) > 1280
    load(fullfile(paraDir, '1920.mat'));
    sampledX1 = newXSample0;
    sampledY1 = newYSample0;
    
    sampledX2 = newXSample;
    sampledY2 = newYSample;
elseif imgSize(2) > 3500
    load(fullfile(paraDir, '4k.mat'));
    sampledX1 = newXSample0;
    sampledY1 = newYSample0;
    
    sampledX2 = newXSample;
    sampledY2 = newYSample;
else
    sbghd = 1;
end




% cheat = 0;

% inputDir0 = pwd;


% imgSize0 = [1080 1920];

% imgSize0 = imgSize;

transVec = cfg.LRBaseline; %stereoParam.transVecRef;
% imgSize = [720 1280];
imgSizeHalf = imgSize;
imgSizeHalf(2) = imgSize(2)/2;


% [~, ~, rotMatL, rotMatR,  intrMatNewL, intrMatNewR] = GetRectifyParam2(stereoParam, imgSize);


%  cfg.lr_param.rectParamL = rectParamL;
%     cfg.lr_param.rectParamR = rectParamR;
if isStereo
    rotMatL = cfg.lr_param.rotMatLeft;
    rotMatR = cfg.lr_param.rotMatRight;
    intrMatNewL = cfg.lr_param.intrMatLeftNew;
    intrMatNewR = cfg.lr_param.intrMatRightNew;
    
    intrMatOldL = [stereoParam.focLeft(1) 0 stereoParam.cenLeft(1); 0 stereoParam.focLeft(2) stereoParam.cenLeft(2); 0 0 1];
    intrMatOldR = [stereoParam.focRight(1) 0 stereoParam.cenRight(1); 0 stereoParam.focRight(2) stereoParam.cenRight(2); 0 0 1];
    
    kcL = stereoParam.kcLeft;
    kcR = stereoParam.kcRight;
    
    intrMatNewL0 = intrMatNewL;
    intrMatNewR0 = intrMatNewR;
    if cheat
        intrMatNewL0 = intrMatNewL;
        intrMatNewL = intrMatNewL./scale;
        intrMatNewL(3,3) = 1;
        
        intrMatNewR0 = intrMatNewR;
        intrMatNewR = intrMatNewR./scale;
        intrMatNewR(3,3) = 1;
        
        intrMatOldL0 = intrMatOldL;
        intrMatOldL = intrMatOldL./scale;
        intrMatOldL(3,3) = 1;
        
        intrMatOldR0 = intrMatOldR;
        intrMatOldR = intrMatOldR./scale;
        intrMatOldR(3,3) = 1;
    end
    
    [matX, matY] = meshgrid(1:imgSize(2), 1:imgSize(1));
    pixRectL = remapRect([matX(:) matY(:)],  intrMatOldL, intrMatNewL ,rotMatL, kcL);
    pixRectR = remapRect([matX(:) matY(:)],  intrMatOldR, intrMatNewR ,rotMatR, kcR);
    reMapLX = reshape(pixRectL(:,1), imgSize);
    reMapLY = reshape(pixRectL(:,2), imgSize);
    reMapRX = reshape(pixRectR(:,1), imgSize);
    reMapRY = reshape(pixRectR(:,2), imgSize);
    
    if 1
        rectImgLtemp = uint8(interp2(matX, matY, double(cfg.imgL(:,:,1)),reMapLX,reMapLY));
        rectImgRtemp = uint8(interp2(matX, matY, double(cfg.imgR(:,:,1)),reMapRX,reMapRY));
    end
    
    if ~cfg.isRotate
        imwrite([imresize(uint8(rectImgLtemp),[240 320]), imresize(uint8(rectImgRtemp),[240 320])], fullfile(inputDir0,sprintf('pinhole_%05d.png',1)));
    else
        imwrite([imrotate(imresize(uint8(rectImgRtemp),[240 320]),-90), imrotate(imresize(uint8(rectImgLtemp),[240 320]),-90)], fullfile(inputDir0,sprintf('pinhole_%05d.png',1)));
    end
    
    if ~cfg.isRotate
        
        [LutVecOrig2RectX, LutVecOrig2RectY, LutVecOrig2RectXNum, LutVecOrig2RectYNum,Orig2RectNameMatX,Orig2RectNameMatY] = genLutData(imgSizeHalf,intrMatOldL,kcL,intrMatNewL,rotMatL, 1, 'L', sampledX1, sampledY1);
        [LutVecRect2OrigX, LutVecRect2OrigY, LutVecRect2OrigXNum, LutVecRect2OrigYNum,Rect2OrigNameMatX,Rect2OrigNameMatY] = genLutData(imgSize,intrMatOldL,kcL,intrMatNewL,rotMatL, 0, 'L', sampledX2, sampledY2);
        makeLut(imgSize, 'L',inputDir0, LutVecOrig2RectXNum,LutVecRect2OrigXNum,LutVecRect2OrigX,LutVecOrig2RectX,Orig2RectNameMatX,Rect2OrigNameMatX,LutVecOrig2RectY,LutVecRect2OrigY,LutVecOrig2RectYNum);
        
        xRect2OrigLT = reMapLX'; yRect2OrigLT = reMapLY';
        mapVecL = [size(reMapLX,2);size(reMapLX,1);xRect2OrigLT(:); yRect2OrigLT(:)];
        mapVecL(3:end) = round(mapVecL(3:end))-1;
        mapVecL(mapVecL < 0) = 0;
        fida = fopen(fullfile(inputDir0,  sprintf('invLutDecL_%d_%d.txt',cfg.img_width, cfg.img_height)),'w');
        fprintf(fida,'%d\n',round(mapVecL));
        fclose(fida);
        
        
        [LutVecOrig2RectXR, LutVecOrig2RectYR, LutVecOrig2RectXNumR, LutVecOrig2RectYNumR,Orig2RectNameMatXR,Orig2RectNameMatYR] = genLutData(imgSizeHalf,intrMatOldR,kcR,intrMatNewR,rotMatR, 1, 'R', sampledX1, sampledY1);
        [LutVecRect2OrigXR, LutVecRect2OrigYR, LutVecRect2OrigXNumR, LutVecRect2OrigYNumR,Rect2OrigNameMatXR,Rect2OrigNameMatYR] = genLutData(imgSize,intrMatOldR,kcR,intrMatNewR,rotMatR, 0, 'R', sampledX2, sampledY2);
        makeLut(imgSize, 'R',inputDir0, LutVecOrig2RectXNumR,LutVecRect2OrigXNumR,LutVecRect2OrigXR,LutVecOrig2RectXR,Orig2RectNameMatXR,Rect2OrigNameMatXR,LutVecOrig2RectYR,LutVecRect2OrigYR,LutVecOrig2RectYNumR);
        
        xRect2OrigRT = reMapRX'; yRect2OrigRT = reMapRY';
        mapVecR = [size(reMapRX,2);size(reMapRX,1);xRect2OrigRT(:); yRect2OrigRT(:)];
        mapVecR(3:end) = round(mapVecR(3:end))-1;
        mapVecR(mapVecR < 0) = 0;
        fida = fopen(fullfile(inputDir0,  sprintf('invLutDecR_%d_%d.txt',cfg.img_width, cfg.img_height)),'w');
        fprintf(fida,'%d\n',round(mapVecR));
        fclose(fida);
        
        fid1 = fopen(fullfile(inputDir0,'cam_param1.txt'),'w');%????
        
        r_stereo = eye(3);
        r_stereo_t = r_stereo';
        r_stereo_t_vec = r_stereo_t(:);
        
        fprintf(fid1,'%d %d %d',imgSize0(2), imgSize0(1), cfg.switch_lr);
        fprintf(fid1, '\n');
        fprintf(fid1,'pinhole');
        fprintf(fid1, '\n');
        %     fprintf(fid1, '%d\n',2);
        fprintf(fid1,'%0.6f %0.6f %0.6f %0.6f\n',intrMatNewL0(1,1), intrMatNewL0(2,2), intrMatNewL0(1,3),intrMatNewL0(2,3));
        fprintf(fid1,'%0.6f %0.6f %0.6f %0.6f\n',intrMatNewR0(1,1), intrMatNewR0(2,2), intrMatNewR0(1,3),intrMatNewR0(2,3));
        %     fprintf(fid1,'%06f %06f %06f %06f %06f %06f %06f %06f %06f %06f %06f %06f\n',r_stereo_t_vec,-norm(transVec),0, 0);
        fprintf(fid1,'%06f %06f %06f %06f %06f %06f\n',0,0,0,-norm(transVec),0, 0);
        
        fclose(fid1);
        if ~cfg.is_kepler
            
            fid1 = fopen(fullfile(inputDir0,'cam_param.txt'),'w');%????
            
            r_stereo = eye(3);
            r_stereo_t = r_stereo';
            r_stereo_t_vec = r_stereo_t(:);
            
            fprintf(fid1,'%d %d %d',imgSize0(2), imgSize0(1), cfg.switch_lr);
            fprintf(fid1, '\n');
            fprintf(fid1,'pinhole');
            fprintf(fid1, '\n');
            fprintf(fid1, '%d\n',2);
            fprintf(fid1,'%0.6f %0.6f %0.6f %0.6f\n',intrMatNewL0(1,1), intrMatNewL0(2,2), intrMatNewL0(1,3),intrMatNewL0(2,3));
            fprintf(fid1,'%0.6f %0.6f %0.6f %0.6f\n',intrMatNewR0(1,1), intrMatNewR0(2,2), intrMatNewR0(1,3),intrMatNewR0(2,3));
            fprintf(fid1,'%06f %06f %06f %06f %06f %06f %06f %06f %06f %06f %06f %06f\n',r_stereo_t_vec,-norm(transVec),0, 0);
            %             fprintf(fid1,'%06f %06f %06f %06f %06f %06f\n',0,0,0,-norm(transVec),0, 0);
            
            fclose(fid1);
            
            
            
            
        end
    else
        
          [LutVecOrig2RectX, LutVecOrig2RectY, LutVecOrig2RectXNum, LutVecOrig2RectYNum,Orig2RectNameMatX,Orig2RectNameMatY] = genLutData(imgSizeHalf,intrMatOldL,kcL,intrMatNewL,rotMatL, 1, 'R', sampledX1, sampledY1);
        [LutVecRect2OrigX, LutVecRect2OrigY, LutVecRect2OrigXNum, LutVecRect2OrigYNum,Rect2OrigNameMatX,Rect2OrigNameMatY] = genLutData(imgSize,intrMatOldL,kcL,intrMatNewL,rotMatL, 0, 'R', sampledX2, sampledY2);
        makeLut(imgSize, 'R',inputDir0, LutVecOrig2RectXNum,LutVecRect2OrigXNum,LutVecRect2OrigX,LutVecOrig2RectX,Orig2RectNameMatX,Rect2OrigNameMatX,LutVecOrig2RectY,LutVecRect2OrigY,LutVecOrig2RectYNum);
        xRect2OrigRT = reMapRX'; yRect2OrigRT = reMapRY';
        mapVecR = [size(reMapRX,2);size(reMapRX,1);xRect2OrigRT(:); yRect2OrigRT(:)];
        mapVecR(3:end) = round(mapVecR(3:end))-1;
        mapVecR(mapVecR < 0) = 0;
        fida = fopen(fullfile(inputDir0,  sprintf('invLutDecR_%d_%d.txt',cfg.img_height, cfg.img_width)),'w');
        fprintf(fida,'%d\n',round(mapVecR));
        fclose(fida);
        
        [LutVecOrig2RectXR, LutVecOrig2RectYR, LutVecOrig2RectXNumR, LutVecOrig2RectYNumR,Orig2RectNameMatXR,Orig2RectNameMatYR] = genLutData(imgSizeHalf,intrMatOldR,kcR,intrMatNewR,rotMatR, 1, 'L', sampledX1, sampledY1);
        [LutVecRect2OrigXR, LutVecRect2OrigYR, LutVecRect2OrigXNumR, LutVecRect2OrigYNumR,Rect2OrigNameMatXR,Rect2OrigNameMatYR] = genLutData(imgSize,intrMatOldR,kcR,intrMatNewR,rotMatR, 0, 'L', sampledX2, sampledY2);
        makeLut(imgSize, 'L',inputDir0, LutVecOrig2RectXNumR,LutVecRect2OrigXNumR,LutVecRect2OrigXR,LutVecOrig2RectXR,Orig2RectNameMatXR,Rect2OrigNameMatXR,LutVecOrig2RectYR,LutVecRect2OrigYR,LutVecOrig2RectYNumR);
        xRect2OrigLT = reMapLX'; yRect2OrigLT = reMapLY';
        mapVecL = [size(reMapLX,2);size(reMapLX,1);xRect2OrigLT(:); yRect2OrigLT(:)];
        mapVecL(3:end) = round(mapVecL(3:end))-1;
        mapVecL(mapVecL < 0) = 0;
        fida = fopen(fullfile(inputDir0,  sprintf('invLutDecL_%d_%d.txt',cfg.img_height, cfg.img_width)),'w');
        fprintf(fida,'%d\n',round(mapVecL));
        fclose(fida);
        
        fid1 = fopen(fullfile(inputDir0,'cam_param1.txt'),'w');%????
        
        r_stereo = eye(3);
        r_stereo_t = r_stereo';
        r_stereo_t_vec = r_stereo_t(:);
        
        fprintf(fid1,'%d %d %d',imgSize0(2), imgSize0(1), cfg.switch_lr);
        fprintf(fid1, '\n');
        fprintf(fid1,'pinhole');
        fprintf(fid1, '\n');
        %     fprintf(fid1, '%d\n',2);
        if 0
            fprintf(fid1,'%0.6f %0.6f %0.6f %0.6f\n',intrMatNewL0(1,1), intrMatNewL0(2,2), intrMatNewL0(1,3),intrMatNewL0(2,3));
            fprintf(fid1,'%0.6f %0.6f %0.6f %0.6f\n',intrMatNewR0(1,1), intrMatNewR0(2,2), intrMatNewR0(1,3),intrMatNewR0(2,3));
        else
            fprintf(fid1,'%06f %06f %06f %06f',intrMatNewL0(2,2),intrMatNewL0(1,1),1+imgSize(1)-intrMatNewL0(2,3),intrMatNewL0(1,3));
            fprintf(fid1, '\n');
            fprintf(fid1,'%06f %06f %06f %06f',intrMatNewL0(2,2),intrMatNewL0(1,1),1+imgSize(1)-intrMatNewL0(2,3),intrMatNewL0(1,3));
            fprintf(fid1, '\n');
        end
        %     fprintf(fid1,'%06f %06f %06f %06f %06f %06f %06f %06f %06f %06f %06f %06f\n',r_stereo_t_vec,-norm(transVec),0, 0);
        fprintf(fid1,'%06f %06f %06f %06f %06f %06f\n',0,0,0,-norm(transVec),0, 0);
        
        fclose(fid1);
        
        if ~cfg.is_kepler
            
            fid1 = fopen(fullfile(inputDir0,'cam_param.txt'),'w');%????
            
            r_stereo = eye(3);
            r_stereo_t = r_stereo';
            r_stereo_t_vec = r_stereo_t(:);
            
            fprintf(fid1,'%d %d %d',imgSize0(2), imgSize0(1), cfg.switch_lr);
            fprintf(fid1, '\n');
            fprintf(fid1,'pinhole');
            fprintf(fid1, '\n');
            fprintf(fid1, '%d\n',2);
            if 0
                fprintf(fid1,'%0.6f %0.6f %0.6f %0.6f\n',intrMatNewL0(1,1), intrMatNewL0(2,2), intrMatNewL0(1,3),intrMatNewL0(2,3));
                fprintf(fid1,'%0.6f %0.6f %0.6f %0.6f\n',intrMatNewR0(1,1), intrMatNewR0(2,2), intrMatNewR0(1,3),intrMatNewR0(2,3));
            else
                fprintf(fid1,'%06f %06f %06f %06f',intrMatNewL0(2,2),intrMatNewL0(1,1),1+imgSize(1)-intrMatNewL0(2,3),intrMatNewL0(1,3));
                fprintf(fid1, '\n');
                fprintf(fid1,'%06f %06f %06f %06f',intrMatNewL0(2,2),intrMatNewL0(1,1),1+imgSize(1)-intrMatNewL0(2,3),intrMatNewL0(1,3));
                fprintf(fid1, '\n');
            end
            fprintf(fid1,'%06f %06f %06f %06f %06f %06f %06f %06f %06f %06f %06f %06f\n',r_stereo_t_vec,-norm(transVec),0, 0);
%             fprintf(fid1,'%06f %06f %06f %06f %06f %06f\n',0,0,0,-norm(transVec),0, 0);
            
            fclose(fid1);
            
            
            
            
        end
        
        
    end
else
    %     rotMatL = eye(3); %cfg.lr_param.rotMatLeft;
    %     rotMatR = cfg.lr_param.rotMatRight;
    %     intrMatNewL = cfg.lr_param.intrMatLeftNew;
    %     intrMatNewR = cfg.lr_param.intrMatRightNew;
    
    
    intrMatRGBNew = [cfg.stereoParamRGB.foc(1) 0 cfg.stereoParamRGB.cen(1); 0 cfg.stereoParamRGB.foc(2) cfg.stereoParamRGB.cen(2); 0 0 1];
    focLeftRGB = cfg.stereoParamRGB.foc;
    cenLeftRGB = cfg.stereoParamRGB.cen;
    kcLeftRGB = cfg.stereoParamRGB.kc;
    alphaLeftRGB = cfg.stereoParamRGB.alpha;
    intrMatRGBNew(1,1) = focLeftRGB(1) + cfg.ext_focal;
    intrMatRGBNew(2,2) = focLeftRGB(2) + cfg.ext_focal;
    rotMatRGB = eye(3);
    
    
    intrMatRGB_Err = cfg.RGBIntrMat - intrMatRGBNew;
    
    
    intrMatRGBOld = [cfg.stereoParamRGB.foc(1) 0 cfg.stereoParamRGB.cen(1); 0 cfg.stereoParamRGB.foc(2) cfg.stereoParamRGB.cen(2); 0 0 1];
    
    
    %     intrMatOldL = [stereoParam.focLeft(1) 0 stereoParam.cenLeft(1); 0 stereoParam.focLeft(2) stereoParam.cenLeft(2); 0 0 1];
    %     intrMatOldR = [stereoParam.focRight(1) 0 stereoParam.cenRight(1); 0 stereoParam.focRight(2) stereoParam.cenRight(2); 0 0 1];
    
    kcRGB = cfg.stereoParamRGB.kc;
    %     kcR = stereoParam.kcRight;
    
    intrMatRGBNew0 = intrMatRGBNew;
    intrMatRGBOld0 = intrMatRGBOld;
    %     intrMatNewR0 = intrMatNewR;
    if cheat
        intrMatRGBNew0 = intrMatRGBNew;
        intrMatRGBNew = intrMatRGBNew0./scale;
        intrMatRGBNew(3,3) = 1;
        
        %         intrMatNewR0 = intrMatNewR;
        %         intrMatNewR = intrMatNewR./scale;
        %         intrMatNewR(3,3) = 1;
        %
        intrMatRGBOld0 = intrMatRGBOld;
        intrMatRGBOld = intrMatRGBOld0./scale;
        intrMatRGBOld(3,3) = 1;
        
        %         intrMatOldR0 = intrMatOldR;
        %         intrMatOldR = intrMatOldR./scale;
        %         intrMatOldR(3,3) = 1;
    end
    
    [matX, matY] = meshgrid(1:imgSize(2), 1:imgSize(1));
    pixRectRGB = remapRect([matX(:) matY(:)],  intrMatRGBOld, intrMatRGBNew ,rotMatRGB, kcRGB);
    %         pixRectR = remapRect([matX(:) matY(:)],  intrMatOldR, intrMatNewR ,rotMatR, kcR);
    reMapRGBX = reshape(pixRectRGB(:,1), imgSize);
    reMapRGBY = reshape(pixRectRGB(:,2), imgSize);
    %         reMapRX = reshape(pixRectR(:,1), imgSize);
    %         reMapRY = reshape(pixRectR(:,2), imgSize);
    
    if 1
        imgOutRGB = uint8(interp2(matX, matY, double(cfg.imgRGB(:,:,1)),reMapRGBX,reMapRGBY));
        %             rectImgRtemp = uint8(interp2(matX, matY, double(imgR(:,:,1)),reMapRX,reMapRY));
    end
    
    [LutVecOrig2RectX, LutVecOrig2RectY, LutVecOrig2RectXNum, LutVecOrig2RectYNum,Orig2RectNameMatX,Orig2RectNameMatY] = genLutData(imgSizeHalf,intrMatRGBOld,kcRGB,intrMatRGBNew,rotMatRGB, 1, 'L', sampledX1, sampledY1);
    [LutVecRect2OrigX, LutVecRect2OrigY, LutVecRect2OrigXNum, LutVecRect2OrigYNum,Rect2OrigNameMatX,Rect2OrigNameMatY] = genLutData(imgSize,intrMatRGBOld,kcRGB,intrMatRGBNew,rotMatRGB, 0, 'L', sampledX2, sampledY2);
    makeLut(imgSize, 'RGB',fullfile(inputDir0,'RGB'), LutVecOrig2RectXNum,LutVecRect2OrigXNum,LutVecRect2OrigX,LutVecOrig2RectX,Orig2RectNameMatX,Rect2OrigNameMatX,LutVecOrig2RectY,LutVecRect2OrigY,LutVecOrig2RectYNum);
    
    xRect2OrigRGBT = reMapRGBX'; yRect2OrigRGBT = reMapRGBY';
        mapVecRGB = [size(reMapRGBX,2);size(reMapRGBX,1);xRect2OrigRGBT(:); yRect2OrigRGBT(:)];
        mapVecRGB(3:end) = round(mapVecRGB(3:end))-1;
        mapVecRGB(mapVecRGB < 0) = 0;
        if ~cfg.isRotate
            fida = fopen(fullfile(inputDir0,'RGB',  sprintf('invLutDecRGB_%d_%d.txt',cfg.img_width, cfg.img_height)),'w');
        else
            fida = fopen(fullfile(inputDir0,'RGB',  sprintf('invLutDecRGB_%d_%d.txt',cfg.img_height, cfg.img_width)),'w');
        end
        fprintf(fida,'%d\n',round(mapVecRGB));
        fclose(fida);
    
    
     fid1 = fopen(fullfile(fullfile(inputDir0,'RGB'),'cam_param1.txt'),'w');%????
     fprintf(fid1,'%d %d %d',imgSize0(2), imgSize0(1), cfg.switch_lr);
     fprintf(fid1, '\n');
     fprintf(fid1,'pinhole');
     fprintf(fid1, '\n');
     if ~cfg.isRotate
         fprintf(fid1,'%06f %06f %06f %06f',intrMatRGBNew0(1,1),intrMatRGBNew0(2,2),intrMatRGBNew0(1,3),intrMatRGBNew0(2,3));
     else
         fprintf(fid1,'%06f %06f %06f %06f',intrMatRGBNew0(2,2),intrMatRGBNew0(1,1),1+imgSize(1)-intrMatRGBNew0(2,3),intrMatRGBNew0(1,3));
     end
     if isfield(cfg, 'L2RGB')
         fprintf(fid1, '\n');
         rotVecRGB = rodrigues(cfg.L2RGB(1:3,1:3));
         tVecRGB = cfg.L2RGB(1:3,4);
         fprintf(fid1,'%06f %06f %06f %06f %06f %06f',rotVecRGB(1),rotVecRGB(2),rotVecRGB(3),tVecRGB(1),tVecRGB(2),tVecRGB(3));
     end
     
     %             fprintf(fid1, '\n');
     %             fprintf(fid1,'%06f %06f %06f %06f',KK_newR(1,1),KK_newR(2,2),KK_newR(1,3),KK_newR(2,3));
     %             fprintf(fid1, '\n');
     %             fprintf(fid1,'%06f %06f %06f %06f %06f %06f',0,0,0,-norm(stereoParam.transVecRef),0,0);
     fprintf(fid1, '\n');
     fclose(fid1);
     if 0
         figure,imshow([imresize(rgb2gray(cfg.imgLRGB),[240 320]) imresize(uint8(imgOutRGB),[240 320])]);
     end
     if 1
         if 0
             saveas(gcf,fullfile(fullfile(inputDir0,'RGB'),sprintf('pinhole_mono__%05d.png',2)));
         else
             imwrite([imresize(rgb2gray(cfg.imgRGB),[240 320]) imresize(uint8(imgOutRGB),[240 320])],fullfile(fullfile(inputDir0,'RGB'),sprintf('pinhole_mono__%05d.png',2)));
         end
     end
     
     
      fid1 = fopen(fullfile(inputDir0,'cam_param.txt'),'w');%????
            
        r_stereo = eye(3);
        r_stereo_t = r_stereo';
        r_stereo_t_vec = r_stereo_t(:);
        
        fprintf(fid1,'%d %d %d',cfg.img_width, cfg.img_height, cfg.switch_lr);
        fprintf(fid1, '\n');
        if ~cfg.wide_angle
            fprintf(fid1,'pinhole');
        else
            fprintf(fid1,'fisheye');
        end
        fprintf(fid1, '\n');
        if cfg.is_kepler
            imgSizeRGB = size(imgOutRGB);
            R_RGB_t = cfg.L2RGB(1:3,1:3)';
            R_RGB_t_vec =R_RGB_t(:);
            
            
            fprintf(fid1, '%d\n',3);
            if ~cfg.isRotate
                fprintf(fid1,'%0.6f %0.6f %0.6f %0.6f\n',cfg.LRIntrMat(1,1), cfg.LRIntrMat(2,2), cfg.LRIntrMat(1,3),cfg.LRIntrMat(2,3));
                fprintf(fid1,'%0.6f %0.6f %0.6f %0.6f\n',cfg.LRIntrMat(1,1), cfg.LRIntrMat(2,2), cfg.LRIntrMat(1,3),cfg.LRIntrMat(2,3));
                fprintf(fid1,'%0.6f %0.6f %0.6f %0.6f\n',cfg.RGBIntrMat(1,1), cfg.RGBIntrMat(2,2), cfg.RGBIntrMat(1,3),cfg.RGBIntrMat(2,3));
            else
                fprintf(fid1,'%0.6f %0.6f %0.6f %0.6f\n',cfg.LRIntrMat(2,2), cfg.LRIntrMat(1,1), 1+cfg.img_height(1)-cfg.LRIntrMat(2,3),cfg.LRIntrMat(1,3));
                fprintf(fid1,'%0.6f %0.6f %0.6f %0.6f\n',cfg.LRIntrMat(2,2), cfg.LRIntrMat(1,1), 1+cfg.img_height(1)-cfg.LRIntrMat(2,3),cfg.LRIntrMat(1,3));
                fprintf(fid1,'%0.6f %0.6f %0.6f %0.6f\n',cfg.RGBIntrMat(2,2), cfg.RGBIntrMat(1,1), 1+imgSizeRGB(1)-cfg.RGBIntrMat(2,3),cfg.RGBIntrMat(1,3));
                
            end
            fprintf(fid1,'%06f %06f %06f %06f %06f %06f %06f %06f %06f %06f %06f %06f\n',r_stereo_t_vec,-cfg.LRBaseline,0, 0);
            fprintf(fid1,'%06f %06f %06f %06f %06f %06f %06f %06f %06f %06f %06f %06f\n',R_RGB_t_vec, cfg.L2RGB(1:3,4));
        else
            fprintf(fid1, '%d\n',2);
            if ~cfg.isRotate
                fprintf(fid1,'%0.6f %0.6f %0.6f %0.6f\n',cfg.LRIntrMat(1,1), cfg.LRIntrMat(2,2), cfg.LRIntrMat(1,3),cfg.LRIntrMat(2,3));
                fprintf(fid1,'%0.6f %0.6f %0.6f %0.6f\n',cfg.LRIntrMat(1,1), cfg.LRIntrMat(2,2), cfg.LRIntrMat(1,3),cfg.LRIntrMat(2,3));
            else
                fprintf(fid1,'%0.6f %0.6f %0.6f %0.6f\n',cfg.LRIntrMat(2,2), cfg.LRIntrMat(1,1), 1+cfg.img_height(1)-cfg.LRIntrMat(2,3),cfg.LRIntrMat(1,3));
                fprintf(fid1,'%0.6f %0.6f %0.6f %0.6f\n',cfg.LRIntrMat(2,2), cfg.LRIntrMat(1,1), 1+cfg.img_height(1)-cfg.LRIntrMat(2,3),cfg.LRIntrMat(1,3));
            end
%             fprintf(fid1,'%0.6f %0.6f %0.6f %0.6f\n',KK_new_rgb(1,1), KK_new_rgb(2,2), KK_new_rgb(1,3),KK_new_rgb(2,3));
            fprintf(fid1,'%06f %06f %06f %06f %06f %06f %06f %06f %06f %06f %06f %06f\n',r_stereo_t_vec,-cfg.LRBaseline,0, 0);
%             fprintf(fid1,'%06f %06f %06f %06f %06f %06f %06f %06f %06f %06f %06f %06f\n',R_RGB_t_vec,stereoParam_.T_RGB);
        end
        
    
    fclose(fid1);
     
     
end


end
function [LutVecCharX, LutVecCharY, lutVecX, lutVecY, nameMatX, nameMatY] = genLutData(imgSize,intrMatOld,kc,intrMatNew,rotMat, reverseMapping, whichCam, sampledX, sampledY)
nc = imgSize(2);
nr = imgSize(1);
lutSize = [length(sampledY) length(sampledX)];




[xMatSampled, yMatSampled] = meshgrid(sampledX, sampledY);
pixSampled = [xMatSampled(:) yMatSampled(:)];

if reverseMapping == 0
    pixOrigSampled = remapRect(pixSampled, intrMatOld, intrMatNew ,rotMat, kc);
else
    pixOrigSampled = Orig2Rect(pixSampled, intrMatOld, intrMatNew, rotMat, kc);
end

xOrigMat = reshape(pixOrigSampled(:,1), lutSize);
yOrigMat = reshape(pixOrigSampled(:,2), lutSize);


if reverseMapping == 0
    deltaLutX2 = xMatSampled - xOrigMat; %MakeOfst2(inValidLutX,deltaLutX,ones(size(validMat)));
    deltaLutY2 = yMatSampled - yOrigMat; %MakeOfst2(inValidLutY,deltaLutY,ones(size(validMat)));
    thr = 2^7;
else
    pixOrigSampled_floor = floor(pixOrigSampled);
    pixOrigSampled_ceil = ceil(pixOrigSampled);
    ind1 = find(pixOrigSampled_floor(:,1) >= 1 & pixOrigSampled_floor(:,1) <= nc & pixOrigSampled_floor(:,2) >= 1 & pixOrigSampled_floor(:,2) <= nr);
    ind2 = find(pixOrigSampled_ceil(:,1) >= 1 & pixOrigSampled_ceil(:,1) <= nc & pixOrigSampled_ceil(:,2) >= 1 & pixOrigSampled_ceil(:,2) <= nr);
    
    validMat = zeros(lutSize);
    validMat([ind1; ind2]) = 1;
    validMat([1 end],:) = 0;
    validMat(:,[1 end]) = 0;
    se = strel('square',[3]);
    validMat = imdilate(validMat,se);
    deltaLutX = xMatSampled - xOrigMat;
    deltaLutY = yMatSampled - yOrigMat;
    inValidLutX = deltaLutX .* (~validMat);
    inValidLutY = deltaLutY .* (~validMat);
    deltaLutX2 = MakeOfst2(inValidLutX,deltaLutX,validMat);
    deltaLutY2 = MakeOfst2(inValidLutY,deltaLutY,validMat);
    thr = 2^8;
end


deltaDeltaLutX2 = cumsum([deltaLutX2(1,:);diff(deltaLutX2)]);
deltaDeltaLutX2(deltaDeltaLutX2 > thr) = (thr-1);
deltaDeltaLutX2(deltaDeltaLutX2 < -thr) = (-thr+1);
xOrigMatRecovered = xMatSampled - deltaDeltaLutX2;


deltaDeltaLutY2 = cumsum([deltaLutY2(1,:);diff(deltaLutY2)]);
deltaDeltaLutY2(deltaDeltaLutY2 > thr) = (thr-1);
deltaDeltaLutY2(deltaDeltaLutY2 < -thr) = (-thr+1);
yOrigMatRecovered = yMatSampled - deltaDeltaLutY2;

[LutVecCharX, LutVecCharY,lutVecX,lutVecY,nameMatX, nameMatY] = ...
    BiLinearInterp(reverseMapping, xMatSampled, yMatSampled, xOrigMatRecovered, yOrigMatRecovered,imgSize);

end

function pixDist = remapRect(pixRect, KDistort, KRect, R, distCoeff)

alpha = 0;
pixRectHomo = [pixRect'; ones(1,size(pixRect,1))];
rays = inv(KRect)*pixRectHomo;


% Rotation: (or affine transformation):

rays2 = R'*rays;

x = [rays2(1,:)./rays2(3,:);rays2(2,:)./rays2(3,:)];


% Add distortion:
xd = apply_distortion(x,distCoeff);


% Reconvert in pixels:

px2_ = KDistort(1,1)*(xd(1,:) + alpha*xd(2,:)) + KDistort(1,3);
py2_ = KDistort(2,2)*xd(2,:) + KDistort(2,3);
pixDist = [px2_;py2_]';

end
function [xd,dxddk] = apply_distortion(x,k)


% Complete the distortion vector if you are using the simple distortion model:
length_k = length(k);
if length_k <5 ,
    k = [k ; zeros(5-length_k,1)];
end;


[m,n] = size(x);

% Add distortion:

r2 = x(1,:).^2 + x(2,:).^2;

r4 = r2.^2;

r6 = r2.^3;


% Radial distortion:

cdist = 1 + k(1) * r2 + k(2) * r4 + k(5) * r6;

if nargout > 1,
    dcdistdk = [ r2' r4' zeros(n,2) r6'];
end;


xd1 = x .* (ones(2,1)*cdist);



if nargout > 1,
    dxd1dk = zeros(2*n,5);
    dxd1dk(1:2:end,:) = (x(1,:)'*ones(1,5)) .* dcdistdk;
    dxd1dk(2:2:end,:) = (x(2,:)'*ones(1,5)) .* dcdistdk;
end;


% tangential distortion:

a1 = 2.*x(1,:).*x(2,:);
a2 = r2 + 2*x(1,:).^2;
a3 = r2 + 2*x(2,:).^2;

delta_x = [k(3)*a1 + k(4)*a2 ;
    k(3) * a3 + k(4)*a1];



if nargout > 1,
    ddelta_xdk = zeros(2*n,5);
    ddelta_xdk(1:2:end,3) = a1';
    ddelta_xdk(1:2:end,4) = a2';
    ddelta_xdk(2:2:end,3) = a3';
    ddelta_xdk(2:2:end,4) = a1';
end;

xd = xd1 + delta_x;

if nargout > 1,
    dxddk = dxd1dk + ddelta_xdk ;
    if length_k < 5,
        dxddk = dxddk(:,1:length_k);
    end;
end;


end

function pixRect = Orig2Rect(pix, intrMatOld, intrMatNew, R, kc)

[pixUndist] = normalize_pixel(pix',[intrMatOld(1,1);intrMatOld(2,2)],[intrMatOld(1,3);intrMatOld(2,3)],kc,0);
pixUndistHomo = [pixUndist; ones(1, size(pixUndist,2))];
pixUndistR = R*pixUndistHomo;
pixRect = intrMatNew*pixUndistR;
pixRect = [pixRect(1,:)./pixRect(3,:); pixRect(2,:)./pixRect(3,:)];
pixRect = pixRect(1:2,:)';

end
function [xn] = normalize_pixel(x_kk,fc,cc,kc,alpha_c)

%normalize
%
%[xn] = normalize_pixel(x_kk,fc,cc,kc,alpha_c)
%
%Computes the normalized coordinates xn given the pixel coordinates x_kk
%and the intrinsic camera parameters fc, cc and kc.
%
%INPUT: x_kk: Feature locations on the images
%       fc: Camera focal length
%       cc: Principal point coordinates
%       kc: Distortion coefficients
%       alpha_c: Skew coefficient
%
%OUTPUT: xn: Normalized feature locations on the image plane (a 2XN matrix)
%
%Important functions called within that program:
%
%comp_distortion_oulu: undistort pixel coordinates.

if nargin < 5,
    alpha_c = 0;
    if nargin < 4;
        kc = [0;0;0;0;0];
        if nargin < 3;
            cc = [0;0];
            if nargin < 2,
                fc = [1;1];
            end;
        end;
    end;
end;


% First: Subtract principal point, and divide by the focal length:
x_distort = [(x_kk(1,:) - cc(1))/fc(1);(x_kk(2,:) - cc(2))/fc(2)];

% Second: undo skew
x_distort(1,:) = x_distort(1,:) - alpha_c * x_distort(2,:);

if norm(kc) ~= 0,
    % Third: Compensate for lens distortion:
    xn = comp_distortion_oulu(x_distort,kc);
else
    xn = x_distort;
end;

end

function [x] = comp_distortion_oulu(xd,k);

%comp_distortion_oulu.m
%
%[x] = comp_distortion_oulu(xd,k)
%
%Compensates for radial and tangential distortion. Model From Oulu university.
%For more informatino about the distortion model, check the forward projection mapping function:
%project_points.m
%
%INPUT: xd: distorted (normalized) point coordinates in the image plane (2xN matrix)
%       k: Distortion coefficients (radial and tangential) (4x1 vector)
%
%OUTPUT: x: undistorted (normalized) point coordinates in the image plane (2xN matrix)
%
%Method: Iterative method for compensation.
%
%NOTE: This compensation has to be done after the subtraction
%      of the principal point, and division by the focal length.


if length(k) == 1,
    
    [x] = comp_distortion(xd,k);
    
else
    
    k1 = k(1);
    k2 = k(2);
    k3 = k(5);
    p1 = k(3);
    p2 = k(4);
    
    x = xd; 				% initial guess
    
    for kk=1:20,
        
        r_2 = sum(x.^2);
        k_radial =  1 + k1 * r_2 + k2 * r_2.^2 + k3 * r_2.^3;
        delta_x = [2*p1*x(1,:).*x(2,:) + p2*(r_2 + 2*x(1,:).^2);
            p1 * (r_2 + 2*x(2,:).^2)+2*p2*x(1,:).*x(2,:)];
        x = (xd - delta_x)./(ones(2,1)*k_radial);
        
    end;
    
end;

end

function deltaLutNew = MakeOfst2(inValidY,deltaLut,validMat)

imValidYUp = inValidY(1:round(size(inValidY,1)/2),:);
imValidYDown = inValidY(round(size(inValidY,1)/2)+1:end,:);

for ii = 1 : size(imValidYUp,2)
    id = find(imValidYUp(:,ii) ~= 0);
    if ~isempty(id)
        imValidYUp(id,ii) = imValidYUp(id(end),ii);
    end
end
for ii = 1 : size(imValidYDown,2)
    id = find(imValidYDown(:,ii) ~= 0);
    if ~isempty(id)
        imValidYDown(id,ii) = imValidYDown(id(1),ii);
    end
end
inValidY = [imValidYUp; imValidYDown];
deltaLut(~validMat) = inValidY(inValidY~=0);
deltaLutNew = deltaLut;

end
function [LutVecCharX, LutVecCharY,lutVecX,lutVecY,nameMatX, nameMatY] = BiLinearInterp(reverseMapping,xMatSampled, yMatSampled, xOrigMatRecovered, yOrigMatRecovered,imgSize)


optFrac = [5 12 12 9 15 15];
optFracCell = {[1:2];[3:4];[7:8];[9 10 13 14 15 16];[17 18 19 20];[21];[22]};




intLen = [0 0 0 0 13 12 7 6 7 6 7 6 0 0 0 0 0 0 0 0 13 13] +1;


[LutVecCharX,lutVecX,nameMatX] = bilinear2x2(1,reverseMapping, xMatSampled, yMatSampled, xOrigMatRecovered,imgSize,intLen,optFracCell,optFrac);
[LutVecCharY,lutVecY,nameMatY] = bilinear2x2(2,reverseMapping, xMatSampled, yMatSampled, yOrigMatRecovered,imgSize,intLen,optFracCell,optFrac);


end
function [LutVecChar,LutVec,nameMat] = bilinear2x2(coordType,reverseMapping, xMatSampled, yMatSampled, Lut, imgSize,intLen,optFracCell,optFrac)

nr = imgSize(1);
nc = imgSize(2);


optFrac = [optFrac(1:3) max(optFrac(2:3)) optFrac(4:end)];

wordInfo(:,1) = intLen;
for i = 1 : length(optFracCell)
    wordInfo(optFracCell{i},2) = optFrac(i);
end



Lut = round(2^optFrac(1).*Lut)./2^optFrac(1);



if coordType == 1
    Lut2 = xMatSampled - Lut;
else
    Lut2 = yMatSampled - Lut;
end


cropPt1 = [1 1]';
cropPt2 = [nc nr]';

wordInfo = [13 0 0;12 0 0; 11 0 0; 10 0 0; 13 0 0; 12 0 0];
if reverseMapping == 1
    % forward case
    wordInfo = [wordInfo;[5 0 0;6 0 0; 13 1 1;12 1 1; 8 optFrac(1) 1; 8 optFrac(1) 1]];
else
    % reverse case
    wordInfo = [wordInfo;[6 0 0;6 0 0; 13 1 1;12 1 1; 7 optFrac(1) 1; 7 optFrac(1) 1]];
    
end
wordInfo(:,1) = wordInfo(:,1) + wordInfo(:,3);

LutVec = [nc; nr;cropPt1;cropPt2; size(Lut,2);size(Lut,1);  xMatSampled(1,:)'.*2^wordInfo(9,2); yMatSampled(:,1).*2^wordInfo(10,2);Lut2(:).*2^wordInfo(11,2)];
lengthMat = [wordInfo(1:2,:);wordInfo([3 4],:);wordInfo([5 6],:);wordInfo([7 8],:);repmat(wordInfo([9],:),size(xMatSampled,2),1);repmat(wordInfo([10],:),size(yMatSampled,1),1);repmat(wordInfo([11],:),size(yMatSampled,1)*size(yMatSampled,2),1) ];

nameMat = {['ic' ] [sum(wordInfo(1,1:2))];['ir'] [sum(wordInfo(2,1:2))];['x1'] [sum(wordInfo(3,1:2))];['y1'] [sum(wordInfo(4,1:2))];['x2'] [sum(wordInfo(5,1:2))];['y2'] [sum(wordInfo(6,1:2))];['tc'] [sum(wordInfo(7,1:2))];['tr'] [sum(wordInfo(8,1:2))];[repmat('controlX',size(xMatSampled,2),1)] [size(xMatSampled,2)*sum(wordInfo(9,1:2))];[repmat('controlY',size(yMatSampled,1),1)] [size(yMatSampled,1)*sum(wordInfo(10,1:2))];[repmat('Lut',size(yMatSampled,1)*size(yMatSampled,2),1)] [size(yMatSampled,1)*size(yMatSampled,2)*sum(wordInfo(11,1:2))]};


LutVecChar = [];
for u = 1 : length(LutVec)
    if LutVec(u) >= 0
        bin_x = dec2bin(LutVec(u),lengthMat(u,1)+lengthMat(u,2));
        LutVecChar = [LutVecChar '' bin_x];
    else
        bin_x = dec2bin(2^(lengthMat(u,1)+lengthMat(u,2)) + LutVec(u),lengthMat(u,1)+lengthMat(u,2));
        LutVecChar = [LutVecChar '' bin_x];
        
    end
end



end
function makeLut(imgSize, whichCam, para,LutVecOrig2RectXNum,LutVecRect2OrigXNum,LutVecRect2OrigX,LutVecOrig2RectX,Orig2RectNameMatX,Rect2OrigNameMatX,LutVecOrig2RectY,LutVecRect2OrigY,LutVecOrig2RectYNum)
global cfg

supposedLen = 25+21+25+11+12+LutVecOrig2RectXNum(7)*15+LutVecOrig2RectXNum(8)*14+LutVecRect2OrigXNum(7)*15+LutVecRect2OrigXNum(8)*14+LutVecOrig2RectXNum(7)*LutVecOrig2RectXNum(8)*28+LutVecRect2OrigXNum(7)*LutVecRect2OrigXNum(8)*26;
finalVec = [LutVecRect2OrigX(1:25)...
    LutVecOrig2RectX(37:46) LutVecOrig2RectX(26:36)...   %crop1 Y(10) X(11)%                 LutVecOrig2RectX(60:71) LutVecOrig2RectX(47:59)...   %crop2 Y(12) X(13)
    LutVecOrig2RectX(60:71) LutVecRect2OrigX(47:59)...   %crop2 Y(12) X(13)
    LutVecOrig2RectX(72:82) LutVecRect2OrigX(72:83)...
    LutVecOrig2RectX(sum(cell2mat(Orig2RectNameMatX(1:8,2)))+1:sum(cell2mat(Orig2RectNameMatX(1:9,2)))) LutVecOrig2RectY(sum(cell2mat(Orig2RectNameMatX(1:9,2)))+1:sum(cell2mat(Orig2RectNameMatX(1:10,2))))...
    LutVecRect2OrigX(sum(cell2mat(Rect2OrigNameMatX(1:8,2)))+1:sum(cell2mat(Rect2OrigNameMatX(1:9,2)))) LutVecRect2OrigY(sum(cell2mat(Rect2OrigNameMatX(1:9,2)))+1:sum(cell2mat(Rect2OrigNameMatX(1:10,2))))];


data = [{dec2bin(163)};...
    {LutVecRect2OrigX(1:13)};...
    {LutVecOrig2RectX(14:25)};...
    {LutVecRect2OrigX(78:83)};...
    {LutVecRect2OrigX(72:77)};...
    {LutVecOrig2RectX(77:82)};...
    {LutVecOrig2RectX(72:76)};...
    {LutVecOrig2RectX(37:46)};...
    {LutVecOrig2RectX(26:36)};...
    LutVecRect2OrigX(47:59);...
    {LutVecOrig2RectX(60:71)}];



vec1 = [LutVecOrig2RectX(sum(cell2mat(Orig2RectNameMatX(1:8,2)))+1:sum(cell2mat(Orig2RectNameMatX(1:9,2))))];
vec1 = reshape(vec1,15,[]); vec1 = vec1';
for oo = 1 : size(vec1,1)
    data = [data;{vec1(oo,:)}];
    
end
vec2 = [LutVecOrig2RectY(sum(cell2mat(Orig2RectNameMatX(1:9,2)))+1:sum(cell2mat(Orig2RectNameMatX(1:10,2))))];
vec2 = reshape(vec2,14,[]); vec2 = vec2';
for oo = 1 : size(vec2,1)
    data = [data;{vec2(oo,:)}];
    
end
vec3 = [LutVecRect2OrigX(sum(cell2mat(Rect2OrigNameMatX(1:8,2)))+1:sum(cell2mat(Rect2OrigNameMatX(1:9,2))))];
vec3 = reshape(vec3,15,[]); vec3 = vec3';
for oo = 1 : size(vec3,1)
    data = [data;{vec3(oo,:)}];
    
end
vec4 = [LutVecRect2OrigY(sum(cell2mat(Rect2OrigNameMatX(1:9,2)))+1:sum(cell2mat(Rect2OrigNameMatX(1:10,2))))];
vec4 = reshape(vec4,14,[]); vec4 = vec4';
for oo = 1 : size(vec4,1)
    data = [data;{vec4(oo,:)}];
    
end
% forward
for jj = 1 : LutVecOrig2RectYNum(7) * LutVecOrig2RectYNum(8)
    unitLength1 = Orig2RectNameMatX{11,2}/size(Orig2RectNameMatX{11,1},1);
    startId1 = sum(cell2mat(Orig2RectNameMatX(1:10,2)));
    tmp =  [LutVecOrig2RectY(startId1+(jj-1)*unitLength1+1:startId1+jj*unitLength1) LutVecOrig2RectX(startId1+(jj-1)*unitLength1+1:startId1+jj*unitLength1)];
    data = [data;{tmp}];
    
    finalVec = [finalVec tmp];
end
% reverse
for kk = 1 : LutVecRect2OrigXNum(7) * LutVecRect2OrigXNum(8)
    unitLength2 = Rect2OrigNameMatX{11,2}/size(Rect2OrigNameMatX{11,1},1);
    startId2 = sum(cell2mat(Rect2OrigNameMatX(1:10,2)));
    tmp =  [LutVecRect2OrigY(startId2+(kk-1)*unitLength2+1:startId2+kk*unitLength2) LutVecRect2OrigX(startId2+(kk-1)*unitLength2+1:startId2+kk*unitLength2)];
    data = [data;{tmp}];
    
    finalVec = [finalVec tmp];
    
end


errLen = supposedLen - length(finalVec);

assert(errLen == 0);


headerLen = 11;

bodyPart = data(12:end,:);
if ~cfg.isRotate
    fid222 = fopen(fullfile(para,strcat('LutDec',whichCam,sprintf('_%d_%d.txt',imgSize(2), imgSize(1)))),'w');
else
    fid222 = fopen(fullfile(para,strcat('LutDec',whichCam,sprintf('_%d_%d.txt',imgSize(1), imgSize(2)))),'w');
end

binMat = [];

for i = 1 : headerLen + 3
    if i < headerLen + 1
        tempDec = bin2dec(data{i});
        tempBin = dec2bin(tempDec);
        binMat = [binMat; sprintf('%032s',tempBin)];
        fprintf(fid222,sprintf('%d\n',tempDec));
        
    end
    if i == headerLen + 1
        tempHex = '5A5A5A5A';
        tempDec = hex2dec(tempHex);
        tempBin = dec2bin(tempDec);
        binMat = [binMat; sprintf('%032s',tempBin)];
        fprintf(fid222,sprintf('%d\n',tempDec));
    end
    if i == headerLen + 2
        tempHex = dec2hex(length(bodyPart));
        tempDec = hex2dec(tempHex);
        tempBin = dec2bin(tempDec);
        binMat = [binMat; sprintf('%032s',tempBin)];
        fprintf(fid222,sprintf('%d\n',tempDec));
    end
    if i == headerLen + 3
        tempHex = 'FFFFFFFF';
        tempDec = hex2dec(tempHex);
        tempBin = dec2bin(tempDec);
        binMat = [binMat; sprintf('%032s',tempBin)];
        fprintf(fid222,sprintf('%d\n',tempDec));
    end
    
end



for j = 1 : length(bodyPart)
    
    if j < length(bodyPart)
        
        tempBin = bodyPart{j};
        tempDec = bin2dec(tempBin);
        binMat = [binMat; sprintf('%032s',tempBin)];
        fprintf(fid222,sprintf('%d\n',tempDec));
        
        
    else
        
        
        tempBin = bodyPart{j};
        tempDec = bin2dec(tempBin);
        binMat = [binMat; sprintf('%032s',tempBin)];
        fprintf(fid222,sprintf('%d',tempDec));
        
        
        
    end
end

fclose(fid222);




BinMat = binMat';
% BinMat = flipud(BinMat);
% BinMat = [BinMat(17:end,:); BinMat(1:16,:)];

BinMat = [BinMat(9:16,:);BinMat(1:8,:); BinMat(25:32,:);BinMat(17:24,:)];
BinMat = flipud(BinMat);
BinMat = [BinMat(17:end,:); BinMat(1:16,:)];
BinMatVec = str2num(BinMat(:));
if 0
    aaa = (load(fullfile(para,strcat('LutTableDec',whichCam,'.txt'))));
    calibRomFid = fopen(fullfile(inputDir,strcat('LutTableBin',whichCam,'.bin')),'wb');
    %     fwrite(calibRomFid, aaa,'ubit1');
    fwrite(calibRomFid, [aaa],'uint32');
    fclose(calibRomFid);
else
    
    %     aaa = (load(fullfile(inputDir,strcat('LutTableDec',whichOne,'.txt'))));
    if ~cfg.isRotate
        calibRomFid = fopen(fullfile(para,strcat('LutBin',whichCam,sprintf('_%d_%d.bin',imgSize(2), imgSize(1)))),'wb');
    else
        calibRomFid = fopen(fullfile(para,strcat('LutBin',whichCam,sprintf('_%d_%d.bin',imgSize(1), imgSize(2)))),'wb');
    end
    %     fwrite(calibRomFid, aaa,'ubit1');
    fwrite(calibRomFid, [BinMatVec],'ubit1');
    fclose(calibRomFid);
    
end
if ~cfg.isRotate
    calibRomFid = fopen(fullfile(para,strcat('LutBin',whichCam,sprintf('_%d_%d.bin',imgSize(2), imgSize(1)))),'rb');
else
    calibRomFid = fopen(fullfile(para,strcat('LutBin',whichCam,sprintf('_%d_%d.bin',imgSize(1), imgSize(2)))),'rb');
end
qqq = fread(calibRomFid,'uint32');
fclose(calibRomFid);
if 0
    figure,plot(aaa-qqq)
end


end



