function [imgL, imgR, rectImgL, rectImgR, meanError, KK_newL, imgLRGB, imgRRGB] = CallArucoDetection(inputDir, calibFuncDir, paraDir, config)
global cfg inputDir0 stepDir paraDir0


funcDir = calibFuncDir;
funcStrSeg = strsplit(funcDir,'\');

funcName = [];
for j = 1 : length(funcStrSeg)
    funcName = [funcName funcStrSeg{j},'/'];
end

funcName = funcName(1:end-1);

imgLRGB = [];
imgRRGB = [];

KK_newL = [];

valid_boarder_margin = cfg.valid_boarder_margin;

draw1 = 0;
draw2 = 0;

if 1
    % 生成标定板控制点
    markerIdList = 0:cfg.aruco_col*cfg.aruco_row-1;
    sameCornerStep = cfg.aruco_size*4/3;
    [xMat1, yMat1] = meshgrid(0:sameCornerStep:(cfg.aruco_col-1)*sameCornerStep, 0:sameCornerStep:(cfg.aruco_row-1)*sameCornerStep);
    [xMat2, yMat2] = meshgrid(cfg.aruco_size:sameCornerStep:(cfg.aruco_col-1)*sameCornerStep+cfg.aruco_size, 0:sameCornerStep:(cfg.aruco_row-1)*sameCornerStep);
    [xMat3, yMat3] = meshgrid(cfg.aruco_size:sameCornerStep:(cfg.aruco_col-1)*sameCornerStep+cfg.aruco_size, cfg.aruco_size:sameCornerStep:(cfg.aruco_row-1)*sameCornerStep+cfg.aruco_size);
    [xMat4, yMat4] = meshgrid(0:sameCornerStep:(cfg.aruco_col-1)*sameCornerStep+0, cfg.aruco_size:sameCornerStep:(cfg.aruco_row-1)*sameCornerStep+cfg.aruco_size);
    
    xMat1_t = xMat1';yMat1_t = yMat1';
    xMat2_t = xMat2';yMat2_t = yMat2';
    xMat3_t = xMat3';yMat3_t = yMat3';
    xMat4_t = xMat4';yMat4_t = yMat4';
    
    xyz1 = [xMat1_t(:) yMat1_t(:) zeros(cfg.aruco_col*cfg.aruco_row, 1)];
    xyz2 = [xMat2_t(:) yMat2_t(:) zeros(cfg.aruco_col*cfg.aruco_row, 1)];
    xyz3 = [xMat3_t(:) yMat3_t(:) zeros(cfg.aruco_col*cfg.aruco_row, 1)];
    xyz4 = [xMat4_t(:) yMat4_t(:) zeros(cfg.aruco_col*cfg.aruco_row, 1)];
    
    xyzCellGlobal = [{xyz1};{xyz2};{xyz3};{xyz4}];
    
end
if draw1
    figure,imshow(0.5.*ones(250,350)); hold on;
    for j = 1 : length(markerIdList)
        text(xyz1(j,1)+3,xyz1(j,2)+3,num2str(markerIdList(j)), 'Color',[1 0 1],'FontSize',15);
        plot(xyz1(:,1), xyz1(:,2),'-or');
        plot(xyz2(:,1), xyz2(:,2),'.g');
        plot(xyz3(:,1), xyz3(:,2),'.b');
        plot(xyz4(:,1), xyz4(:,2),'.y');
    end
    
end

if ~cfg.is_kepler
    
    strSeg = strsplit(inputDir,'\');
    
    inputName = [];
    for j = 1 : length(strSeg)
        inputName = [inputName strSeg{j},'/'];
        
    end
    inputName = inputName(1:end-1);
    
    cmdStr = [strcat(funcName,'/ConsoleApplication4.exe') ' ',num2str(cfg.aruco_col),' ',num2str(cfg.aruco_row),' ','1', ' ',inputName];
    % 调用opencv检测aruco角点
    system(cmdStr);
    
    dirInfo = dir(fullfile(inputDir,'img*.png'));
    frameNum = length(dirInfo)/2;
    cornerInfo = dir(fullfile(inputDir,'corners*.txt'));
    leftCorner = [];
    rightCorner = [];
    goodIdL = [];
    goodIdR = [];
    for i = 1 : size(cornerInfo,1)
        cornerName = cornerInfo(i).name;
        corner_ = load(fullfile(inputDir, cornerInfo(i).name));
        [~,ind] = sort(corner_(:,1));
        corner = corner_(ind,:);
        if 1
            corner(:,2:9) = corner(:,2:9)+1;
        else
            %% 
            corner(:,2:9) = corner(:,2:9)+0;
        end
        pixCell00 = [{corner(:,2:3)};{corner(:,4:5)};{corner(:,6:7)};{corner(:,8:9)}];
        corner00 = corner;
        
        % 过滤靠近图像边缘的角点
        valid_corner1 = corner(:,2) > valid_boarder_margin & corner(:,2) < cfg.img_width - valid_boarder_margin & corner(:,3) > valid_boarder_margin & corner(:,3) < cfg.img_height - valid_boarder_margin ;
        valid_corner2 = corner(:,4) > valid_boarder_margin & corner(:,4) < cfg.img_width - valid_boarder_margin & corner(:,5) > valid_boarder_margin & corner(:,5) < cfg.img_height - valid_boarder_margin ;
        valid_corner3 = corner(:,6) > valid_boarder_margin & corner(:,6) < cfg.img_width - valid_boarder_margin & corner(:,7) > valid_boarder_margin & corner(:,7) < cfg.img_height - valid_boarder_margin ;
        valid_corner4 = corner(:,8) > valid_boarder_margin & corner(:,8) < cfg.img_width - valid_boarder_margin & corner(:,9) > valid_boarder_margin & corner(:,9) < cfg.img_height - valid_boarder_margin ;
        valid_corner = valid_corner1 & valid_corner2 & valid_corner3 & valid_corner4;
        corner = corner00(valid_corner,:);
        corner0 = corner;
        
        pixCell = [{corner(:,2:3)};{corner(:,4:5)};{corner(:,6:7)};{corner(:,8:9)}];
        pixCell0 = pixCell;
        
        
        [~,distL1] = NormalizeVector(pixCell{1} - pixCell{2});
        [~,distL2] = NormalizeVector(pixCell{2} - pixCell{3});
        [~,distL3] = NormalizeVector(pixCell{3} - pixCell{4});
        [~,distL4] = NormalizeVector(pixCell{4} - pixCell{1});
        [~,distL5] = NormalizeVector(pixCell{1} - pixCell{3});
        [~,distL6] = NormalizeVector(pixCell{2} - pixCell{4});
        distMatL = [distL1 distL2 distL3 distL4 distL5 distL6];
        distMatL_min = min(distMatL')';
        % 过滤尺寸太小的aruco id
        validMarker = find(distMatL_min > cfg.min_aruco_len);
        corner = corner0(validMarker,:);
        pixCell = [{corner(:,2:3)};{corner(:,4:5)};{corner(:,6:7)};{corner(:,8:9)}];
        
        % 仅当图片中有效的aruco数目大于某阈值时，才认为这张图片有效，否则跳过。具体角点检测情况可通过令“draw2 = 1”查看。
        if size(corner,1) > cfg.min_valid_aruco
            pix = cell2mat(pixCell);
            if strcmp(cornerName(8),'L')
                goodIdL = [goodIdL;str2double(cornerName(end-7:end-4))];
                imgL = imread(fullfile(inputDir, dirInfo(goodIdL(end)).name));
                if 0
                    pixList = cell2mat(pixCell);
                    pixList_ = cornerfinder([pixList]',double(rgb2gray(imgL)),5,5);
                end
                if draw2
                    figure(10),imshow(imgL);hold on;
                    for k = 1 : size(corner00, 1)
                        text(pixCell00{1}(k,1)+3,pixCell00{1}(k,2)+3,num2str(corner00(k,1)), 'Color',[1 1 1],'FontSize',15);  %,'FontWeight','bold');
                        plot(pixCell00{1}(:,1), pixCell00{1}(:,2),'or');
                        plot(pixCell00{2}(:,1), pixCell00{2}(:,2),'og');
                        plot(pixCell00{3}(:,1), pixCell00{3}(:,2),'ob');
                        plot(pixCell00{4}(:,1), pixCell00{4}(:,2),'oy');
                    end
                    for k = 1 : size(corner, 1)
                        text(pixCell{1}(k,1)+3,pixCell{1}(k,2)+3,num2str(corner(k,1)), 'Color',[1 0 1],'FontSize',15);  %,'FontWeight','bold');
                        plot(pixCell{1}(:,1), pixCell{1}(:,2),'xr');
                        plot(pixCell{2}(:,1), pixCell{2}(:,2),'xg');
                        plot(pixCell{3}(:,1), pixCell{3}(:,2),'xb');
                        plot(pixCell{4}(:,1), pixCell{4}(:,2),'xy');
                    end
                end
                %            [~,distL1] = NormalizeVector(pixCell{1} - pixCell{2});
                %            [~,distL2] = NormalizeVector(pixCell{2} - pixCell{3});
                %            [~,distL3] = NormalizeVector(pixCell{3} - pixCell{4});
                %            [~,distL4] = NormalizeVector(pixCell{4} - pixCell{1});
                %            [~,distL5] = NormalizeVector(pixCell{1} - pixCell{3});
                %            [~,distL6] = NormalizeVector(pixCell{2} - pixCell{4});
                %            distMatL = [distL1 distL2 distL3 distL4 distL5 distL6];
                %            distMatL_min = min(distMatL')';
                leftCorner = [leftCorner;{[{corner(:,1)} {pixCell}]}];
            end
            if strcmp(cornerName(8),'R')
                goodIdR = [goodIdR;str2double(cornerName(end-7:end-4))];
                imgR = imread(fullfile(inputDir, dirInfo(goodIdR(end)+frameNum).name));
                
                if draw2
                    figure(10),imshow(imgR);hold on;
                    for k = 1 : size(corner00, 1)
                        text(pixCell00{1}(k,1)+3,pixCell00{1}(k,2)+3,num2str(corner00(k,1)), 'Color',[1 1 1],'FontSize',15);  %,'FontWeight','bold');
                        plot(pixCell00{1}(:,1), pixCell00{1}(:,2),'or');
                        plot(pixCell00{2}(:,1), pixCell00{2}(:,2),'og');
                        plot(pixCell00{3}(:,1), pixCell00{3}(:,2),'ob');
                        plot(pixCell00{4}(:,1), pixCell00{4}(:,2),'oy');
                    end
                    for k = 1 : size(corner, 1)
                        text(pixCell{1}(k,1)+3,pixCell{1}(k,2)+3,num2str(corner(k,1)), 'Color',[1 0 1],'FontSize',15);  %,'FontWeight','bold');
                        plot(pixCell{1}(:,1), pixCell{1}(:,2),'xr');
                        plot(pixCell{2}(:,1), pixCell{2}(:,2),'xg');
                        plot(pixCell{3}(:,1), pixCell{3}(:,2),'xb');
                        plot(pixCell{4}(:,1), pixCell{4}(:,2),'xy');
                    end
                end
                
                rightCorner = [rightCorner;{[{corner(:,1)} {pixCell}]}];
            end
        end
    end
    
    leftCorner0 = leftCorner;
    rightCorner0 = rightCorner;
    % 找到左右有效图片中公共的id
    commId = intersect(goodIdL, goodIdR);
    
    leftCorner = leftCorner0(ismember(goodIdL, commId));
    rightCorner = rightCorner0(ismember(goodIdR, commId));
    
    imgListL = {}; imgListR = {};
    leftTrack = {};
    rightTrack = {};
    % 找到与图片中每个2d点对应的3d点
    for i = 1 : length(commId)
        left_corner_id = leftCorner{i}{1};
        right_corner_id = rightCorner{i}{1};
        comm_id = intersect(left_corner_id, right_corner_id);
        
        leftCorner_temp0 = leftCorner{i}{2};
        rightCorner_temp0 = rightCorner{i}{2};
        
        globalInd = ismember(markerIdList', comm_id);
        
        leftCorner_temp = {}; rightCorner_temp = {};
        for jj = 1 : 4
            leftCorner_temp{jj,1} = [leftCorner_temp0{jj}(ismember(left_corner_id, comm_id),:) xyzCellGlobal{jj}(globalInd,:)];
            rightCorner_temp{jj,1} = [rightCorner_temp0{jj}(ismember(right_corner_id, comm_id),:) xyzCellGlobal{jj}(globalInd,:)];
        end
        if ~cfg.switch_lr
            leftTrack = [leftTrack; {[{comm_id} {leftCorner_temp}]}];
            rightTrack = [rightTrack; {[{comm_id} {rightCorner_temp}]}];
            
            imgListL = [imgListL; fullfile(inputDir, dirInfo(commId(i)).name)];
            imgListR = [imgListR; fullfile(inputDir, dirInfo(commId(i) + frameNum).name)];
        else
            rightTrack = [rightTrack; {[{comm_id} {leftCorner_temp}]}];
            leftTrack = [leftTrack; {[{comm_id} {rightCorner_temp}]}];
            
            imgListR = [imgListR; fullfile(inputDir, dirInfo(commId(i)).name)];
            imgListL = [imgListL; fullfile(inputDir, dirInfo(commId(i) + frameNum).name)];
            
        end
    end
    
    imgSize = size(imgL(:,:,1));
    leftTrack0 = leftTrack;
    rightTrack0 = rightTrack;
    % 分别单独标定左相机和右相机
    [camParamL0, cbcXYL0, cbGridL0, configL0] = CalibArucoSingle(imgListL, leftTrack0, imgSize, config);
    [camParamR0, cbcXYR0, cbGridR0, configR0] = CalibArucoSingle(imgListR, rightTrack0, imgSize, config);
    if cfg.reproj_pix_thr < 20
        % 做一次deoutlier，剔除掉重投影大的点后重新标定左右相机
        [leftTrack, rightTrack] = Deoutlier(leftTrack0, camParamL0, cbcXYL0, cbGridL0, rightTrack0, camParamR0, cbcXYR0, cbGridR0);
        
        [camParamL, cbcXYL, cbGridL, configL] = CalibArucoSingle(imgListL, leftTrack, imgSize, config);
        [camParamR, cbcXYR, cbGridR, configR] = CalibArucoSingle(imgListR, rightTrack, imgSize, config);
    else
        % 实验发现，deoutlier对精度提升不明显，为了提高速度考虑，当“config.txt”中配置的“reproj_pix_thr”大于20pixel时，直接跳过deoutlier这一步
        camParamL = camParamL0;
        cbcXYL = cbcXYL0;
        cbGridL = cbGridL0;
        configL = configL0;
        
        camParamR = camParamR0;
        cbcXYR = cbcXYR0;
        cbGridR = cbGridR0;
        configR = configR0;
    end
    % 标定左右相机外参，同时优化各自内参
    stereoParam = CbCalibStereo(camParamL, camParamR, cbcXYL, cbcXYR, cbGridL, configL, configR);
    
    
    
    
    % 量化双目标定精度，计算左目到右目的重投影误差
    [KL, KR, rectParamL, rectParamR, rotMatLeft, rotMatRight,  intrMatLeftNew, intrMatRightNew, meanError, KK_newL] = CheckStereoCalibReprojection(stereoParam, imgSize, imgListL, imgListR, camParamL, camParamR, cbcXYL, cbcXYR, cbGridL, cbGridR, markerIdList, xyzCellGlobal);
    intrMatNewL = intrMatLeftNew;
    intrMatNewR = intrMatRightNew;
    rotMatL = rotMatLeft;
    rotMatR = rotMatRight;
    save(fullfile(inputDir0, 'calib.mat'), 'stereoParam', 'intrMatNewL', 'intrMatNewR', 'rotMatL', 'rotMatR');
    % 判断相机是否装反，若反了则在report.txt中写入0，后退出程序。
    if ~cfg.isRotate
        if stereoParam.transVecRef(1) > 0
            fid1 = fopen(fullfile(inputDir,'report.txt'),'w');%????
            fprintf(fid1,sprintf('%d\n',0));
            fclose(fid1);
            copyfile(fullfile(paraDir0,'config.txt'), fullfile(inputDir));
            exit(0)
        end
    else
        if stereoParam.transVecRef(2) > 0
            fid1 = fopen(fullfile(inputDir,'report.txt'),'w');%????
            fprintf(fid1,sprintf('%d\n',0));
            fclose(fid1);
            copyfile(fullfile(paraDir0,'config.txt'), fullfile(inputDir));
            exit(0)
        end
    end
    
    
    cfg.lr_param.rectParamL = rectParamL;
    cfg.lr_param.rectParamR = rectParamR;
    cfg.lr_param.rotMatLeft = rotMatLeft;
    cfg.lr_param.rotMatRight = rotMatRight;
    cfg.lr_param.intrMatLeftNew = intrMatLeftNew;
    cfg.lr_param.intrMatRightNew = intrMatRightNew;
    
    imgL = imread(imgListL{1});
    imgR = imread(imgListR{1});
    
    [matX, matY] = meshgrid(1:imgSize(2), 1:imgSize(1));
    pixRectL = remapRect([matX(:) matY(:)]', intrMatLeftNew, KL, stereoParam.kcLeft, rotMatLeft);
    pixRectR = remapRect([matX(:) matY(:)]', intrMatRightNew, KR, stereoParam.kcRight, rotMatRight);
    
    reMapLX = reshape(pixRectL(:,1), imgSize);
    reMapLY = reshape(pixRectL(:,2), imgSize);
    reMapRX = reshape(pixRectR(:,1), imgSize);
    reMapRY = reshape(pixRectR(:,2), imgSize);
    
    % 测试rectify效果
    if 1 
        rectImgL = uint8(interp2(matX, matY, double(imgL(:,:,1)),reMapLX,reMapLY));
        rectImgR = uint8(interp2(matX, matY, double(imgR(:,:,1)),reMapRX,reMapRY));
    else
        rectImgL = uint8(interp2(matX, matY, double(imgR(:,:,1)),reMapLX,reMapLY));
        rectImgR = uint8(interp2(matX, matY, double(imgL(:,:,1)),reMapRX,reMapRY));
    end
    
    save(fullfile(inputDir0, 'rectMap.mat'), 'matX', 'matY', 'reMapLX', 'reMapLY', 'reMapRX', 'reMapRY');
 
    cfg.imgL = imgL;
    cfg.imgR = imgR;
    % 生成芯片吃的参数
    ConvertParam2(inputDir0, stepDir, stereoParam, imgSize, 1);
    
    if 0
        [rectParamL, rectParamR, rotMatLeft, rotMatRight] = GetRectifyParam2(stereoParam, imgSize);
        [rectImgL, rectImgR] = RectifyImagePair(stereoParam, imgL, imgR);
    end
    
else
%     funcDir = calibFuncDir;
%     funcStrSeg = strsplit(funcDir,'\');
%     
%     funcName = [];
%     for j = 1 : length(funcStrSeg)
%         funcName = [funcName funcStrSeg{j},'/'];
%     end
    
    
    strSeg = strsplit(inputDir,'\');
    inputName = [];
    for j = 1 : length(strSeg)
        inputName = [inputName strSeg{j},'/'];
        
    end
    inputName = inputName(1:end-1);
    
    cmdStr1 = [strcat(funcName,'/ConsoleApplication4.exe') ' ',num2str(cfg.aruco_col),' ',num2str(cfg.aruco_row),' ','1', ' ',inputName];
    system(cmdStr1);
    cmdStr2 = [strcat(funcName,'/ConsoleApplication4.exe') ' ',num2str(cfg.aruco_col),' ',num2str(cfg.aruco_row),' ','0', ' ',strcat(inputName,'/RGB')];
    system(cmdStr2);
    
    
    dirInfo = dir(fullfile(inputDir,'img*.png'));
    dirInfoRGB = dir(fullfile(inputDir,'RGB','img*.png'));
    frameNum = length(dirInfo)/2;
    cornerInfo = dir(fullfile(inputDir,'corners*.txt'));
    cornerInfoRGB = dir(fullfile(fullfile(inputDir,'RGB'),'corner_*.txt'));
    leftCorner = [];
    rightCorner = [];
    rgbCorner = [];
    goodIdL = [];
    goodIdR = [];
    goodIdRGB = [];
    for i = 1 : size(cornerInfo,1)
        
        if i <= min(frameNum, size(cornerInfoRGB,1))
            cornerNameRGB = cornerInfoRGB(i).name;
            cornerRGB_ = load(fullfile(inputDir,'RGB', cornerInfoRGB(i).name));
            [~,indRGB] = sort(cornerRGB_(:,1));
            cornerRGB = cornerRGB_(indRGB,:);
            cornerRGB(:,2:9) = cornerRGB(:,2:9)+1;
            pixCellRGB00 = [{cornerRGB(:,2:3)};{cornerRGB(:,4:5)};{cornerRGB(:,6:7)};{cornerRGB(:,8:9)}];
            cornerRGB00 = cornerRGB;
            
            valid_cornerRGB1 = cornerRGB(:,2) > valid_boarder_margin & cornerRGB(:,2) < cfg.img_width_rgb - valid_boarder_margin & cornerRGB(:,3) > valid_boarder_margin & cornerRGB(:,3) < cfg.img_height_rgb - valid_boarder_margin ;
            valid_cornerRGB2 = cornerRGB(:,4) > valid_boarder_margin & cornerRGB(:,4) < cfg.img_width_rgb - valid_boarder_margin & cornerRGB(:,5) > valid_boarder_margin & cornerRGB(:,5) < cfg.img_height_rgb - valid_boarder_margin ;
            valid_cornerRGB3 = cornerRGB(:,6) > valid_boarder_margin & cornerRGB(:,6) < cfg.img_width_rgb - valid_boarder_margin & cornerRGB(:,7) > valid_boarder_margin & cornerRGB(:,7) < cfg.img_height_rgb - valid_boarder_margin ;
            valid_cornerRGB4 = cornerRGB(:,8) > valid_boarder_margin & cornerRGB(:,8) < cfg.img_width_rgb - valid_boarder_margin & cornerRGB(:,9) > valid_boarder_margin & cornerRGB(:,9) < cfg.img_height_rgb - valid_boarder_margin ;
            valid_cornerRGB = valid_cornerRGB1 & valid_cornerRGB2 & valid_cornerRGB3 & valid_cornerRGB4;
            cornerRGB = cornerRGB00(valid_cornerRGB,:);
            cornerRGB0 = cornerRGB;
            
            pixCellRGB = [{cornerRGB(:,2:3)};{cornerRGB(:,4:5)};{cornerRGB(:,6:7)};{cornerRGB(:,8:9)}];
            pixCellRGB0 = pixCellRGB;
            
            
            [~,distLRGB1] = NormalizeVector(pixCellRGB{1} - pixCellRGB{2});
            [~,distLRGB2] = NormalizeVector(pixCellRGB{2} - pixCellRGB{3});
            [~,distLRGB3] = NormalizeVector(pixCellRGB{3} - pixCellRGB{4});
            [~,distLRGB4] = NormalizeVector(pixCellRGB{4} - pixCellRGB{1});
            [~,distLRGB5] = NormalizeVector(pixCellRGB{1} - pixCellRGB{3});
            [~,distLRGB6] = NormalizeVector(pixCellRGB{2} - pixCellRGB{4});
            distMatLRGB = [distLRGB1 distLRGB2 distLRGB3 distLRGB4 distLRGB5 distLRGB6];
            distMatLRGB_min = min(distMatLRGB')';
            validMarkerRGB = find(distMatLRGB_min > cfg.min_aruco_len);
            cornerRGB = cornerRGB0(validMarkerRGB,:);
            pixCellRGB = [{cornerRGB(:,2:3)};{cornerRGB(:,4:5)};{cornerRGB(:,6:7)};{cornerRGB(:,8:9)}];
            
            if size(cornerRGB,1) > cfg.min_valid_aruco
                goodIdRGB = [goodIdRGB;str2double(cornerNameRGB(end-7:end-4))];
                imgRGB = imread(fullfile(inputDir, 'RGB', dirInfoRGB(goodIdRGB(end)).name));
                if draw2
                    figure(10),imshow(imgRGB);hold on;
                    for k = 1 : size(cornerRGB00, 1)
                        text(pixCellRGB00{1}(k,1)+3,pixCellRGB00{1}(k,2)+3,num2str(cornerRGB00(k,1)), 'Color',[1 1 1],'FontSize',15);  %,'FontWeight','bold');
                        plot(pixCellRGB00{1}(:,1), pixCellRGB00{1}(:,2),'or');
                        plot(pixCellRGB00{2}(:,1), pixCellRGB00{2}(:,2),'og');
                        plot(pixCellRGB00{3}(:,1), pixCellRGB00{3}(:,2),'ob');
                        plot(pixCellRGB00{4}(:,1), pixCellRGB00{4}(:,2),'oy');
                    end
                    for k = 1 : size(cornerRGB, 1)
                        text(pixCellRGB{1}(k,1)+3,pixCellRGB{1}(k,2)+3,num2str(cornerRGB(k,1)), 'Color',[1 0 1],'FontSize',15);  %,'FontWeight','bold');
                        plot(pixCellRGB{1}(:,1), pixCellRGB{1}(:,2),'xr');
                        plot(pixCellRGB{2}(:,1), pixCellRGB{2}(:,2),'xg');
                        plot(pixCellRGB{3}(:,1), pixCellRGB{3}(:,2),'xb');
                        plot(pixCellRGB{4}(:,1), pixCellRGB{4}(:,2),'xy');
                    end
                end
                rgbCorner = [rgbCorner;{[{cornerRGB(:,1)} {pixCellRGB}]}];
            
            end
            
        end
        
        
        
        cornerName = cornerInfo(i).name;
        corner_ = load(fullfile(inputDir, cornerInfo(i).name));
        [~,ind] = sort(corner_(:,1));
        corner = corner_(ind,:);
        corner(:,2:9) = corner(:,2:9)+1;
        pixCell00 = [{corner(:,2:3)};{corner(:,4:5)};{corner(:,6:7)};{corner(:,8:9)}];
        corner00 = corner;
        
        valid_corner1 = corner(:,2) > valid_boarder_margin & corner(:,2) < cfg.img_width - valid_boarder_margin & corner(:,3) > valid_boarder_margin & corner(:,3) < cfg.img_height - valid_boarder_margin ;
        valid_corner2 = corner(:,4) > valid_boarder_margin & corner(:,4) < cfg.img_width - valid_boarder_margin & corner(:,5) > valid_boarder_margin & corner(:,5) < cfg.img_height - valid_boarder_margin ;
        valid_corner3 = corner(:,6) > valid_boarder_margin & corner(:,6) < cfg.img_width - valid_boarder_margin & corner(:,7) > valid_boarder_margin & corner(:,7) < cfg.img_height - valid_boarder_margin ;
        valid_corner4 = corner(:,8) > valid_boarder_margin & corner(:,8) < cfg.img_width - valid_boarder_margin & corner(:,9) > valid_boarder_margin & corner(:,9) < cfg.img_height - valid_boarder_margin ;
        valid_corner = valid_corner1 & valid_corner2 & valid_corner3 & valid_corner4;
        corner = corner00(valid_corner,:);
        corner0 = corner;
        
        pixCell = [{corner(:,2:3)};{corner(:,4:5)};{corner(:,6:7)};{corner(:,8:9)}];
        pixCell0 = pixCell;
        
        
        [~,distL1] = NormalizeVector(pixCell{1} - pixCell{2});
        [~,distL2] = NormalizeVector(pixCell{2} - pixCell{3});
        [~,distL3] = NormalizeVector(pixCell{3} - pixCell{4});
        [~,distL4] = NormalizeVector(pixCell{4} - pixCell{1});
        [~,distL5] = NormalizeVector(pixCell{1} - pixCell{3});
        [~,distL6] = NormalizeVector(pixCell{2} - pixCell{4});
        distMatL = [distL1 distL2 distL3 distL4 distL5 distL6];
        distMatL_min = min(distMatL')';
        validMarker = find(distMatL_min > cfg.min_aruco_len);
        corner = corner0(validMarker,:);
        pixCell = [{corner(:,2:3)};{corner(:,4:5)};{corner(:,6:7)};{corner(:,8:9)}];
        
        
        if size(corner,1) > cfg.min_valid_aruco
            pix = cell2mat(pixCell);
            if strcmp(cornerName(8),'L')
                goodIdL = [goodIdL;str2double(cornerName(end-7:end-4))];
                imgL = imread(fullfile(inputDir, dirInfo(goodIdL(end)).name));
                
                if draw2
                    figure(10),imshow(imgL);hold on;
                    for k = 1 : size(corner00, 1)
                        text(pixCell00{1}(k,1)+3,pixCell00{1}(k,2)+3,num2str(corner00(k,1)), 'Color',[1 1 1],'FontSize',15);  %,'FontWeight','bold');
                        plot(pixCell00{1}(:,1), pixCell00{1}(:,2),'or');
                        plot(pixCell00{2}(:,1), pixCell00{2}(:,2),'og');
                        plot(pixCell00{3}(:,1), pixCell00{3}(:,2),'ob');
                        plot(pixCell00{4}(:,1), pixCell00{4}(:,2),'oy');
                    end
                    for k = 1 : size(corner, 1)
                        text(pixCell{1}(k,1)+3,pixCell{1}(k,2)+3,num2str(corner(k,1)), 'Color',[1 0 1],'FontSize',15);  %,'FontWeight','bold');
                        plot(pixCell{1}(:,1), pixCell{1}(:,2),'xr');
                        plot(pixCell{2}(:,1), pixCell{2}(:,2),'xg');
                        plot(pixCell{3}(:,1), pixCell{3}(:,2),'xb');
                        plot(pixCell{4}(:,1), pixCell{4}(:,2),'xy');
                    end
                end
                %            [~,distL1] = NormalizeVector(pixCell{1} - pixCell{2});
                %            [~,distL2] = NormalizeVector(pixCell{2} - pixCell{3});
                %            [~,distL3] = NormalizeVector(pixCell{3} - pixCell{4});
                %            [~,distL4] = NormalizeVector(pixCell{4} - pixCell{1});
                %            [~,distL5] = NormalizeVector(pixCell{1} - pixCell{3});
                %            [~,distL6] = NormalizeVector(pixCell{2} - pixCell{4});
                %            distMatL = [distL1 distL2 distL3 distL4 distL5 distL6];
                %            distMatL_min = min(distMatL')';
                leftCorner = [leftCorner;{[{corner(:,1)} {pixCell}]}];
            end
            if strcmp(cornerName(8),'R')
                goodIdR = [goodIdR;str2double(cornerName(end-7:end-4))];
                imgR = imread(fullfile(inputDir, dirInfo(goodIdR(end)+frameNum).name));
                
                if draw2
                    figure(10),imshow(imgR);hold on;
                    for k = 1 : size(corner00, 1)
                        text(pixCell00{1}(k,1)+3,pixCell00{1}(k,2)+3,num2str(corner00(k,1)), 'Color',[1 1 1],'FontSize',15);  %,'FontWeight','bold');
                        plot(pixCell00{1}(:,1), pixCell00{1}(:,2),'or');
                        plot(pixCell00{2}(:,1), pixCell00{2}(:,2),'og');
                        plot(pixCell00{3}(:,1), pixCell00{3}(:,2),'ob');
                        plot(pixCell00{4}(:,1), pixCell00{4}(:,2),'oy');
                    end
                    for k = 1 : size(corner, 1)
                        text(pixCell{1}(k,1)+3,pixCell{1}(k,2)+3,num2str(corner(k,1)), 'Color',[1 0 1],'FontSize',15);  %,'FontWeight','bold');
                        plot(pixCell{1}(:,1), pixCell{1}(:,2),'xr');
                        plot(pixCell{2}(:,1), pixCell{2}(:,2),'xg');
                        plot(pixCell{3}(:,1), pixCell{3}(:,2),'xb');
                        plot(pixCell{4}(:,1), pixCell{4}(:,2),'xy');
                    end
                end
                
                rightCorner = [rightCorner;{[{corner(:,1)} {pixCell}]}];
            end
        end
    end
    
    leftCorner0 = leftCorner;
    rightCorner0 = rightCorner;
    rgbCorner0 = rgbCorner;
    
    commId_ = intersect(goodIdL, goodIdR);
    commId = intersect(commId_, goodIdRGB);
    
    leftCorner = leftCorner0(ismember(goodIdL, commId));
    rightCorner = rightCorner0(ismember(goodIdR, commId));
    rgbCorner = rgbCorner0(ismember(goodIdRGB, commId));
    
    imgListL = {}; imgListR = {}; imgListRGB = {};
    leftTrack = {};
    rightTrack = {};
    rgbTrack = {};
    for i = 1 : length(commId)
        left_corner_id = leftCorner{i}{1};
        right_corner_id = rightCorner{i}{1};
        rgb_corner_id = rgbCorner{i}{1};
        comm_id_ = intersect(left_corner_id, right_corner_id);
        comm_id = intersect(comm_id_, rgb_corner_id);
        
        leftCorner_temp0 = leftCorner{i}{2};
        rightCorner_temp0 = rightCorner{i}{2};
        rgbCorner_temp0 = rgbCorner{i}{2};
        
        globalInd = ismember(markerIdList', comm_id);
        
        leftCorner_temp = {}; rightCorner_temp = {}; rgbCorner_temp = {};
        for jj = 1 : 4
            leftCorner_temp{jj,1} = [leftCorner_temp0{jj}(ismember(left_corner_id, comm_id),:) xyzCellGlobal{jj}(globalInd,:)];
            rightCorner_temp{jj,1} = [rightCorner_temp0{jj}(ismember(right_corner_id, comm_id),:) xyzCellGlobal{jj}(globalInd,:)];
            rgbCorner_temp{jj,1} = [rgbCorner_temp0{jj}(ismember(rgb_corner_id, comm_id),:) xyzCellGlobal{jj}(globalInd,:)];
        end
        rgbTrack = [rgbTrack; {[{comm_id} {rgbCorner_temp}]}];
        imgListRGB = [imgListRGB; fullfile(inputDir, 'RGB',dirInfoRGB(commId(i)).name)];
        if ~cfg.switch_lr
            leftTrack = [leftTrack; {[{comm_id} {leftCorner_temp}]}];
            rightTrack = [rightTrack; {[{comm_id} {rightCorner_temp}]}];
            
            imgListL = [imgListL; fullfile(inputDir, dirInfo(commId(i)).name)];
            imgListR = [imgListR; fullfile(inputDir, dirInfo(commId(i) + frameNum).name)];
        else
            rightTrack = [rightTrack; {[{comm_id} {leftCorner_temp}]}];
            leftTrack = [leftTrack; {[{comm_id} {rightCorner_temp}]}];
            
            imgListR = [imgListR; fullfile(inputDir, dirInfo(commId(i)).name)];
            imgListL = [imgListL; fullfile(inputDir, dirInfo(commId(i) + frameNum).name)];
            
        end
    end
    
    imgSize = size(imgL(:,:,1));
    imgSizeRGB = size(imgRGB(:,:,1));
    leftTrack0 = leftTrack;
    rightTrack0 = rightTrack;
    rgbTrack0 = rgbTrack;
    
    [camParamL0, cbcXYL0, cbGridL0, configL0] = CalibArucoSingle(imgListL, leftTrack0, imgSize, config);
    [camParamR0, cbcXYR0, cbGridR0, configR0] = CalibArucoSingle(imgListR, rightTrack0, imgSize, config);
    [camParamRGB0, cbcXYRGB0, cbGridRGB0, configRGB0] = CalibArucoSingle(imgListRGB, rgbTrack0, imgSizeRGB, config);
    
    if cfg.reproj_pix_thr < 20
        [leftTrack, rightTrack, rgbTrack] = Deoutlier(leftTrack0, camParamL0, cbcXYL0, cbGridL0, rightTrack0, camParamR0, cbcXYR0, cbGridR0, rgbTrack0, camParamRGB0, cbcXYRGB0, cbGridRGB0);
        
        [camParamL, cbcXYL, cbGridL, configL] = CalibArucoSingle(imgListL, leftTrack, imgSize, config);
        [camParamR, cbcXYR, cbGridR, configR] = CalibArucoSingle(imgListR, rightTrack, imgSize, config);
        [camParamRGB, cbcXYRGB, cbGridRGB, configRGB] = CalibArucoSingle(imgListRGB, rgbTrack, imgSizeRGB, config);
    else
        camParamL = camParamL0;
        cbcXYL = cbcXYL0;
        cbGridL = cbGridL0;
        configL = configL0;
        
        camParamR = camParamR0;
        cbcXYR = cbcXYR0;
        cbGridR = cbGridR0;
        configR = configR0;
        
        camParamRGB = camParamRGB0;
        cbcXYRGB = cbcXYRGB0;
        cbGridRGB = cbGridRGB0;
        configRGB = configRGB0;
        
    end
    stereoParam = CbCalibStereo(camParamL, camParamR, cbcXYL, cbcXYR, cbGridL, configL, configR);
    
     if ~cfg.isRotate
        if stereoParam.transVecRef(1) > 0
            fid1 = fopen(fullfile(inputDir,'report.txt'),'w');%????
            fprintf(fid1,sprintf('%d\n',0));
            fclose(fid1);
            copyfile(fullfile(paraDir0,'config.txt'), fullfile(inputDir));
            exit(0)
        end
    else
        if stereoParam.transVecRef(2) > 0
            fid1 = fopen(fullfile(inputDir,'report.txt'),'w');%????
            fprintf(fid1,sprintf('%d\n',0));
            fclose(fid1);
            copyfile(fullfile(paraDir0,'config.txt'), fullfile(inputDir));
            exit(0)
        end
    end
    
    
    [KL, KR, rectParamL, rectParamR, rotMatLeft, rotMatRight,  intrMatLeftNew, intrMatRightNew, meanError, KK_newL] = CheckStereoCalibReprojection(stereoParam, imgSize, imgListL, imgListR, camParamL, camParamR, cbcXYL, cbcXYR, cbGridL, cbGridR, markerIdList, xyzCellGlobal);
    cfg.KK_newL = KK_newL;
%     cfg.imgLRGB = imgL;
    cfg.lr_param.rectParamL = rectParamL;
    cfg.lr_param.rectParamR = rectParamR;
    cfg.lr_param.rotMatLeft = rotMatLeft;
    cfg.lr_param.rotMatRight = rotMatRight;
    cfg.lr_param.intrMatLeftNew = intrMatLeftNew;
    cfg.lr_param.intrMatRightNew = intrMatRightNew;
    
    imgL = imread(imgListL{1});
    imgR = imread(imgListR{1});
    imgRGB = imread(imgListRGB{1});
    
    [matX, matY] = meshgrid(1:imgSize(2), 1:imgSize(1));
    pixRectL = remapRect([matX(:) matY(:)]', intrMatLeftNew, KL, stereoParam.kcLeft, rotMatLeft);
    pixRectR = remapRect([matX(:) matY(:)]', intrMatRightNew, KR, stereoParam.kcRight, rotMatRight);
    
    reMapLX = reshape(pixRectL(:,1), imgSize);
    reMapLY = reshape(pixRectL(:,2), imgSize);
    reMapRX = reshape(pixRectR(:,1), imgSize);
    reMapRY = reshape(pixRectR(:,2), imgSize);
    
    
     if 1 %~cfg.switch_lr
        rectImgL = uint8(interp2(matX, matY, double(imgL(:,:,1)),reMapLX,reMapLY));
        rectImgR = uint8(interp2(matX, matY, double(imgR(:,:,1)),reMapRX,reMapRY));
    else
        rectImgL = uint8(interp2(matX, matY, double(imgR(:,:,1)),reMapLX,reMapLY));
        rectImgR = uint8(interp2(matX, matY, double(imgL(:,:,1)),reMapRX,reMapRY));
     end
    
    
    save(fullfile(inputDir0, 'rectMap.mat'), 'matX', 'matY', 'reMapLX', 'reMapLY', 'reMapRX', 'reMapRY');
    % 量化左目到RGB的标定精度，计算左目到RGB的重投影误差
    stereoParamRGB = CheckL2RGBReprojection(imgListL, imgListRGB, camParamRGB, cbcXYRGB, cbGridL, cbGridRGB, markerIdList, xyzCellGlobal, configL, configRGB);
    
    inreMatRGB_old = [stereoParamRGB.focRight(1) 0 stereoParamRGB.cenRight(1);0 stereoParamRGB.focRight(2) stereoParamRGB.cenRight(2);0 0 1];
     [matRGBX, matRGBY] = meshgrid(1:imgSizeRGB(2), 1:imgSizeRGB(1));
    
    pixRectRGB = remapRect([matRGBX(:) matRGBY(:)]', cfg.RGBIntrMat, inreMatRGB_old, stereoParamRGB.kcRight, eye(3));
    reMapRGBX = reshape(pixRectRGB(:,1), imgSizeRGB);
    reMapRGBY = reshape(pixRectRGB(:,2), imgSizeRGB);
    save(fullfile(inputDir0, 'rectMap.mat'), 'matX', 'matY', 'reMapLX', 'reMapLY', 'reMapRX', 'reMapRY','imgSize', 'imgSizeRGB','matRGBX','matRGBY','reMapRGBX','reMapRGBY','cfg');
    
    cfg.imgL = imgL;
    cfg.imgR = imgR;
    cfg.imgRGB = imgRGB;
    cfg.stereoParamRGB = camParamRGB;
    ConvertParam2(inputDir0, stepDir, stereoParam, imgSize, 1);
%     imgSizeRGB = size(imgRGB(:,:,1));
    if 0
        ConvertParam2(fullfile(inputDir0,'RGB'), stepDir, camParamRGB, imgSize,0);
    else
        ConvertParam2(inputDir0, stepDir, camParamRGB, imgSizeRGB,0);
    end
end


% 生成report.txt，文件中有3个数，则第1个数是重投影误差（像素，越小越好），第2，3个数是棋盘格x，y方向的误差（毫米，越小越好）。
if ~exist('meanError', 'var')
    meanError = [0.1 0.1 0.1];
end
fid1 = fopen(fullfile(inputDir,'report.txt'),'w');%????
if 0
    fprintf(fid1,sprintf('%0.5f %0.5f %0.5f\n',meanError(1),meanError(2), meanError(3)));
else
    for v = 1 : length(meanError)
        if v < length(meanError)
            fprintf(fid1,sprintf('%0.5f ',meanError(v)));
        else
            fprintf(fid1,sprintf('%0.5f\n',meanError(v)));
        end
    end
end
fclose(fid1);

copyfile(fullfile(paraDir0,'config.txt'), fullfile(inputDir));

end


function varargout = Deoutlier(leftTrack, camParamL, cbcXYL, cbGridL, rightTrack, camParamR, cbcXYR, cbGridR, varargin)
global cfg

if (nargin > 8)
    rgbTrack =  varargin{1};
    camParamRGB =  varargin{2};
    cbcXYRGB =  varargin{3};
    cbGridRGB =  varargin{4};

else
    rgbTrack = [];
end

if isempty(rgbTrack)
    
    K_pinholeL = [camParamL.foc(1) 0 camParamL.cen(1);0 camParamL.foc(2) camParamL.cen(2);0 0 1];
    K_pinholeR = [camParamR.foc(1) 0 camParamR.cen(1);0 camParamR.foc(2) camParamR.cen(2);0 0 1];
    for u = 1 : length(cbcXYL)
        undistPixNormL = normalize_pixel(cbcXYL{u},camParamL.foc,camParamL.cen,camParamL.kc,camParamL.alpha);
        undistPixL = pflat((K_pinholeL)*pextend(undistPixNormL));
        undistPixL = undistPixL(1:2,:)';
        
        undistPixNormR = normalize_pixel(cbcXYR{u},camParamR.foc,camParamR.cen,camParamR.kc,camParamR.alpha);
        undistPixR = pflat((K_pinholeR)*pextend(undistPixNormR));
        undistPixR = undistPixR(1:2,:)';
        
        poseL = [rodrigues(camParamL.rotVec(:, u)) camParamL.tranVec(:, u); 0 0 0 1];
        poseR = [rodrigues(camParamR.rotVec(:, u)) camParamR.tranVec(:, u); 0 0 0 1];
        
        pt3dL = cbGridL{u}';
        pt3dL(:,3) = 0;
        pt3dR = cbGridR{u}';
        pt3dR(:,3) = 0;
        
        ptIcsL = TransformAndProject(pt3dL, K_pinholeL, poseL(1:3,1:3), poseL(1:3,4));
        ptIcsR = TransformAndProject(pt3dR, K_pinholeR, poseR(1:3,1:3), poseR(1:3,4));
        
        [~, errorL] = NormalizeVector(ptIcsL - undistPixL);
        [~, errorR] = NormalizeVector(ptIcsR - undistPixR);
        
        validDetect_ = find(errorL < cfg.reproj_pix_thr & errorR < cfg.reproj_pix_thr);
        
        goodIdAll = reshape(1:length(errorL),4,[]);
        validDetect = [];
        goodId = [];
        for j = 1 : size(goodIdAll,2)
           inId = ismember(validDetect_, goodIdAll(:,j));
            if sum(inId) == 4;
                goodId = [goodId; j];
                validDetect = [validDetect; validDetect_(inId)];
            end
        end
        
        
        leftTrack{u,1}{1} = leftTrack{u,1}{1}(goodId,:);
        leftTrack{u,1}{2}{1} = leftTrack{u,1}{2}{1}(goodId,:);
        leftTrack{u,1}{2}{2} = leftTrack{u,1}{2}{2}(goodId,:);
        leftTrack{u,1}{2}{3} = leftTrack{u,1}{2}{3}(goodId,:);
        leftTrack{u,1}{2}{4} = leftTrack{u,1}{2}{4}(goodId,:);
        
        rightTrack{u,1}{1} = rightTrack{u,1}{1}(goodId,:);
        rightTrack{u,1}{2}{1} = rightTrack{u,1}{2}{1}(goodId,:);
        rightTrack{u,1}{2}{2} = rightTrack{u,1}{2}{2}(goodId,:);
        rightTrack{u,1}{2}{3} = rightTrack{u,1}{2}{3}(goodId,:);
        rightTrack{u,1}{2}{4} = rightTrack{u,1}{2}{4}(goodId,:);
%         err{u, 1} = errorL;
    end
    varargout{1} = leftTrack;
    varargout{2} = rightTrack;
else
    K_pinholeL = [camParamL.foc(1) 0 camParamL.cen(1);0 camParamL.foc(2) camParamL.cen(2);0 0 1];
    K_pinholeR = [camParamR.foc(1) 0 camParamR.cen(1);0 camParamR.foc(2) camParamR.cen(2);0 0 1];
    K_pinholeRGB = [camParamRGB.foc(1) 0 camParamRGB.cen(1);0 camParamRGB.foc(2) camParamRGB.cen(2);0 0 1];
    for u = 1 : length(cbcXYL)
        undistPixNormL = normalize_pixel(cbcXYL{u},camParamL.foc,camParamL.cen,camParamL.kc,camParamL.alpha);
        undistPixL = pflat((K_pinholeL)*pextend(undistPixNormL));
        undistPixL = undistPixL(1:2,:)';
        
        undistPixNormR = normalize_pixel(cbcXYR{u},camParamR.foc,camParamR.cen,camParamR.kc,camParamR.alpha);
        undistPixR = pflat((K_pinholeR)*pextend(undistPixNormR));
        undistPixR = undistPixR(1:2,:)';
        
        undistPixNormRGB = normalize_pixel(cbcXYRGB{u},camParamRGB.foc,camParamRGB.cen,camParamRGB.kc,camParamRGB.alpha);
        undistPixRGB = pflat((K_pinholeRGB)*pextend(undistPixNormRGB));
        undistPixRGB = undistPixRGB(1:2,:)';
        
        poseL = [rodrigues(camParamL.rotVec(:, u)) camParamL.tranVec(:, u); 0 0 0 1];
        poseR = [rodrigues(camParamR.rotVec(:, u)) camParamR.tranVec(:, u); 0 0 0 1];
        poseRGB = [rodrigues(camParamRGB.rotVec(:, u)) camParamRGB.tranVec(:, u); 0 0 0 1];
        
        pt3dL = cbGridL{u}';
        pt3dL(:,3) = 0;
        pt3dR = cbGridR{u}';
        pt3dR(:,3) = 0;
        pt3dRGB = cbGridRGB{u}';
        pt3dRGB(:,3) = 0;
        
        ptIcsL = TransformAndProject(pt3dL, K_pinholeL, poseL(1:3,1:3), poseL(1:3,4));
        ptIcsR = TransformAndProject(pt3dR, K_pinholeR, poseR(1:3,1:3), poseR(1:3,4));
        ptIcsRGB = TransformAndProject(pt3dRGB, K_pinholeRGB, poseRGB(1:3,1:3), poseRGB(1:3,4));
        
        [~, errorL] = NormalizeVector(ptIcsL - undistPixL);
        [~, errorR] = NormalizeVector(ptIcsR - undistPixR);
        [~, errorRGB] = NormalizeVector(ptIcsRGB - undistPixRGB);
        validDetect_ = find(errorL < cfg.reproj_pix_thr & errorR < cfg.reproj_pix_thr & errorRGB < cfg.reproj_pix_thr);
        
        goodIdAll = reshape(1:length(errorL),4,[]);
        validDetect = [];
        goodId = [];
        for j = 1 : size(goodIdAll,2)
           inId = ismember(validDetect_, goodIdAll(:,j));
            if sum(inId) == 4;
                goodId = [goodId; j];
                validDetect = [validDetect; validDetect_(inId)];
            end
        end
        
        
        leftTrack{u,1}{1} = leftTrack{u,1}{1}(goodId,:);
        leftTrack{u,1}{2}{1} = leftTrack{u,1}{2}{1}(goodId,:);
        leftTrack{u,1}{2}{2} = leftTrack{u,1}{2}{2}(goodId,:);
        leftTrack{u,1}{2}{3} = leftTrack{u,1}{2}{3}(goodId,:);
        leftTrack{u,1}{2}{4} = leftTrack{u,1}{2}{4}(goodId,:);
        
        rightTrack{u,1}{1} = rightTrack{u,1}{1}(goodId,:);
        rightTrack{u,1}{2}{1} = rightTrack{u,1}{2}{1}(goodId,:);
        rightTrack{u,1}{2}{2} = rightTrack{u,1}{2}{2}(goodId,:);
        rightTrack{u,1}{2}{3} = rightTrack{u,1}{2}{3}(goodId,:);
        rightTrack{u,1}{2}{4} = rightTrack{u,1}{2}{4}(goodId,:);
        
        rgbTrack{u,1}{1} = rgbTrack{u,1}{1}(goodId,:);
        rgbTrack{u,1}{2}{1} = rgbTrack{u,1}{2}{1}(goodId,:);
        rgbTrack{u,1}{2}{2} = rgbTrack{u,1}{2}{2}(goodId,:);
        rgbTrack{u,1}{2}{3} = rgbTrack{u,1}{2}{3}(goodId,:);
        rgbTrack{u,1}{2}{4} = rgbTrack{u,1}{2}{4}(goodId,:);
%         err{u, 1} = errorL;
    end
    
    varargout{1} = leftTrack;
    varargout{2} = rightTrack;
    varargout{3} = rgbTrack;
end
    
end