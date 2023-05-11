function LoopClosure(intrMatL,t1, keyframeInfo)
global num_levels numIter kfSplit farZ minInlierNum icpInlierThr showICP useICP mergeDist minOvlpRatio icpErrThr mergeSize reprojErr...
pnpErrThr useSURF

reprojErr = {};

2.5; 5; % meter
minInlierNum = 20;

icpInlierThr = [0.4 0.2];
showICP = 0;
forceBreak = 0;
if 1
    maxTrans = 7; 5.0; 3; % [1.5 3]; % meter
    maxAng = 120; 50; 3;5;10; % degree
    keyframeApart = [1 10000]; [1 4]; [1 200000];2;
%     farZ = [1.2 4];
    farZ = [0.3 4];
    icpInlierThr = [0.85 0.85];
    loopStart = 290;
    startt = 1;
    mergeSize = 0.01;
    
    useSURF = 1;
    
    icpErrThr = 0.07; 0.05; 0.05; % meter
    pnpErrThr = 4; % pixel;
    
    mergeDist = 0.1; 0.2;  0.1; 0.2;  0.05; %meter
    minOvlpRatio = 0.7; 0.6; 0.7; 0.80; %0.5;
    useICP = 1;
elseif 0
    maxTrans = 0.3; 1.5; % [1.5 3]; % meter
    maxAng = 10; 5; 3;5;10; % degree
    keyframeApart = [1 10]; 2;
    farZ = [0.2 1.5];
else
    maxTrans = 2.0; 0.4; 1.5; % [1.5 3]; % meter
    maxAng = 30; 5; 3;5;10; % degree
    if 0
        keyframeApart = [1 30]; 2;
    else
        keyframeApart = [1 30000000]; 2;
    end
%     farZ = [0.2 2.0];
    farZ = [0.9 3];
    icpInlierThr = [0.85 0.85];
    
    startt = 1;
    mergeSize = 0.02;
    
    icpErrThr = 0.05; % meter
    pnpErrThr = 3;  % pixel
    
    mergeDist = 0.05; %meter
    minOvlpRatio = 0.80; %0.5;
    useICP = 1;
end
minFrameInterval = 100;
% keyframeApart = 5;2;
num_levels = 3;
numIter = 50; 10; 3;
kfSplit = 1;

max_level=5;
win_width = 31; 
win_height = 31;

global OUTPUT_ROOT RECORDED_BIN_ROOT DUMP_ROOT titleStr thetaListGT ShiftByCoord probPath doubleShift mean_x_err rou...
    shiftPixThr doubleShift2 err_x switchMethod rou_temp_list Rou UseDoubleShift...
    UseEx

ShiftByCoord = false; true;
doubleShift = false; true; false; true; false;true;
mean_x_err = [0 0];
rou = 10; 7; 10; 7; 10; 7; 10;7; 10; 7; 10;7; 8; 14;9; 14; 9; 7; 12; 10; 7; 10; 7; 6; 10; 7; 10; 7; 5; 8; 9; 10.1;
Rou = 0;  0.5; 0; 1; 0; 1; 1.5;

shiftPixThr = inf; -2; 3; 2;
doubleShift2 = true;
err_x = [0 0];
switchMethod = true;
rou_temp_list = [];
UseDoubleShift = true; false; true; false; true;false; true; false; true;
UseEx = false; true; false; true; false; true;

baseline = norm(t1);

if 0
    load('D:\Temp\zed_calib\result\stereo_param.mat');
    imgSize = [1080 1920]';
else
    load('D:\Temp\zed_calib\720\result\stereo_param.mat');
    imgSize = [720 1280]';
end


assert(~isempty(OUTPUT_ROOT) && ~isempty(RECORDED_BIN_ROOT) && ~isempty(DUMP_ROOT), 'The directories are not correctly set');
csaParamFilePath = fullfile(DUMP_ROOT, 'coordsys_alignment.txt');
csa = CoordSysAligner;
 Load(csa, csaParamFilePath);
robotCfgFilePath = fullfile(OUTPUT_ROOT, 'robot_config.prototxt');
robotConfig = ReadRobotConfig(robotCfgFilePath);

camParamFilePath = fullfile(DUMP_ROOT, 'camera_param.txt');
cam = BCStereoCameraModel(camParamFilePath);
CvtToRectifiedModel(cam);
vsl = VisualLocalizer(cam, csa, robotConfig.tracking_config);
if isempty(vsl.scaleLvl)
    vsl.scaleLvl = 0;
end

intrMat = intrMatL;

if 0
    R = eye(3);
    T =  t1';
    KL = intrMat; %[stereoParam.focLeft(1) 0 stereoParam.cenLeft(1); 0 stereoParam.focLeft(2) stereoParam.cenLeft(2); 0 0 1];
    KR = intrMat; %[stereoParam.focRight(1) 0 stereoParam.cenRight(1); 0 stereoParam.focRight(2) stereoParam.cenRight(2); 0 0 1];
    kcL = zeros(5,1); %stereoParam.kcLeft;
    kcR = zeros(5,1); %stereoParam.kcRight;
    [rotMat1ToCP, rotMat2ToCP, projMat1, projMat2] = OcvStereoRectify(KL, kcL, KR, kcR, imgSize, rodrigues(R), T);
end

R_ = rodrigues(stereoParam.rotVecRef);
T =  stereoParam.transVecRef./1000;
KL = [stereoParam.focLeft(1) 0 stereoParam.cenLeft(1); 0 stereoParam.focLeft(2) stereoParam.cenLeft(2); 0 0 1];
KR = [stereoParam.focRight(1) 0 stereoParam.cenRight(1); 0 stereoParam.focRight(2) stereoParam.cenRight(2); 0 0 1];
kcL = stereoParam.kcLeft;
kcR = stereoParam.kcRight;
[rotMat1ToCP, rotMat2ToCP, projMat1, projMat2] = OcvStereoRectify(KL, kcL, KR, kcR, imgSize, rodrigues(R_), T);



vsl.camModel.monoCamModel1.focLen = [KL(1,1); KL(2,2)];
vsl.camModel.monoCamModel1.princpPt = [KL(1,3); KL(2,3)];
vsl.camModel.monoCamModel1.distCoeff = kcL;
vsl.camModel.monoCamModel2.focLen = [KR(1,1); KR(2,2)];
vsl.camModel.monoCamModel2.princpPt = [KR(1,3); KR(2,3)];
vsl.camModel.monoCamModel2.distCoeff = kcR;
SetPinholeCamParam(vsl.camModel.monoCamModel2, [projMat2(1,1); projMat2(2,2)], projMat2(1:2,3), projMat2(1,2));
SetPinholeCamParam(vsl.camModel.monoCamModel1, [projMat1(1,1); projMat1(2,2)], projMat1(1:2,3), projMat1(1,2));
vsl.camModel.imgSize = imgSize;
[mapX1, mapY1] = OcvInitUndistortRectifyMap(Get(vsl.camModel, 'IntrMat', 1), Get(vsl.camModel, 'DistCoeff', 1), rotMat1ToCP, Get(vsl.camModel, 'PinholeIntrMat', 1), vsl.camModel.imgSize);
[mapX2, mapY2] = OcvInitUndistortRectifyMap(Get(vsl.camModel, 'IntrMat', 2), Get(vsl.camModel, 'DistCoeff', 2), rotMat2ToCP, Get(vsl.camModel, 'PinholeIntrMat', 2), vsl.camModel.imgSize);
SetMap(vsl.camModel.monoCamModel1, mapX1, mapY1);
SetMap(vsl.camModel.monoCamModel2, mapX2, mapY2);
vsl.camModel.rotVec1To2 = rodrigues(R_);
vsl.camModel.transVec1To2 = T;

intrMat_vsl = Get(vsl.camModel, 'PinholeIntrMat', 1,vsl.scaleLvl);


vsl.featPtTracker.configParam.max_level = max_level;
vsl.featPtTracker.configParam.win_width = win_width; 
vsl.featPtTracker.configParam.win_height = win_height;



keyframeNum = size(keyframeInfo, 1);
% % maxTrans = 0.5; % [1.5 3]; % meter
% % maxAng = 3; % degree
% % minFrameInterval = 100;
% % keyframeApart = 2;

loop_closure = {};

for i = startt : keyframeNum  % 200
    if i > loopStart
        keyframeApart = [1 10000000];
    end
    %    for j = max(1,i + keyframeApart(3)) : i
    if forceBreak
        break;
    end
%     for j = max(1,-10) : i
    for j = max(startt ,i - keyframeApart(2)-3) : i
       first_pose = keyframeInfo{j,2}(:,1:4);
       second_pose = keyframeInfo{i,2}(:,1:4);
       ref_frame = keyframeInfo{j,1}(1);
       cur_frame = keyframeInfo{i,1}(1);
       ref_pc = keyframeInfo{j,3};
       cur_pc = keyframeInfo{i,3};
       ref_disp = keyframeInfo{j,6};
       cur_disp = keyframeInfo{i,6};
       
       ref_img = keyframeInfo{j,4};
       cur_img = keyframeInfo{i,4};
       ref_pt2d = keyframeInfo{j,5};
       cur_pt2d = keyframeInfo{i,5};
       trans_dis = norm(first_pose(1:3,4) - second_pose(1:3,4));
       if 0
           first_z = first_pose(1:3, 3)';
           second_z = second_pose(1:3, 3)';
           ang = real(CalcDegree(first_z, second_z));
       else
           first_z = first_pose(1:3, 3)';
           second_z = second_pose(1:3, 3)';
           deltaRot = rodrigues(first_pose(1:3,1:3)'*second_pose(1:3,1:3));
           ang = rad2deg(norm(deltaRot));
       end
       
       if 1 % abs(i - j) < 500
           %            if ang < maxAng && trans_dis < maxTrans && dot(first_z, second_z) > 0 && abs(i - j) > keyframeApart % && abs(cur_frame - ref_frame) > minFrameInterval
           if ang < maxAng && trans_dis < maxTrans  && abs(i - j) >= keyframeApart(1) && abs(i - j) <= keyframeApart(2) % &&&& dot(first_z, second_z) > 0 && abs(cur_frame - ref_frame) > minFrameInterval
              
               posePrv = first_pose;
               poseCur = second_pose;
               
               imgPrv = ref_img;
               imgCur = cur_img;
               depthPrv = ref_pc.X';
               depthCur = cur_pc.X';
               pt2dPrv = ref_pt2d;
               pt2dCur = cur_pt2d;
               
               
               depPrv = intrMat(1,1)*baseline./ref_disp;
               depCur = intrMat(1,1)*baseline./cur_disp;
               
               poseNew = imageAlignment(vsl, [ref_frame cur_frame], imgPrv, imgCur, depthPrv, depthCur, intrMat, pt2dPrv, pt2dCur, (poseCur \ posePrv), depPrv, depCur);
               
%                if (abs(i-j) <= 2) && isempty(poseNew)
               if (abs(i-j) <= keyframeApart(1)+1) && isempty(poseNew)
                  poseNew = poseCur \ posePrv;
               end
               if ~isempty(poseNew)
                   loop_closure = [loop_closure; [{[ref_frame cur_frame ang trans_dis]} {first_pose} {second_pose} {ref_pc} {cur_pc} {ref_img} {cur_img} {ref_pt2d} {cur_pt2d} {ref_disp} {cur_disp} {poseNew}]];
               else
                   tfafgd = 1;
               end
               if 0
                   figure,imshow([ref_img cur_img]);title(sprintf('keyframeId: [%d, %d]',j, i));% title(sprintf('%d %d',ref_frame, cur_frame));
               end
           end
       else

       end
   end
end


figure,plot(cell2mat(reprojErr(:,2)))


loop_closure_bak = loop_closure;


[poseMat0, ptCloudScene0, poseMat, ptCloudScene] = processing(vsl, intrMat, baseline, loop_closure);




figure,pcshow(ptCloudScene0.Location,repmat([1 0 0],size(ptCloudScene0.Location,1) ,1));hold on;pcshow(ptCloudScene.Location,repmat([0 0 1],size(ptCloudScene.Location,1) ,1));legend('before','after');

save('loopData.mat', 'ptCloudScene0','ptCloudScene','poseMat0','poseMat','-v7.3')


if 0
    figure,pcshow(ptCloudScene); plotPath(poseMat);
end


end

function [poseMat, ptCloudScene] = ToPoseMat(vsl, jointPose, ptCloud)
global mergeSize
% mergeSize = 0.01;

poseMat = [];
for i = 1 : length(jointPose)
    if i == 1
        ptCloudScene = pointCloud(ptCloud{i,1}.X','Color',ptCloud{i,1}.C');
    else
        ptCloudCur = pointCloud(ptCloud{i,1}.X','Color',ptCloud{i,1}.C');
        L2G = jointPose{i,1};
        [UU,SS,VV] = svd(L2G(1:3,1:3));
        L2G(1:3,1:3) = UU*VV';
        L2G(4,:) = [0 0 0 1];
        
        ptCloudAligned = pctransform(ptCloudCur, affine3d(L2G'));
        ptCloudScene = pcmerge(ptCloudScene, ptCloudAligned, mergeSize);
    end
    poseMat = [poseMat; [reshape(jointPose{i,1}(1:3,1:3),1,9) [jointPose{i,1}(1:3,4)']]]; 
end

end
function poseNew = imageAlignment(obj, frameIdPair, imgPrv, imgCur, depthPrv, depthCur, intrMat, pt2dPrv, pt2dCur, initPose, depPrv, depCur)
global num_levels farZ minInlierNum icpInlierThr showICP useICP minOvlpRatio icpErrThr reprojErr pnpErrThr
para = 5;
initial_sigma = 0;
default_dof = 0;
% num_levels =5; 3; 1;

farz = farZ.*1000;

depthPrv = depthPrv.*1000;
depthCur = depthCur.*1000;


depPrv = depPrv.*1000;
depCur = depCur.*1000;


initPose(1:3,4) = initPose(1:3,4).*1000;




img_prev = (im2double(imgPrv));
img_curr = (im2double(imgCur));

imgSize = size(imgPrv);

depthPrv_temp2 = nan(imgSize);
indPrv = sub2ind(imgSize, round(pt2dPrv(:,2)), round(pt2dPrv(:,1)));
depthPrv_temp2(indPrv) = depthPrv(:,3);

[XYZ1] = GetXYZFromDepth(intrMat, pt2dPrv,depthPrv(:,3));
err1 = XYZ1 - depthPrv;


depthCur_temp2 = nan(imgSize);
indCur = sub2ind(imgSize, round(pt2dCur(:,2)), round(pt2dCur(:,1)));
depthCur_temp2(indCur) = depthCur(:,3);

[XYZ2] = GetXYZFromDepth(intrMat, pt2dCur,depthCur(:,3));
err2 = XYZ2 - depthCur;

depthPrv_temp = depPrv;
depthCur_temp = depCur;

try
%     dhshzhf;
    [rt_init, inId, matches] = TD(obj, imgPrv,imgCur, depthPrv_temp, depthCur_temp, intrMat, pt2dPrv, depthPrv(:,3));
    rt_init0 = rt_init;
    rt_init0(4,:) = [0 0 0 1];
    if length(inId) > minInlierNum
        poseNew = rt_init;
    else
        jfjsdjsd;
        poseNew = [];
        return;
    end
    % % % % % % % % % % % % % % % % % % % % %     try
    if 1
        depthPrv_bak = depthPrv;
        depthPrv0 = (rt_init(1:3,1:3)*depthPrv' + repmat(rt_init(1:3,4),1,size(depthPrv,1)))';
        %         tform = pcregrigid(pointCloud(depthPrv(depthPrv(:,3)<farz,:)), pointCloud(depthCur(depthCur(:,3)<farz,:)),'InlierRatio',0.6,'MaxIterations',100,'InitialTransform',affine3d(rt_init0'),'Extrapolate',true);
        depPrv00 = depthPrv0(depthPrv0(:,3)<farz(2)&depthPrv0(:,3)>farz(1),:);
        depCur00 = depthCur(depthCur(:,3)<farz(2)&depthCur(:,3)>farz(1),:);
        [depPrv000, depCur000, ovlpRatio] = FindMutualPart(depPrv00 ,depCur00);
        if ovlpRatio < minOvlpRatio
% %             ashgdg;
            poseNew = [];
            return;
        end
        [tform,~,rmse] = pcregrigid(pointCloud(depPrv000), pointCloud(depCur000),'InlierRatio',icpInlierThr(1),'MaxIterations',100,'Extrapolate',true); %,'Metric','pointToPlane');
        
        if rmse > 1000*icpErrThr
% %             sdlkjgf;
            poseNew = [];
            return;
        end
        tt = tform.T'*rt_init;
        poseNew = tt;
        pt2dReproj = TransformAndProject(matches(:,3:5), intrMat, poseNew(1:3,1:3), poseNew(1:3,4));
        error = matches(:,6:7) - pt2dReproj;
        
        if mean(mean(abs(error),1)) > pnpErrThr
                poseNew = [];
                return;
        end
        
        reprojErr = [reprojErr; {frameIdPair} {[error]}];
        %%
        if ~useICP
            poseNew = rt_init;
        end
        ptCurNew = (tt(1:3,1:3)*depthPrv' + repmat(tt(1:3,4),1,size(depthPrv,1)))';
         tt1 = tform.T';
%         poseNew = tt;
        ptCurNew1 = (tt1(1:3,1:3)*depthPrv0(depthPrv0(:,3)<farz(2)&depthPrv0(:,3)>farz(1),:)' + repmat(tt1(1:3,4),1,size(depthPrv0(depthPrv0(:,3)<farz(2)&depthPrv0(:,3)>farz(1),:),1)))';
        err = ptCurNew(depthPrv0(:,3)<farz(2)&depthPrv0(:,3)>farz(1),:) - ptCurNew1;
        if showICP
%             figure,subplot(1,2,1),pcshow(depthCur,repmat([1 0 0],size(depthCur,1),1));hold on;pcshow(depthPrv,repmat([0 0 1],size(depthPrv,1),1));legend('cur','prv');subplot(1,2,2),pcshow(depthCur,repmat([1 0 0],size(depthCur,1),1));hold on;pcshow(ptCurNew,repmat([0 0 1],size(ptCurNew,1),1));legend('cur','prv 2 cur');
            
            figure,subplot(1,2,1),imshow(imgPrv);title('prv');subplot(1,2,2),imshow(imgCur);title('cur');
            
           figure,subplot(2,2,1),pcshow(depthCur(depthCur(:,3)<farz(2)&depthCur(:,3)>farz(1),:),repmat([1 0 0],size(depthCur(depthCur(:,3)<farz(2)&depthCur(:,3)>farz(1),:),1),1));hold on;pcshow(depthPrv(depthPrv0(:,3)<farz(2)&depthPrv0(:,3)>farz(1),:),repmat([0 0 1],size(depthPrv(depthPrv0(:,3)<farz(2)&depthPrv0(:,3)>farz(1),:),1),1));legend('cur','prv');subplot(2,2,2),pcshow(depthCur(depthCur(:,3)<farz(2)&depthCur(:,3)>farz(1),:),repmat([1 0 0],size(depthCur(depthCur(:,3)<farz(2)&depthCur(:,3)>farz(1),:),1),1));hold on;pcshow(ptCurNew1,repmat([0 0 1],size(ptCurNew1,1),1));legend('cur','prv 2 cur');
            subplot(2,2,3),pcshow(depthCur(depthCur(:,3)<farz(2)&depthCur(:,3)>farz(1),:),repmat([1 0 0],size(depthCur(depthCur(:,3)<farz(2)&depthCur(:,3)>farz(1),:),1),1));hold on;pcshow(depthPrv0(depthPrv0(:,3)<farz(2)&depthPrv0(:,3)>farz(1),:),repmat([0 0 1],size(depthPrv0(depthPrv0(:,3)<farz(2)&depthPrv0(:,3)>farz(1),:),1),1));legend('cur','init prv 2 cur');title(sprintf('ovlp: %0.3f\n reprojErr: %0.3f',ovlpRatio, mean(mean(abs(error),1))));
            subplot(2,2,4),pcshow(depCur000,repmat([1 0 0],size(depCur000,1),1));hold on;pcshow(depPrv000,repmat([0 0 1],size(depPrv000,1),1));legend('cur','init prv 2 cur');title('mutual part');
        end
        
        hzdhsr = 1;
        
        
    end
% % % % % % % % % % % % % % % % % % % % %     catch
% % % % % % % % % % % % % % % % % % % % %         fgj=1;
% % % % % % % % % % % % % % % % % % % % %     end
catch
    try
        if 0
            poseNew = [];
            return;
        end
        rt_init = initPose;
        rt_init0 = initPose;
        rt_init0(4,:) = [0 0 0 1];
        depthPrv_bak = depthPrv;
        depthPrv0 = (rt_init(1:3,1:3)*depthPrv' + repmat(rt_init(1:3,4),1,size(depthPrv,1)))';
        depPrv00 = depthPrv0(depthPrv0(:,3)<farz(2)&depthPrv0(:,3)>farz(1),:);
        depCur00 = depthCur(depthCur(:,3)<farz(2)&depthCur(:,3)>farz(1),:);
        %     try
        [depPrv000, depCur000, ovlpRatio] = FindMutualPart(depPrv00 ,depCur00);
        if ovlpRatio < minOvlpRatio
            poseNew = [];
            return;
        end
        %     catch
        %         poseNew = [];
        %         return;
        %     end
        if 1
            %         tform = pcregrigid(pointCloud(depthPrv(depthPrv(:,3)<farz,:)), pointCloud(depthCur(depthCur(:,3)<farz,:)),'InlierRatio',0.2,'MaxIterations',100,'InitialTransform',affine3d(rt_init0'),'Extrapolate',true);
            [tform,~,  rmse] = pcregrigid(pointCloud(depPrv000), pointCloud(depCur000),'InlierRatio',icpInlierThr(2),'MaxIterations',100,'Extrapolate',true);
            
            if rmse > 1000*icpErrThr
                poseNew = [];
                return;
            end
            tt = tform.T'*rt_init;
            poseNew = tt;
            
            pt2dReproj = TransformAndProject(matches(:,3:5), intrMat, poseNew(1:3,1:3), poseNew(1:3,4));
            error = matches(:,6:7) - pt2dReproj;
            
            if mean(mean(abs(error),1)) > pnpErrThr
                poseNew = [];
                return;
            end
            
            reprojErr = [reprojErr; {frameIdPair} {[error]}];
            
            
            ptCurNew = (tt(1:3,1:3)*depthPrv' + repmat(tt(1:3,4),1,size(depthPrv,1)))';
            tt1 = tform.T';
            %         poseNew = tt;
            ptCurNew1 = (tt1(1:3,1:3)*depthPrv0(depthPrv0(:,3)<farz(2)&depthPrv0(:,3)>farz(1),:)' + repmat(tt1(1:3,4),1,size(depthPrv0(depthPrv0(:,3)<farz(2)&depthPrv0(:,3)>farz(1),:),1)))';
            err = ptCurNew(depthPrv0(:,3)<farz(2)&depthPrv0(:,3)>farz(1),:) - ptCurNew1;
            
            if showICP
                %             figure,subplot(1,2,1),pcshow(depthCur,repmat([1 0 0],size(depthCur,1),1));hold on;pcshow(depthPrv,repmat([0 0 1],size(depthPrv,1),1));legend('cur','prv');subplot(1,2,2),pcshow(depthCur,repmat([1 0 0],size(depthCur,1),1));hold on;pcshow(ptCurNew,repmat([0 0 1],size(ptCurNew,1),1));legend('cur','prv 2 cur');
                
                figure,subplot(1,2,1),imshow(imgPrv);title('prv');subplot(1,2,2),imshow(imgCur);title('cur');
                
                figure,subplot(2,2,1),pcshow(depthCur(depthCur(:,3)<farz(2)&depthCur(:,3)>farz(1),:),repmat([1 0 0],size(depthCur(depthCur(:,3)<farz(2)&depthCur(:,3)>farz(1),:),1),1));hold on;pcshow(depthPrv(depthPrv0(:,3)<farz(2)&depthPrv0(:,3)>farz(1),:),repmat([0 0 1],size(depthPrv(depthPrv0(:,3)<farz(2)&depthPrv0(:,3)>farz(1),:),1),1));legend('cur','prv');subplot(2,2,2),pcshow(depthCur(depthCur(:,3)<farz(2)&depthCur(:,3)>farz(1),:),repmat([1 0 0],size(depthCur(depthCur(:,3)<farz(2)&depthCur(:,3)>farz(1),:),1),1));hold on;pcshow(ptCurNew1,repmat([0 0 1],size(ptCurNew1,1),1));legend('cur','prv 2 cur');
                subplot(2,2,3),pcshow(depthCur(depthCur(:,3)<farz(2)&depthCur(:,3)>farz(1),:),repmat([1 0 0],size(depthCur(depthCur(:,3)<farz(2)&depthCur(:,3)>farz(1),:),1),1));hold on;pcshow(depthPrv0(depthPrv0(:,3)<farz(2)&depthPrv0(:,3)>farz(1),:),repmat([0 0 1],size(depthPrv0(depthPrv0(:,3)<farz(2)&depthPrv0(:,3)>farz(1),:),1),1));legend('cur','init prv 2 cur');title(sprintf('ovlp: %0.3f\n reprojErr: %0.3f',ovlpRatio, mean(mean(abs(error),1))));  % title(num2str(ovlpRatio));
                subplot(2,2,4),pcshow(depCur000,repmat([1 0 0],size(depCur000,1),1));hold on;pcshow(depPrv000,repmat([0 0 1],size(depPrv000,1),1));legend('cur','init prv 2 cur');title('mutual part');
            end
            
            dfxjjzd = 1;
            
            
        else
            
            [pose_rel, score, warped_image, valid_mask, weightMap,  warped_z, warped_z_mask] = estimateVisualOdometry(img_curr, img_prev, depthCur_temp, depthPrv_temp, intrMat, num_levels, rt_init, 0, initial_sigma, default_dof, para);
            [warped_image1, valid_mask1, warped_z1, warped_z_mask1] = warpImage(img_curr, depthPrv_temp, rt_init, intrMat);
            score1 = mean((warped_image1(valid_mask1) - img_prev(valid_mask1)).^2);
            poseNew = pose_rel;
        end
    catch
        poseNew = [];
        return;
    end
end
poseNew(1:3,4)=poseNew(1:3,4)./1000;


if 0
    figure(4),clf;subplot(2,2,1),imshowpair(warped_image, img_prev);title('opt');
%     hold on;plot(pt2dPrv(:,1), pt2dPrv(:,2),'.r');
    subplot(2,2,2),imshowpair(warped_image1, img_prev);title('init'); 
%     hold on;plot(pt2dPrv(:,1), pt2dPrv(:,2),'.r');
    subplot(2,2,3),imshow(abs(warped_image - img_prev).*valid_mask, []);title(num2str(score));subplot(2,2,4),imshow(abs(warped_image1 - img_prev).*valid_mask1);title(num2str(score1));drawnow;
    figure,imshow(img_prev, []);hold on;plot(pt2dPrv(:,1), pt2dPrv(:,2),'.r');
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
function [inId, rtTemp_all_in] = PnPDeoutlier(id,prvPt, curPt, zList, pixLenThr, intrMat, pnp_ang_est_max_margin)
[~, traceLenPix_all] = NormalizeVector(prvPt - curPt);
validPnP_all = find(traceLenPix_all > pixLenThr);
[XYZ_Prv_all] = GetXYZFromDepth(intrMat, prvPt, zList);

inlierRate = length(validPnP_all)/length(traceLenPix_all);
inlierRateThr = 0.7;
if inlierRate > inlierRateThr
    [rtTemp_all, outlierId_all] = posest(double(curPt(validPnP_all,:)), double(XYZ_Prv_all(validPnP_all,:)), 0.8, intrMat, 'repr_err');
    
else

    try
        validPnP_all = find(traceLenPix_all > 0);
        [rtTemp_all, outlierId_all] = posest(double(curPt(validPnP_all,:)), double(XYZ_Prv_all(validPnP_all,:)), 0.8, intrMat, 'repr_err');
        
        
    catch
        validPnP_all = find(traceLenPix_all > 0);
        [rtTemp_all, outlierId_all] = posest(double(curPt(validPnP_all,:)), double(XYZ_Prv_all(validPnP_all,:)), 0.8, intrMat, 'repr_err');
        inId = id;
        asdgfk = 1;
        rtTemp_all_in = rtTemp_all;
        return;
    end
end
ptIcsTemp_gt_all = TransformAndProject(XYZ_Prv_all, intrMat, rodrigues(rtTemp_all(1:3)), rtTemp_all(4:6));
c2cErr_all = ptIcsTemp_gt_all - curPt;
[~, trackingEr_all] = NormalizeVector(c2cErr_all);
inTrackId_all = find(trackingEr_all < pnp_ang_est_max_margin(3));

if (nargout > 1)
    [rtTemp_all_in, outlierId_all_in] = posest(double(curPt(inTrackId_all,:)), double(XYZ_Prv_all(inTrackId_all,:)), 0.8, intrMat, 'repr_err');
else
    rtTemp_all_in = [];
end

inId = id(inTrackId_all);
end
function [rt_init, inId, matches] = TD(obj, imgPrv,imgCur, depthPrv, depthCur, intrMat, pixPrv, zList)
global minInlierNum useSURF
% [marg, ~, ~, ~] = GetPredictMargin(obj.featPtTracker);
pixPrv0 = pixPrv;
pixPrv0 = round(pixPrv0);

marg = 2;
pixLenThr = 5; 25; 1;
deOutlierThr = 5;
pnp_ang_est_max_margin  = [deg2rad([-1 1]) deOutlierThr];
b2cPmat = GetPctBody2Cam(obj.coordSysAligner, 1);
b2cPmat.transformMat = eye(4);
if 0
    [pixPrv, zList] = DetectValidFeature(imgPrv, depthPrv, intrMat, 0.02, marg);
end

%  zList = depthPrv(:,3);
if ~useSURF
    [predPtList, inTrackFlag, ~] = TrackFeaturePoints(obj.featPtTracker, imgPrv, imgCur, pixPrv, [], intrMat, b2cPmat);
else
    
    if 0
        [pixCur, zListCur] = DetectValidFeature(imgCur, depthCur, intrMat, 0.02, marg);
        
        [features1,valid_points1] = extractFeatures(imgPrv,pixPrv);
        [features2,valid_points2] = extractFeatures(imgCur,pixCur);
        indexPairs = matchFeatures(features1,features2,'MaxRatio', 0.5);
        
        predPtList = pixPrv;
        predPtList(indexPairs(:,1),:) = pixCur(indexPairs(:,2),:);
        inTrackFlag = false(size(pixPrv, 1), 1);
        inTrackFlag(indexPairs(:,1),:) = true;
        
    else
 
        cameraParams = cameraParameters('WorldUnits', 'm',...
            'IntrinsicMatrix', intrMat', ...
            'RotationVectors', zeros(1, 3), ...
            'TranslationVectors', [0 0 0]);
        
        
        [features_prev, validPoints_prev] = extract_features(imgPrv);
        [features_curr, validPoints_curr] = extract_features(imgCur);
        
        matchedIdx = matchFeatures(features_prev, features_curr, 'Unique', true, ...
            'Method', 'Approximate', 'MatchThreshold', .8);
        
        matchedPoints1 = validPoints_prev(matchedIdx(:, 1));
        matchedPoints2 = validPoints_curr(matchedIdx(:, 2));
        
        validPrv = FindValidDepth(matchedPoints1.Location, depthPrv, marg);
        pixPrv = matchedPoints1.Location(validPrv > 0, :);
        predPtList = matchedPoints2.Location(validPrv > 0, :);
        zList = validPrv(validPrv > 0);
        inTrackFlag = true(size(pixPrv,1),1);
    end
    
    
    if 0
        figure,showMatchedFeatures(imgPrv, imgCur, pixPrv(inTrackFlag,:), predPtList(inTrackFlag,:))
        figure,showMatchedFeatures(imgPrv, imgCur, pixPrv(indexPairs(:,1),:), pixCur(indexPairs(:,2),:))
    end
end
if 0
    [~,id1,id2] = intersect(pixPrv, pixPrv_);
end
    inTrackFlagId = find(inTrackFlag);
    if ~useSURF
        [fundMat, inn] = estimateFundamentalMatrix(pixPrv(inTrackFlag,:), predPtList(inTrackFlag,:), 'Method', 'RANSAC', 'InlierPercentage', 0.7);
        %         [fundMat, inn] = estimateFundamentalMatrix(pixPrv(inTrackFlag,:), predPtList(inTrackFlag,:), 'Method', 'LMedS');
        
    else
        
     
        
        [relativeOrient, relativeLoc, inn, status] = estimate_relative_motion(...
            matchedPoints1(validPrv > 0), matchedPoints2(validPrv > 0), cameraParams);
    end
    if 1
        inId = inTrackFlagId(inn);
    else
        inId_ = inTrackFlagId(inn);
        inId = intersect(inId_, id1);
    end
    if 0
        [~, traceLenPix] = NormalizeVector(pixPrv(inId,:) - predPtList(inId,:));
        validPnP = find(traceLenPix > pixLenThr);
    end
    if 0
        figure,showMatchedFeatures(imgPrv, imgCur, pixPrv(inId,:), predPtList(inId,:))
    end
    
%     [k2cBodyRotAng_2,Err_2,inId, reprojErrList] = VisualLocalizer.RotationAngleEstimate_PNP_(pnp_ang_est_max_margin,intrMat,pixPrv(inTrackFlag,:),predPtList(inTrackFlag,:),zList(inTrackFlag), find(inTrackFlag),true(sum(inTrackFlag),1),b2c,double(thetaP2C));
    
    inId0 = inId;
    vldPtCur = FindValidDepth(predPtList(inId0,:), depthCur, marg);
    
    inId = inId0(vldPtCur > 0,:);
    
    
    [~, traceLenPix] = NormalizeVector(pixPrv(inId,:) - predPtList(inId,:));
    validPnP = find(traceLenPix > pixLenThr);
    
    vldPtCurValid = vldPtCur(ismember(inId0, inId),:);
    
    pixPrvValid = pixPrv(inId,:);
    pixCurValid = predPtList(inId,:);
    zListValid = zList(inId,:);
    %     [pixCurReproj] = VisualLocalizer.GetGtTrace2(b2cPmat, k2cBodyRotAng_2, [pixPrvValid],zListValid,intrMat);
    
    [XYZ_Prv] = GetXYZFromDepth(intrMat, pixPrvValid, zListValid);
    try
        [rtTemp, outlierId] = posest(double(pixCurValid(validPnP,:)), double(XYZ_Prv(validPnP,:)), 0.8, intrMat, 'repr_err');
    catch
        validPnP = find(traceLenPix > 0);
        [rtTemp, outlierId] = posest(double(pixCurValid(validPnP,:)), double(XYZ_Prv(validPnP,:)), 0.8, intrMat, 'repr_err');
        
    end
    ptIcsTemp_gt = TransformAndProject(XYZ_Prv, intrMat, rodrigues(rtTemp(1:3)), rtTemp(4:6));
    c2cErr = ptIcsTemp_gt - pixCurValid;
    [~, trackingEr] = NormalizeVector(c2cErr);
    inTrackId = find(trackingEr < pnp_ang_est_max_margin(3) & traceLenPix >pixLenThr);
    inId_bak = inId;
    inId = inId(inTrackId);
    if ~useSURF
        rt_init = [rodrigues(rtTemp(1:3)) rtTemp(4:6); 0 0 0 1];
    else
        
        tNorm = rtTemp(4:6)./norm(rtTemp(4:6));
        [~, id1] = max(abs(tNorm));
        relativeLoc_ = sign(tNorm(id1))*relativeLoc;
        relativeLoc_s = relativeLoc_*norm(rtTemp(4:6));
        rt_init = [relativeOrient relativeLoc_s'; 0 0 0 1];
        
        
        ptIcsTemp_surf = TransformAndProject(XYZ_Prv, intrMat, rt_init(1:3,1:3), rt_init(1:3,4));
        c2cErr_surf = ptIcsTemp_surf - pixCurValid;
        
    end
    
    
    if length(inId) > minInlierNum
        if 0
            figure,showMatchedFeatures(imgPrv, imgCur, pixPrv(inId,:), predPtList(inId,:))
        end
    end
    [XYZ_Prv_inId] = GetXYZFromDepth(intrMat, pixPrv(inId,:), zList(inId,:));
    matches = [pixPrv(inId,:) XYZ_Prv_inId predPtList(inId,:)];
end
function vldPtCur = FindValidDepth(predPtList000, depthMapCurGT, marg)

vldPtCur = nan(size(predPtList000,1),1);
vldPrTrace = find(predPtList000(:,1) > marg & predPtList000(:,2) > marg & predPtList000(:,1) < size(depthMapCurGT,2) - marg & predPtList000(:,2) < size(depthMapCurGT,1) - marg);

vldPtCur_ = sub2ind(size(depthMapCurGT), round(predPtList000(vldPrTrace,2)), round(predPtList000(vldPrTrace,1)));
vldPtCur_depth = depthMapCurGT(vldPtCur_);

vldPtCur(vldPrTrace) = vldPtCur_depth;


end
function [pixKey, zList] = DetectValidFeature(imgPrv, depthMap, intrMat, thr, marg)

if ndims(imgPrv) == 3
    imgPrv = rgb2gray(imgPrv);
end
ptIcsKey = detectFASTFeatures(imgPrv, 'MinQuality', thr, 'MinContrast', thr);
pixKey0 = ptIcsKey.Location;
mask = ~isinf(depthMap) & depthMap > 0 & ~isnan(depthMap);
% se = strel('square',[7]);
se = strel('square',[1]);
mask_ = imerode(mask,se);
if 0
    ind0 = find(mask_ > 0);
else
    ind0 = find(mask > 0);
end
[y, x] = ind2sub(size(depthMap), ind0);
[xyzKey_all] = GetXYZFromDepth(intrMat, [x y],depthMap(ind0));


ind1 = sub2ind(size(depthMap), round(pixKey0(:,2)), round(pixKey0(:,1)));
ind = intersect(ind0, ind1);
pixKey = pixKey0(ismember(ind1, ind), :);

inBndFlag = pixKey(:, 1) >= marg + 1 & ...
    pixKey(:, 1) <= Cols(imgPrv) - marg & ...
    pixKey(:, 2) >= marg + 1 & ...
    pixKey(:, 2) <= Rows(imgPrv) - marg;

pixKey = pixKey(inBndFlag,:);
ind_check = sub2ind(size(depthMap), round(pixKey(:,2)), round(pixKey(:,1)));
zList = depthMap(ind_check);
end
function jointPose = pgo(poses, odom, Keys)
global numIter kfSplit
numPoses = size(poses, 1);
numConnections = size(odom, 1);
connection_pair = cell2mat(Keys');
poses_opt = poses;
A = sparse(7 * (numConnections + 1), 7 * numPoses);
A((7 * numConnections + 1):end, 1:7) = eye(7);
b = zeros(7 * (numConnections + 1), 1);
% numIter = 2;
% disp('Optimizing ...')
 fprintf('Optimizing ...[%2d]\n',  numIter)
for i = 1:numIter
%     fprintf('[%2d/%2d]\n', i, numIter)
    for k = 1:numConnections
        %idx1 = find(vs.Views.ViewId==vs.Connections.ViewId1(k));
        %idx2 = find(vs.Views.ViewId==vs.Connections.ViewId2(k));
        
        idx1 = connection_pair(k,1);
        idx2 = connection_pair(k,2);

        if idx1 == 1
            p1 = randn(7, 1) * 1e-8;
        else
            p1 = poses_opt(idx1, :)';
        end
        p2 = poses_opt(idx2, :)';

        J = calc_measurement_jacob(p1, p2);

        A((7 * k - 6):(7 * k), (7 * idx1 - 6):(7 * idx1)) = ...
            J(:, 1:7);
        A((7 * k - 6):(7 * k), (7 * idx2 - 6):(7 * idx2)) = ...
           J(:, 8:14);

        odom_hat = calc_odom(p1, p2);
        
        d = abs(int64(idx1) - int64(idx2)); 
        if d > 1 && d < kfSplit
            l = norm(odom_hat(5:7), 2);
            odom_hat(5:7) = odom_hat(5:7) / l;
        end
% 
% %         l = norm(odom_hat(5:7), 2);
% %         odom_hat(5:7) = odom_hat(5:7) / l;
              
        % disp([odom_hat'; odom(k, :)])
            
        b((7 * k - 6):(7 * k)) = odom(k, :)' - odom_hat;
    end
    
    delta = A \ b;
    delta = reshape(delta, [7, numPoses])';
    poses_opt = poses_opt + 0.5 * delta;
end



for i = 1 : size(poses_opt,1)
   jointPose{i,1} =  [rodrigues(poses_opt(i,2:4))' poses_opt(i,5:7)';0 0 0 1];
    
end

%%
order = colamd(A);
L_reordered = chol(A(:, order)' * A(:, order));
L_original = chol(A' * A);
if 0
figure(1)
% set(gcf, 'position', [100, 100, 800, 400])
subplot(1, 2, 1)
spy(L_original)
subplot(1, 2, 2)
spy(L_reordered)

end
% print(num2str(seq, 'seq%02d_R.png'), '-r300', '-dpng')


figure,
pcshow(poses(:, 5:7), 'VerticalAxis', 'Y', 'VerticalAxisDir', 'Down');
set(gca,  'CameraUpVector',[0 -1 0],'DataAspectRatio',[1 1 1]);
hold on;
%for k = 1:numConnections 
%idx1 = find(vs.Views.ViewId==vs.Connections.ViewId1(k));
%idx2 = find(vs.Views.ViewId==vs.Connections.ViewId2(k));
idx1 = connection_pair(:,1);
idx2 = connection_pair(:,2);
% plot(poses(:, 5), poses(:, 7), 'kx');hold on;
plot3(...
    [poses(idx1, 5)'; poses(idx2, 5)'], ...
    [poses(idx1, 6)'; poses(idx2, 6)'], ...
    [poses(idx1, 7)'; poses(idx2, 7)'], 'r')

plot3(poses(:, 5), poses(:, 6), poses(:, 7), '*k-','LineWidth',2,'MarkerSize',8);
plot3(poses(:, 5), poses(:, 6), poses(:, 7), 'ok-','LineWidth',2,'MarkerSize',6);
plot3(poses(1, 5), poses(1, 6), poses(1, 7), 'sg','LineWidth',10,'MarkerSize',20);
plot3(poses(end, 5), poses(end, 6), poses(end, 7), 'sy','LineWidth',10,'MarkerSize',20);
axis equal


end

function u = rod(R)
 u = [R(3, 2) - R(2, 3)
        R(1, 3) - R(3, 1)
        R(2, 1) - R(1, 2)];
    l_hat = norm(u); 
    if l_hat
        theta = acos((trace(R) - 1) / 2);
        u = theta * u / l_hat;
    else
        u = zeros(3, 1);
    end
end
function [pt00, pt11, ovlpRatio] = FindMutualPart(pt0 ,pt1)
global mergeDist
scal = 1000* mergeDist;
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
   if 0
       figure,pcshow(pt0(samePtFlag2,:),repmat([1 0 0],sum(samePtFlag2),1),'MarkerSize',100);hold on;pcshow(pt1(samePtFlag1,:),repmat([0 0 1],sum(samePtFlag1),1),'MarkerSize',100);
   end
   ovlpRatio = max(sum(samePtFlag2)/length(samePtFlag2), sum(samePtFlag1)/length(samePtFlag1));
end

%%
function [poseMat0, ptCloudScene0, poseMat, ptCloudScene] = processing(vsl, intrMat, baseline, loop_closure)
loop_closure_bak = loop_closure;
aa = cell2mat(loop_closure(:,1));
aaa = aa(:,1:2);
connections = (unique(aaa(:)));
ofst = 10000;
id_matches = [connections ofst + [1:length(connections)]'];
ptCloud ={};

poses = zeros(length(connections), 7);
for k = 1 : length(connections)
    id_1 = find(aa(:,1) == connections(k));
    if ~isempty(id_1)
        ptCloud{k,1} = loop_closure_bak{id_1(1), 4};
        jointPose0{k,1} = loop_closure_bak{id_1(1), 2};
    else
        id_2 = find(aa(:,2) == connections(k));
        ptCloud{k,1} = loop_closure_bak{id_2(1), 5};
        jointPose0{k,1} = loop_closure_bak{id_2(1), 3};
    end
    R = jointPose0{k,1}(1:3,1:3)';
    poses(k, 1) = log(det(R));
    poses(k, 2:4) = rod(R);
    poses(k, 5:7) =  jointPose0{k,1}(1:3,4);
end

Z = {};
Keys = {};
ct = 1;
% chosenFrameList = 1 :30: size(vslMat,1)-1;
odom = zeros(size(loop_closure,1), 7);
for g = 1 : size(loop_closure,1)
    if 0
        posePrv = [reshape(vslMat(chosenFrameList(g),1:9),3,3) vslMat(chosenFrameList(g),10:12)';0 0 0 1];
        poseCur = [reshape(vslMat(chosenFrameList(g+1),1:9),3,3) vslMat(chosenFrameList(g+1),10:12)';0 0 0 1];
    end
    posePrv = loop_closure{g, 2};
    poseCur = loop_closure{g, 3};
    
    if 0
        Keys{1,ct} = [loop_closure{g,1}(1:2)];
    else
        id_prv = loop_closure{g,1}(1);
        id_cur = loop_closure{g,1}(2);
        id1 = find(id_matches(:,1) == id_prv);
        id2 = find(id_matches(:,1) == id_cur);
        Keys{1,ct} = [(id_matches(id1,2) - ofst) (id_matches(id2,2) - ofst)];
    end
    if 0
        Z{ct,1} = inv(poseCur \ posePrv);
    else
        imgPrv = loop_closure{g,6};
        imgCur = loop_closure{g,7};
        depthPrv = loop_closure{g,4}.X';
        depthCur = loop_closure{g,5}.X';
        pt2dPrv = loop_closure{g,8};
        pt2dCur = loop_closure{g,9};
        
        depPrv = intrMat(1,1)*baseline./loop_closure{g,10};
        depCur = intrMat(1,1)*baseline./loop_closure{g,11};
        
        % prev to curr   /  golbal to local
        if 0
            try
                poseNew = imageAlignment(vsl, imgPrv, imgCur, depthPrv, depthCur, intrMat, pt2dPrv, pt2dCur, (poseCur \ posePrv), depPrv, depCur);
            catch
                continue
            end
        else
            poseNew = loop_closure{g,12};
        end
        Z{ct,1} = inv(poseNew);
        
        poseNew_ = inv(poseNew); % local to global
        
        RR = poseNew_(1:3,1:3)';
        tt = poseNew_(1:3,4)';
        odom(ct, 1) = log(det(RR));
        odom(ct, 2:4) = rod(RR);
        odom(ct, 5:7) = tt./norm(1);
        
    end
        ct = ct+1;
end

bb =cell2mat(Keys');
[~,idx] = sort(bb(:,1));
cc = bb(idx,:);

[poseMat0, ptCloudScene0] = ToPoseMat(vsl, jointPose0, ptCloud);

if 0
    jointPose = joint_optimization(Z, Keys, length(connections));
    [poseMat, ptCloudScene] = ToPoseMat(vsl, jointPose, ptCloud);
else
    jointPose = pgo(poses, odom, Keys);
    [poseMat, ptCloudScene] = ToPoseMat(vsl, jointPose, ptCloud);
end



figure,plot(aaa(:,1));


figure,subplot(1,2,1);pcshow(ptCloudScene0);hold on;
plot3(poseMat0(:,10), poseMat0(:,11), poseMat0(:,12),'-r','LineWidth',4,'MarkerSize',4);
plot3(poseMat0(1, 10), poseMat0(1, 11), poseMat0(1, 12), 'sg','LineWidth',10,'MarkerSize',10);
plot3(poseMat0(end, 10), poseMat0(end, 11), poseMat0(end, 12), 'sy','LineWidth',10,'MarkerSize',10);
% legend('before');
title('before');

subplot(1,2,2);pcshow(ptCloudScene);hold on;
plot3(poseMat0(:,10), poseMat0(:,11), poseMat0(:,12),'-r','LineWidth',4,'MarkerSize',4);
plot3(poseMat(:,10), poseMat(:,11), poseMat(:,12),'-b','LineWidth',4,'MarkerSize',4);
plot3(poseMat(1, 10), poseMat(1, 11), poseMat(1, 12), 'sg','LineWidth',10,'MarkerSize',10);
plot3(poseMat(end, 10), poseMat(end, 11), poseMat(end, 12), 'sy','LineWidth',10,'MarkerSize',10);
% legend('before','after');
title('after');


if 0
    plot3(poses(end, 5), poses(end, 6), poses(end, 7), 'sg','LineWidth',10,'MarkerSize',20);
end

end