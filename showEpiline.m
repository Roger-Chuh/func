function showEpiline()


close all

inputDir = 'C:\Users\Roger\Desktop\vins-fusion-fisheye\yvr\april\controller_2022-04-12-05-46-16';
inputDir = 'G:\matlab\data\tag\controller_2022-05-31-06-33-32';
inputDir = 'G:\matlab\data\hao\controller_2022-07-28-20-48-53';
inputDir = 'G:\matlab\data\hao\controller_2022-07-29-20-02-09';
inputDir = 'G:\matlab\data\hao\Device_D2HD232328D9000859_2022_08_01_23_20_34\controller_2022-08-01-23-18-28';
inputDir = 'G:\matlab\data\hao\controller_2022-08-05-15-33-54';
inputDir = 'G:\matlab\data\hui\Device_D2HD232328D9000957_2022_08_08_20_42_19\controller_2022-08-08-20-40-03';
inputDir = 'G:\matlab\data\close\controller_2022-04-12-05-33-20';
inputDir = 'G:\matlab\data\close\new_tag\hui\Device_D2HD432713D9000609_2022_08_16_17_07_46\controller_2022-08-15-20-09-14';

inputDir = 'G:\matlab\data\close\Device_ZTS26R1500AE_2022_09_16_18_13_31\controller_2022-06-16-03-31-11';

% inputDir = 'G:\matlab\data\hao\test\controller_2022-04-12-05-39-11';
% inputDir = 'G:\matlab\data\hao\test\controller_2022-04-12-05-40-31';
% inputDir = 'G:\matlab\data\hao\test\controller_2022-04-12-05-41-31';


% inputDir = 'G:\matlab\data\close\Device_D2HD131C17D9000719_2022_08_01_16_01_11\controller_2022-04-12-05-29-14';
% 
% inputDir = 'G:\matlab\data\close\ning\Device_D2HD232328D9000957_2022_08_01_17_56_51\controller_2022-08-01-07-21-35';
% inputDir = 'G:\matlab\data\close\hao\Device_D2HD232328D9000859_2022_08_01_19_10_57\controller_2022-08-01-19-08-13';


scale = 1;
camInfo = dir(fullfile(inputDir, 'Camera*'));

paraMat = load(fullfile(inputDir,'results', 'cam.txt'));
% paraMat(:,13:16) = paraMat(:,13:16)+1;
% if scale == 2
data1 = load(fullfile(inputDir, 'results','outlier_img','matlab_apriltag_corners.txt'));
name1 = (fullfile(inputDir, 'results','outlier_img','matlab_image_names.txt'));
% else
%     data1 = load(fullfile(inputDir, 'results','outlier_img1','matlab_apriltag_corners.txt'));
%     name1 = (fullfile(inputDir, 'results','outlier_img1','matlab_image_names.txt'));
% end
[FrameCam1, names1] = readNames(name1);
% [FrameCam2, names2] = readNames(name2);
frameList = (unique(FrameCam1(:,1)));
coviCamsEachTimestamp = {};
thr = 0;
err_ = [];
for i = 241 : 15 : length(frameList) 
%     if i == 468
%         sk = 1;
%     end
    frameId = frameList(i);
    frameData1 = find(data1(:,1) == frameId);
    %     frameData2 = find(data2(:,1) == frameId);
    camList = (unique(data1(frameData1,2)));
    coviCamsEachTimestamp{i,1} = i;
    coviCamsEachTimestamp{i,2} = camList;
    close all;
    cam_data = {};
    plate_common = [];
    for j = 1 : length(camList)
        camId = camList(j);
        frameCamData1 = find(data1(:,1) == frameId & data1(:,2) == camId);
        %         frameCamData2 = find(data2(:,1) == frameId & data2(:,2) == camId);
        frameCamImageData = find(FrameCam1(:,1) == frameId & FrameCam1(:,2) == camId);
        imageFileName = names1{frameCamImageData};
        img1 = imread(fullfile(inputDir, camInfo(camId+1).name,imageFileName));
        %         img2 = imresize(imread(fullfile(inputDir, camInfo(camId+1).name,imageFileName)),2);
        pixels1 = data1(frameCamData1,4:5);
        gridId = data1(frameCamData1,6);
        %         pixels2 = data2(frameCamData2,4:5);
        patternIdList1 = unique(data1(frameCamData1,3))';
        camIdList1 = unique(data1(frameCamData1,2));
        pattern_pixel = {};
        plate_common = [plate_common; patternIdList1'];
        for k = 1 : length(patternIdList1)
            pattern_pixel{k,1} = patternIdList1(k);
            pattern_pixel{k,2} = gridId(find(data1(frameCamData1,3) == patternIdList1(k)));
            pattern_pixel{k,3} = pixels1(find(data1(frameCamData1,3) == patternIdList1(k)),:);
        end
        cam_data{j,1} = camList(j);
        cam_data{j,2} =pattern_pixel;
        cam_data{j,3} =img1;
        %         patternIdList2 = unique(data2(frameCamData1,3));
        
        %         camIdList2 = unique(data2(frameCamData1,2));
        %         pixels21 = [(pixels2(:,1)+1.5)/2,(pixels2(:,2)+1.5)/2];
        %         [D] = pdist2(pixels1+1,pixels21,'euclidean');
        %         idx = find(D < 2);
        %         [y,x] = ind2sub(size(D),idx);
        %
        %         match = [pixels1(y,:)+1 pixels21(x,:)];
        % %         figure(2),clf;quiver(match(:,1),match(:,2),match(:,1)-match(:,3),match(:,2)-match(:,4));axis equal
        %         figure(j),clf;
        %         subplot(1,2,1);imshow(img1);hold on;plot(pixels1(:,1)+1,pixels1(:,2)+1,'.r');plot((pixels2(:,1)+1.5)/2,(pixels2(:,2)+1.5)/2,'.g');title(sprintf('frameid: %d, camid: %d, 640 x 480',frameId, camIdList1));legend('640 x 480', '1280 -> 640');drawnow;
        %         subplot(1,2,2);imshow(img2);hold on;plot(pixels2(:,1)+1,pixels2(:,2)+1,'.r');plot((pixels1(:,1)+1)*2-0.5,(pixels1(:,2)+1)*2-0.5,'.g');title(strcat('1280 x 960 ,pattern list: ',num2str(patternIdList1)));legend('1280 x 960', '640 -> 1280');drawnow;
        %
    end
    plate_common_ = plate_common;
    plate_common = unique(plate_common);
    cam_data;
    
    for plateId = 1 : length(plate_common)
        plate_id = plate_common(plateId);
        
        if (plate_id == 3)
            
        end
        
        cam_bundle = {};
        cam_bundle_cnt = 1;
        for iterCam = 1 : size(cam_data,1)
            try
                cam_id = cam_data{iterCam,1};
            catch
                sgkuh = 1;
            end
            for m = 1 : size(cam_data{iterCam,2},1)
                plate_ = cam_data{iterCam,2}{m,1};
                if(plate_ == plate_id)
                    cam_bundle{cam_bundle_cnt,1} = cam_id;
                    cam_bundle{cam_bundle_cnt,2} = cam_data{iterCam,2}{m,3};
                    cam_bundle{cam_bundle_cnt,3} = cam_data{iterCam,3};
                    cam_bundle_cnt = cam_bundle_cnt+1;
                end
            end
        end
        close all;
        camList = cell2mat(cam_bundle(:,1));
        checkVec = [0 2];
        if ~(sum((ismember(checkVec, camList)))>=2)
            continue;
        end
        if(size(cam_bundle,1) > 1)
            
            for cam1 = 1 : size(cam_bundle,1)-1
                cam_ind1 = cam_bundle{cam1,1}+1;
                cam_param1 = paraMat(cam_ind1,:);
                camEx1 = inv([reshape(cam_param1(1:9),3,3), cam_param1(10:12)'; 0 0 0 1]);
                K1 = [cam_param1(13) 0 cam_param1(15); 0 cam_param1(14) cam_param1(16); 0 0 1];
                dist1 = cam_param1(17:end);
                pix1 = cam_bundle{cam1,2};
%                 figure(cam1),imshow(cam_bundle{cam1,3});hold on;plot(pix1(:,1)+1, pix1(:,2)+1,'.r');title(sprintf('plate id: %d, cam id %d',plate_id,cam_bundle{cam1,1}));
                for cam2 = cam1+1 : size(cam_bundle,1)
                    
                    cam_ind2 = cam_bundle{cam2,1}+1;
%                     if ~(cam_ind1 == 2 && cam_ind2 == 3)
%                         continue;
%                     end

                    if ~(ismember(cam_ind1-1,checkVec) && ismember(cam_ind2-1,checkVec))
                        continue;
                    end
                    cam_param2 = paraMat(cam_ind2,:);
                    camEx2 = inv([reshape(cam_param2(1:9),3,3), cam_param2(10:12)'; 0 0 0 1]);
                    K2 = [cam_param2(13) 0 cam_param2(15); 0 cam_param2(14) cam_param2(16); 0 0 1];
                    dist2 = cam_param2(17:end);
                    pix2 = cam_bundle{cam2,2};
                    
                    for pointId = 1 : 1 : size(pix1,1)
                        pt3d = unprojectKB8( K1(1,1),  K1(2,2),  K1(1,3),  K1(2,3),  dist1(1),  dist1(2),  dist1(3),  dist1(4), pix1(pointId,:));
                        [pt2d, inImageFlag] = projectKB8([pt3d' 1], K1, dist1, eye(3), [0;0;0]);
                        [~,err] = NormalizeVector(pt2d - pix1(pointId,:));
                        if(pt2d(1) >thr && pt2d(1) < 640-thr && pt2d(2) > thr && pt2d(2) < 480-thr )
                            err_ =[err_;err ];
                        end
                        if(err>1.5)
                            sdfk = 1;
                        end
                        
                        T_rel = inv(camEx2) * camEx1;
                        PT = [];
                        XYZ_I = [];
                        for depthId = 0.1 : 0.05 : 1.5
                            XYZ = depthId * pt3d;
                            
                            XYZ_i = T_rel(1:3,1:3) * XYZ + T_rel(1:3,4);
                            [pt, inImageFlag] = projectKB8([XYZ_i' 1], K2, dist2, eye(3), [0;0;0]);
                            PT = [PT; pt];
                            XYZ_I = [XYZ_I; XYZ_i'];
                        end
                    figure(cam1),subplot(1,2,1);imshow(cam_bundle{cam1,3});hold on;plot(pix1(:,1)+1, pix1(:,2)+1,'.r');title(sprintf('plate id: %d, cam id %d',plate_id,cam_bundle{cam1,1}));plot(pix1(pointId,1)+1, pix1(pointId,2)+1,'xg');
                    subplot(1,2,2);imshow(cam_bundle{cam2,3});hold on;plot(pix2(:,1)+1, pix2(:,2)+1,'.r');plot(PT(:,1)+1,PT(:,2)+1,'.g');title(sprintf('plate id: %d, cam id %d',plate_id,cam_bundle{cam2,1}));                        
                    end

                end
                
            end
        end
    end
end




end
function [FrameCam, names] = readNames(fileName)


fid=fopen(fileName);       %首先打开文本文件coordinate.txt
temp = [];

names = {};
FrameCam = [];
while ~feof(fid)    % while循环表示文件指针没到达末尾，则继续
    % 每次读取一行, str是字符串格式
    str = fgetl(fid);
    idx = find(str == '/');
    %     names = [names; {str(end-21:end)}];
    names = [names; {str(idx(end-1)+1:end)}];
    FrameCam = [FrameCam;[str2num(str(1:5)) str2num(str(8))] ];
end
fclose(fid);





end
function [pt2d, inImageFlag] = projectKB8(pt3d, K, dis, R, t)
pt3d = ([R t; 0 0 0 1] * pt3d')';
pt3d = pt3d(:,1:3);

inImageFlag = pt3d(:,3) > 0.005;

fx = K(1,1);
fy = K(2,2);
cx = K(1,3);
cy = K(2,3);
k1 = dis(1);
k2 = dis(2);
k3 = dis(3);
k4 = dis(4);

x = pt3d(:,1);
y = pt3d(:,2);
z = pt3d(:,3);

r2 = x .* x + y .* y;
r = sqrt(r2);

theta = atan2(r, z);
theta2 = theta .* theta;
theta4 = theta2 .* theta2;
theta6 = theta4 .* theta2;
theta8 = theta6 .* theta2;

r_theta = theta .* (1 + k1 .* theta2 + k2 .* theta4 + k3 .* theta6 + k4 .* theta8);
if r > 1e-8
    norm_inv = 1./r;
else
    norm_inv = ones(length(r), 1);
end
% norm_inv = r > 1e-8 ? double(1.0) / r : 1;

mx = r_theta .* x .* norm_inv;
my = r_theta .* y .* norm_inv;

pt2d(:,1) = fx .* mx + cx;
pt2d(:,2) = fy .* my + cy;



end

function [theta,d_func_d_theta1] = solveTheta( k1,  k2,  k3,  k4,   r_theta ,d_func_d_theta, ITER)


theta = r_theta;
for  i = ITER:-1:1
    theta2 = theta * theta;
    
    func = k4 * theta2;
    func = func + k3;
    func = func * theta2;
    func = func + k2;
    func = func * theta2;
    func = func + k1;
    func = func * theta2;
    func = func + 1.0;
    func = func * theta;
    
    d_func_d_theta = (9) * k4 * theta2;
    d_func_d_theta = d_func_d_theta + (7) * k3;
    d_func_d_theta = d_func_d_theta * theta2;
    d_func_d_theta = d_func_d_theta + (5) * k2;
    d_func_d_theta = d_func_d_theta * theta2;
    d_func_d_theta = d_func_d_theta + (3) * k1;
    d_func_d_theta = d_func_d_theta * theta2;
    d_func_d_theta = d_func_d_theta + (1);
    
    
    %theta = theta + (r_theta - func) / d_func_d_theta;
    theta_fix = (r_theta - func) / d_func_d_theta;
    theta = theta+ theta_fix;
    if (abs(theta_fix) < 1e-8)
        break;
    end
end
d_func_d_theta1 = d_func_d_theta;
end
function p3d = unprojectKB8( fx,  fy,  cx,  cy,  k1,  k2,  k3,  k4, proj)


mx = (proj(1) - cx) / fx;
my = (proj(2) - cy) / fy;




theta = 0;
sin_theta = 0;
cos_theta = 1;
thetad = sqrt(mx * mx + my * my);

thetad = min(max(thetad,-3.141592653/2),3.141592653/2);
scaling = 1.0;
d_func_d_theta = 0;
if (thetad > 1e-8)
    [theta,d_func_d_theta1] = solveTheta(k1, k2, k3, k4, thetad, d_func_d_theta,3);
    d_func_d_theta = d_func_d_theta1;
    sin_theta = sin(theta);
    cos_theta = cos(theta);
    scaling = sin_theta / thetad;
end

p3d = [mx * scaling; my * scaling;cos_theta];


end
function [theta,d_func_d_theta1] = solveTheta_2( k1,  k2,  k3,  k4,   r_theta ,d_func_d_theta, ITER)


theta = r_theta;

for  i = ITER:-1:1
    theta2 = theta * theta;
    
    func = k4 * theta2;
    func = func + k3;
    func = func * theta2;
    func = func + k2;
    func = func * theta2;
    func = func + k1;
    func = func * theta2;
    func = func + 1.0;
    func = func * theta;
    
    d_func_d_theta = (9) * k4 * theta2;
    d_func_d_theta = d_func_d_theta + (7) * k3;
    d_func_d_theta = d_func_d_theta * theta2;
    d_func_d_theta = d_func_d_theta + (5) * k2;
    d_func_d_theta = d_func_d_theta * theta2;
    d_func_d_theta = d_func_d_theta + (3) * k1;
    d_func_d_theta = d_func_d_theta * theta2;
    d_func_d_theta = d_func_d_theta + (1);
    
    
    %theta = theta + (r_theta - func) / d_func_d_theta;
    theta_fix = (r_theta - func) / d_func_d_theta;
    theta = theta+ theta_fix;
    %     if ((theta_fix) < 1e-8)
    %         break;
    %     end
end
d_func_d_theta1 = d_func_d_theta;
end
function p3d = unprojectKB8_2( fx,  fy,  cx,  cy,  k1,  k2,  k3,  k4, proj)

mx = (proj(1) - cx) / fx;
my = (proj(2) - cy) / fy;

theta = 0;
sin_theta = 0;
cos_theta = 1;
thetad = sqrt(mx * mx + my * my);
scaling = 1.0;
d_func_d_theta = 0;

thetad = min(max(thetad,-3.141592653/2),3.141592653/2);

if (thetad > 1e-5)
    [theta,d_func_d_theta1] = solveTheta_2(k1, k2, k3, k4, thetad, d_func_d_theta,5);
    d_func_d_theta = d_func_d_theta1;
    sin_theta = sin(theta);
    cos_theta = cos(theta);
    scaling = sin_theta / thetad;
end
p3d = [mx * scaling; my * scaling;cos_theta];
end