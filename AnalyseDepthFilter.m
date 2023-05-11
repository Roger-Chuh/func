function AnalyseDepthFilter()

close all;

inputDir = 'C:\Users\Roger\Desktop\vins-fusion-fisheye\yvr\april\controller_2022-04-12-05-46-16';
inputDir = 'G:\matlab\data\tag\controller_2022-05-31-06-33-32';
inputDir = 'G:\matlab\data\direct\gt\1';
scale = 1;
camInfo = dir(fullfile(inputDir, 'Camera*'));

depth_filter_data = load(fullfile(inputDir, 'depth_filter_statistics.txt'));

depth_filter_data(:,4:12) = depth_filter_data(:,4:12) + 1;


dirCam0 = dir(fullfile(inputDir, camInfo(1).name,'images','*.bmp'));
dirCam3 = dir(fullfile(inputDir, camInfo(4).name,'images','*.bmp'));
assert(length(dirCam0) == length(dirCam3));
timestamp = zeros(length(dirCam0),1);
for i = 1 : length(dirCam0)
    timestamp(i,1) = str2double(dirCam0(i).name(1:end-4));
end

frameNum = unique(depth_filter_data(:,1));

featIndAll = find(abs(depth_filter_data(:,1) - timestamp(1)) < 1000);
featIdAll = depth_filter_data(featIndAll, 4);
% hostCoord = depth_filter_data(featIndAll, 5:6);
% hostCoord = zeros
featMatX = zeros(max(depth_filter_data(featIndAll,4)), length(frameNum));
featMatY = zeros(max(depth_filter_data(featIndAll,4)), length(frameNum));
featMatAX = zeros(max(depth_filter_data(featIndAll,4)), length(frameNum));
featMatBX = zeros(max(depth_filter_data(featIndAll,4)), length(frameNum));
% featMatX(featIdAll,1) = hostCoord(:,1);
% featMatY(featIdAll,1) = hostCoord(:,2);
% featMatAX(featIdAll,1) = hostCoord(:,1);
% featMatAY(featIdAll,1) = hostCoord(:,2);
% featMatBX(featIdAll,1) = hostCoord(:,1);
% featMatBY(featIdAll,1) = hostCoord(:,2);

for imgId = 1 : length(frameNum) %length(dirCam0)
%     img0 = imread(fullfile(inputDir, camInfo(1).name, 'images',dirCam0(imgId).name));
%     img3 = imread(fullfile(inputDir, camInfo(4).name, 'images',dirCam3(imgId).name));
    feats = find(abs(depth_filter_data(:,1) - timestamp(imgId)) < 1000);
    featid = depth_filter_data(feats, 4);
    hostCoord(featid,:) = depth_filter_data(feats, 5:6);
    px_cur = depth_filter_data(feats, 11:12);
    px_A = depth_filter_data(feats, 7:8);
    px_B = depth_filter_data(feats, 9:10);
    zList = depth_filter_data(feats, 13);
    featMatX(featid,imgId) = px_cur(:,1);
    featMatY(featid,imgId) = px_cur(:,2);
    featMatAX(featid,imgId) = px_A(:,1);
    featMatAY(featid,imgId) = px_A(:,2);
    featMatBX(featid,imgId) = px_B(:,1);
    featMatBY(featid,imgId) = px_B(:,2);
    zMat(featid,imgId) = zList;
    %     for feat_id = 1 : feats
    %     end
    
    
    
end


validId = find(sum(featMatX > 0,2) > 0);
img0 = imread(fullfile(inputDir, camInfo(1).name, 'images',dirCam0(1).name));
img3 = imread(fullfile(inputDir, camInfo(4).name, 'images',dirCam3(1).name));
figure(1),subplot(2,2,1);imshow(img0);hold on;subplot(2,2,2);
for feat_id = (validId')
    
    figure(1),subplot(2,2,1);cla;imshow(img0);plot(hostCoord(feat_id,1), hostCoord(feat_id,2),'or', 'MarkerSize',3,'LineWidth',3);
    figure(1),subplot(2,2,2);cla;imshow(img3);hold on;
    figure(1),subplot(2,2,[3 4]);cla;
    z_buf = [];
    for imgid = 1 : size(featMatX,2)
        if featMatX(feat_id, imgid)~=0
            z_buf = [z_buf; zMat(feat_id, imgid)];
            figure(1),subplot(2,2,2);plot(featMatAX(feat_id,imgid), featMatAY(feat_id,imgid),'og','MarkerSize',3,'LineWidth',3');
            plot(featMatBX(feat_id,imgid), featMatBY(feat_id,imgid),'ob','MarkerSize',3,'LineWidth',3');
            plot(featMatX(feat_id,imgid), featMatY(feat_id,imgid),'or','MarkerSize',3,'LineWidth',3');title(sprintf('trial id: %d\n', imgid));
            figure(1),subplot(2,2,[3 4]);plot(z_buf);drawnow;
            ashk = 1;
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