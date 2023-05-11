function showAprilTag()
inputDir = 'C:\Users\Roger\Desktop\vins-fusion-fisheye\yvr\april\controller_2022-04-12-05-46-16';
inputDir = 'G:\matlab\data\tag\controller_2022-05-31-06-33-32';
scale = 1;
camInfo = dir(fullfile(inputDir, 'Camera*'));

% if scale == 2
    data2 = load(fullfile(inputDir, 'results','outlier_img','matlab_apriltag_corners.txt'));
    name2 = (fullfile(inputDir, 'results','outlier_img','matlab_image_names.txt'));
% else
    data1 = load(fullfile(inputDir, 'results','outlier_img1','matlab_apriltag_corners.txt'));
    name1 = (fullfile(inputDir, 'results','outlier_img1','matlab_image_names.txt'));
% end
[FrameCam1, names1] = readNames(name1);
[FrameCam2, names2] = readNames(name2);
frameList = intersect(unique(FrameCam1(:,1)), unique(FrameCam2(:,1)));
coviCamsEachTimestamp = {};

for i = 1 : length(frameList)
    frameId = frameList(i);
    frameData1 = find(data1(:,1) == frameId);
    frameData2 = find(data2(:,1) == frameId);
    camList = intersect(unique(data1(frameData1,2)),unique(data2(frameData2,2)));
    coviCamsEachTimestamp{i,1} = i;
    coviCamsEachTimestamp{i,2} = camList;
    close all;
    for j = 1 : length(camList)
        camId = camList(j);
        frameCamData1 = find(data1(:,1) == frameId & data1(:,2) == camId);
        frameCamData2 = find(data2(:,1) == frameId & data2(:,2) == camId);
        frameCamImageData = find(FrameCam1(:,1) == frameId & FrameCam1(:,2) == camId);
        imageFileName = names1{frameCamImageData};
        img1 = imread(fullfile(inputDir, camInfo(camId+1).name,imageFileName));
        img2 = imresize(imread(fullfile(inputDir, camInfo(camId+1).name,imageFileName)),2);
        pixels1 = data1(frameCamData1,4:5);
        pixels2 = data2(frameCamData2,4:5);
        patternIdList1 = unique(data1(frameCamData1,3))';
        patternIdList2 = unique(data2(frameCamData1,3));
        camIdList1 = unique(data1(frameCamData1,2));
        camIdList2 = unique(data2(frameCamData1,2));
        pixels21 = [(pixels2(:,1)+1.5)/2,(pixels2(:,2)+1.5)/2];
        [D] = pdist2(pixels1+1,pixels21,'euclidean');
        idx = find(D < 2);
        [y,x] = ind2sub(size(D),idx);
        
        match = [pixels1(y,:)+1 pixels21(x,:)];
%         figure(2),clf;quiver(match(:,1),match(:,2),match(:,1)-match(:,3),match(:,2)-match(:,4));axis equal
        figure(j),clf;
        subplot(1,2,1);imshow(img1);hold on;plot(pixels1(:,1)+1,pixels1(:,2)+1,'.r');plot((pixels2(:,1)+1.5)/2,(pixels2(:,2)+1.5)/2,'.g');title(sprintf('frameid: %d, camid: %d, 640 x 480',frameId, camIdList1));legend('640 x 480', '1280 -> 640');drawnow;
        subplot(1,2,2);imshow(img2);hold on;plot(pixels2(:,1)+1,pixels2(:,2)+1,'.r');plot((pixels1(:,1)+1)*2-0.5,(pixels1(:,2)+1)*2-0.5,'.g');title(strcat('1280 x 960 ,pattern list: ',num2str(patternIdList1)));legend('1280 x 960', '640 -> 1280');drawnow;
%         
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