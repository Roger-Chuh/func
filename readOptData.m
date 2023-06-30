function readOptData()
inputDir = 'G:\matlab\data\direct\sim';

[Data, valid] = readData(fullfile(inputDir, 'opt_data1.txt'));
if ~valid
    [Data, valid] = readData(fullfile(inputDir, 'opt_data_inertial.txt'));
end
[info, inertial_mode] = arrangeData(Data);


Pose = {};
for frameId = 1 : size(info,2)
    Pose{frameId,1} = [];
    for iter =   1 :  size(info,1)
        pose = info{iter,frameId};
        Pose{frameId,1} = [Pose{frameId,1}; pose];
        
    end
    
end

if inertial_mode
    PoseMat = {};
    for i = 1 : length(Pose)
        pose_mat = [];
        for j = 1 : length(Pose{i})
            pose_mat = [pose_mat; Pose{i}(j).pose];
        end
        PoseMat{i,1} = pose_mat;
    end
    bgbagravity = [];
    for i =  1:length(Pose{1})
        bgbagravity = [bgbagravity; [Pose{1}(i).bg Pose{1}(i).ba Pose{1}(i).gravity]];
    end
    
    figure,plot(cell2mat(PoseMat))
    
    figure,subplot(1,3,1);plot(bgbagravity(:,1:3));grid on;title('bg');
    subplot(1,3,2),plot(bgbagravity(:,4:6));grid on;title('ba');
    subplot(1,3,3),plot(bgbagravity(:,7:9));grid on;title('gravity');
    asbj = 1;
else
        PoseMat = {};
        BiasMat = [];
    for i = 1 : length(Pose)
        pose_mat = [];
        bias_mat = [];
        for j = 1 : length(Pose{i})
            pose_mat = [pose_mat; Pose{i}(j).pose];
            bias_mat = [bias_mat; [Pose{i}(j).bg Pose{i}(j).ba]];
        end
        PoseMat{i,1}.pose = pose_mat;
        PoseMat{i,1}.bias = bias_mat;
        BiasMat = [BiasMat; bias_mat];
    end
   
figure,plot(PoseMat{4, 1}.bias(:,4:6))    
end




end

function [Data, valid] = readData(fileName)

fid=fopen(fileName);       %首先打开文本文件coordinate.txt

if fid == -1
    Data = 0;
    valid = false;
    return;
end
valid = true;
temp = [];

names = {};
FrameCam = [];
Data = {};
cnt = 1;
while ~feof(fid)    % while循环表示文件指针没到达末尾，则继续
    % 每次读取一行, str是字符串格式
    str = fgetl(fid);     
    a = strsplit(str,' ');
%     idx = find(str == '/');
% %     names = [names; {str(end-21:end)}];
%     names = [names; {str(idx(end-1)+1:end)}];
%    FrameCam = [FrameCam;[str2num(str(1:5)) str2num(str(8))] ];
data.pose = [];
data.bg = [];
data.ba = [];
data.gravity = [];
data.iter = str2double(a(2));
data.fid = str2double(a{4});
if (str2double(a{4}) == -100)
    data.is_inertial_mode = 1;
else
    data.is_inertial_mode = 0;
end
if strcmp(a{5}, 'pose_vec')
    data.pose = str2double(a(6:11));
end
if data.is_inertial_mode
    if strcmp(a{5}, 'bg')
        data.bg = str2double(a(6:8));
    end
    if strcmp(a{9}, 'ba')
        data.ba = str2double(a(10:12));
    end
    if strcmp(a{13}, 'gravity')
        data.gravity = str2double(a(14:16));
    end
else
    if strcmp(a{5}, 'bg')
        data.bg = str2double(a(6:8));
    end
    if strcmp(a{9}, 'ba')
        data.ba = str2double(a(10:12));
    end
end
Data{cnt,1} = data;
cnt = cnt + 1;
end
fclose(fid);



end
function [info, inertial_mode] = arrangeData(Data)

info = {};
for i = 1 : length(Data)
    
    dlt = 0;
    if ~isempty(Data{i,1}.pose)
        info{Data{i,1}.iter+1,Data{i,1}.fid+dlt}.pose = Data{i,1}.pose;
    end
    if Data{i,1}.is_inertial_mode == 1
        inertial_mode = 1;
        info{Data{i,1}.iter+1,1}.bg = Data{i,1}.bg;
        info{Data{i,1}.iter+1,1}.ba = Data{i,1}.ba;
        info{Data{i,1}.iter+1,1}.gravity = Data{i,1}.gravity;
    else
        inertial_mode = 0;
        info{Data{i,1}.iter+1,Data{i,1}.fid+dlt}.bg = Data{i,1}.bg;
        info{Data{i,1}.iter+1,Data{i,1}.fid+dlt}.ba = Data{i,1}.ba;
    end
end


end