function plateEx_ = showPath()

if 0
    plateEx_1 = showPath();
    plateEx_2 = showPath();
    plateEx_ = [plateEx_1 plateEx_2];
    diffVec = [];
    for j = 1 : size(plateEx_,3)
        x = 1000.*sort(plateEx_(10,:,j));
        y = 1000.*sort(plateEx_(11,:,j));
        z = 1000.*sort(plateEx_(12,:,j));
        diffVec = [diffVec [(x - x(1)) (y - y(1)) (z - z(1))]];
        figure(j),subplot(1,3,1);plot(x - x(1));
        subplot(1,3,2);plot(y - y(1));
        subplot(1,3,3);plot(z - z(1));
    end
    figure(100),hist(diffVec,50);title('plate extrinsic change distribution (mm)');
end






close all;
% inputDir = 'G:\matlab\data\dump\path0';
inputDir = 'G:\matlab\data\dump\path1';
inputDir = 'G:\matlab\data\dump\path2';% 不固定内参，做单目多板

inputDir = 'G:\matlab\data\dump\dump1';
inputDir = 'G:\matlab\data\dump\dump2'; % 固定内参，做单目多板

% inputDir = 'G:\matlab\data\dump\dump11';% 固定内参，不做单目多板
% inputDir = 'G:\matlab\data\dump\dump22';
%
%
%
inputDir = 'G:\matlab\data\dump\004\dump1'; % 不加weight，不固定内参，做单目多板
inputDir = 'G:\matlab\data\dump\004\dump2';
%
%
% inputDir = 'G:\matlab\data\dump\005\dump1'; % 加weight，不固定内参，做单目多板
% inputDir = 'G:\matlab\data\dump\005\dump2';
%
inputDir = 'G:\matlab\data\dump\006\dump1'; % 方向向量，加weight，不固定内参，做单目多板
inputDir = 'G:\matlab\data\dump\006\dump2';
% %
% inputDir = 'G:\matlab\data\dump\007\dump1'; % 方向向量，不加weight，不固定内参，做单目多板
% inputDir = 'G:\matlab\data\dump\007\dump2';
%
% inputDir = 'G:\matlab\data\dump\008\dump1'; % 方向向量，加weight(不平方)，不固定内参，做单目多板
% inputDir = 'G:\matlab\data\dump\008\dump2';
%
% inputDir = 'G:\matlab\data\dump\009\dump1'; %重投影，不加weight，固定内参，做单目多板
% inputDir = 'G:\matlab\data\dump\009\dump2';
%
% inputDir = 'G:\matlab\data\dump\010\dump1'; % 方向向量，加weight(不平方)，不固定内参，做单目多板
% inputDir = 'G:\matlab\data\dump\010\dump2';

inputDir = 'G:\matlab\data\dump\011\dump1'; % 方向向量，加weight(三次方)，不固定内参，做单目多板
inputDir = 'G:\matlab\data\dump\011\dump2';


inputDir = 'G:\matlab\data\dump\012\dump1'; % 方向向量，加weight(四次方)，不固定内参，做单目多板
inputDir = 'G:\matlab\data\dump\012\dump2';


% inputDir = 'G:\matlab\data\dump\013\dump1'; % 重投影误差，加weight(四次方)，不固定内参，做单目多板
% inputDir = 'G:\matlab\data\dump\013\dump2';


inputDir = 'G:\matlab\data\dump\014\dump1'; % 方向向量，不加权，开启内参优化，4edge
inputDir = 'G:\matlab\data\dump\014\dump2'; % 方向向量，不加权，开启内参优化，4edge
inputDir = 'G:\matlab\data\dump\014\dump11'; % 方向向量，不加权，开启内参优化，6edge
inputDir = 'G:\matlab\data\dump\014\dump22'; % 方向向量，不加权，开启内参优化，6edge



%20240307
inputDir = 'G:\matlab\data\dump\new_cam_model_001';


dirInfo = dir(fullfile(inputDir, 'path*.txt'));
dirInfo2 = dir(fullfile(inputDir, 'plate*.txt'));
try
    dirInfoCam = dir(fullfile(inputDir, 'cam*.txt'));
catch
    dirInfoCam = [];
end


dirInfo = dirInfo(4:6);
dirInfo2 = dirInfo2(4:6);
dirInfoCam = dirInfoCam(4:6);
% plateEx = [];
for i = 1 : length(dirInfo)
    
    data = load(fullfile(inputDir, dirInfo(i).name));
    data_ = load(fullfile(inputDir, dirInfo2(i).name));
    try
        data_cam = load(fullfile(inputDir, dirInfoCam(i).name));
    catch
        data_cam = [];
    end
    data2 = [];
    for j = 1 : size(data,1)
        T = [reshape(data(j,1:9),3,3), data(j,10:12)'; 0 0 0 1];
        T_inv = inv(T);
        data2 = [data2; [reshape(T_inv(1:3,1:3),1,9) T_inv(1:3,4)']];
    end
    figure,plotPath(data);
    %         data2 = [];
    %    for j = 1 : size(data_,1)
    %        T = [reshape(data_(j,1:9),3,3), data_(j,10:12)'; 0 0 0 1];
    % %        T_inv = inv(T);
    % %        data2 = [data2; [reshape(T_inv(1:3,1:3),1,9) T_inv(1:3,4)']];
    %    end
    plateEx(:,:,i) = data_;
    plotPath(data_);
    camEx(:,:,i) = data_cam;
end
plateEx_ = permute(plateEx, [2,3,1]);

diffVec = [];
for j = 1 : size(plateEx_,3)
    
    x = 1000.*sort(plateEx_(10,:,j));
    y = 1000.*sort(plateEx_(11,:,j));
    z = 1000.*sort(plateEx_(12,:,j));
    diffVec = [diffVec [(x - x(1)) (y - y(1)) (z - z(1))]];
    figure(j),subplot(1,3,1);plot(x - x(1));
    subplot(1,3,2);plot(y - y(1));
    subplot(1,3,3);plot(z - z(1));
end
figure(100),hist(diffVec(diffVec~=0),10);title('plate extrinsic change distribution (mm)');

try
    camEx_ = permute(camEx(:,1:12,:), [2,3,1]);
    camIntr_ = permute(camEx(:,13:end,:), [2,3,1]);
    diffVecCam = [];
    diffVecIntr = [];
    for j = 1 : size(camEx_,3)
        
        x = 1000.*sort(camEx_(10,:,j));
        y = 1000.*sort(camEx_(11,:,j));
        z = 1000.*sort(camEx_(12,:,j));
        diffVecCam = [diffVecCam [(x - x(1)) (y - y(1)) (z - z(1))]];
        figure(j),subplot(1,3,1);plot(x - x(1));
        subplot(1,3,2);plot(y - y(1));
        subplot(1,3,3);plot(z - z(1));
        fx = sort(camIntr_(1,:,j));
        fy = sort(camIntr_(2,:,j));
        cx = sort(camIntr_(3,:,j));
        cy = sort(camIntr_(4,:,j));
        diffVecIntr = [diffVecIntr [(fx - fx(1)) (fy - fy(1)) (cx - cx(1)) (cy - cy(1))]];
        figure(10+j),subplot(2,2,1);plot(fx - fx(1));title('fx');
        subplot(2,2,2);plot(fy - fy(1));title('fy');
        subplot(2,2,3);plot(cx - cx(1));title('cx');
        subplot(2,2,4);plot(cy - cy(1));title('cy');
    end
    figure(101),hist(diffVecCam(diffVecCam~=0),10);title('camera extrinsic change distribution (mm)');
    figure(102),hist(diffVecIntr(diffVecIntr~=0),10);title('camera intrinsic change distribution (pixel)');
catch
    sdfkgj = 1;
end

end