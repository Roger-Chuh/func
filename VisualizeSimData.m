function VisualizeSimData()
inputDir = 'G:\matlab\data\direct\sim';
Data = load(fullfile(inputDir, 'sim_data.txt'));

delta = 1;

Data = [Data(:,1+delta:end) Data(:,1:delta)];

%% fid  |  trans  |   qxyzw   |     trans   |   qxyzw        |  cid  |  pid  |  landmark  |  bearing  
%%  1   |  2 3 4  |  5 6 7 8  |   9 10 11   |  12 13 14 15   |   16  |  17   |  18 19 20  |  21 22 23


[comm, id1, id2] = intersect(Data(:,1), Data(:,1));
camId = Data(id1, 16);
[comm_landmark, id_landmark1, id_landmark2] = intersect(Data(:,18:20), Data(:,18:20), 'rows');
pids = Data(id_landmark1, 17);
% assert(unique(diff(pids)) == 1);

poseMat = [];
poseMatCam = [];
CamTrace = cell(4,1);
cnt = 1;
for i = 1 :length(id1)
    
    
    data_list = Data(find(Data(:,1) == Data(id1(i),1)),:);
    camera_list = unique(data_list(:,16));
    
    for j = 1 : length(camera_list)
       camera_data = data_list(find(data_list(:,16) == camera_list(j)),:); 
       
        xyzw_cam1 = camera_data(1,12:15);
        trans_cam1 = camera_data(1,9:11);
        R_cam1 = quat2rotm(xyzw_cam1([4 1 2 3]));
        Tcw = inv([R_cam1 trans_cam1';0 0 0 1]);
        
        for(ii = 1 : size(camera_data,1))
           xyz_host = camera_data(ii, 21:23)';
           try
               xyz_check = Tcw(1:3,1:3)*comm_landmark(camera_data(ii,17)+1,:)' + Tcw(1:3,4);
           catch
               continue;
           end
           xyz_check = xyz_check./norm(xyz_check);
           err(cnt,1) = norm(xyz_host - xyz_check);
           cnt = cnt + 1;
        end
        
       camera_data(:,21) = camera_data(:,21)./camera_data(:,23);
       camera_data(:,22) = camera_data(:,22)./camera_data(:,23);
       camera_data(:,23) = 1;
       CamTrace{camera_list(j)+1} = [CamTrace{camera_list(j)+1}; {data_list(1,1), [camera_data]}];
    end
    
    data = Data(id1(i),:);
    xyzw = data(5:8);
    trans = data(2:4);
    
%     R2 = quatern2rotMat(xyzw([4 1 2 3]));
    R = quat2rotm(xyzw([4 1 2 3]));
    poseMat = [poseMat; [reshape(R,1,9), trans]];
    
    
    if (Data(id1(i),16) == 0)
        xyzw_cam = data(12:15);
        trans_cam = data(9:11);
        
        %     R2 = quatern2rotMat(xyzw([4 1 2 3]));
        R_cam = quat2rotm(xyzw_cam([4 1 2 3]));
        poseMatCam = [poseMatCam;  [reshape(R_cam,1,9), trans_cam]];
    end
    
    
%     rotm2quat(R)
%     rotMat2quatern(R2)
end

if 0
    figure,plotPath(poseMat(1:1:end,:));hold on; pcshow(Data(:,18:20),[1 0 0],'MarkerSize',10); axis equal;
    figure,plotPath(poseMatCam(1:1:end,:));hold on; pcshow(Data(:,18:20),[1 0 0],'MarkerSize',10); axis equal;
end


assert(size(CamTrace{1,1},1) == size(CamTrace{2,1},1));
assert(size(CamTrace{1,1},1) == size(CamTrace{3,1},1));
assert(size(CamTrace{1,1},1) == size(CamTrace{4,1},1));

if 1
    figure(100),
    for k = 1 : size(CamTrace{1,1},1)
        
        subplot(2,2,1),plot(CamTrace{2,1}{k,2}(:,21), CamTrace{2,1}{k,2}(:,22),'or','MarkerSize',3, 'LineWidth', 3); title(sprintf('cam1, frame: %d',k)); axis([-2 2 -2 2]);
        subplot(2,2,2),plot(CamTrace{3,1}{k,2}(:,21), CamTrace{3,1}{k,2}(:,22),'or','MarkerSize',3, 'LineWidth', 3); title('cam2'); axis([-2 2 -2 2]);
        subplot(2,2,3),plot(CamTrace{1,1}{k,2}(:,21), CamTrace{1,1}{k,2}(:,22),'or','MarkerSize',3, 'LineWidth', 3); title('cam0'); axis([-2 2 -2 2]);
        subplot(2,2,4),plot(CamTrace{4,1}{k,2}(:,21), CamTrace{4,1}{k,2}(:,22),'or','MarkerSize',3, 'LineWidth', 3); title('cam3'); axis([-2 2 -2 2]);
        drawnow;
        
    end
end

imuData = load(fullfile(inputDir,'imu_data.txt'));
imuPoseMat = [];
g = [-9.7964, 0, 0]';
for i = 1 : size(imuData,1)
    imu_data = imuData(i,:);
    xyzw = imu_data(11:14);
    trans = imu_data(8:10);
    
%     R2 = quatern2rotMat(xyzw([4 1 2 3]));
    R = quat2rotm(xyzw([4 1 2 3]));
    imuPoseMat = [imuPoseMat; [reshape(R,1,9), trans]];
end

% figure,plotPath(imuPoseMat(1:100:end,:));
[imuPoseMatInt, G, pose] = integrateImu(imuData(:,5:7), imuData(:,2:4), g, eye(3), eye(3), zeros(3,1), zeros(3,1),imuData(1,15:17)',reshape(imuPoseMat(1,1:9),3,3), imuPoseMat(1,10:12)', imuData(:,1).*1e-9);


ang_trans_error = [];
for i = 1 : size(imuData,1)
    angErr = rad2deg(norm(rodrigues(reshape(imuPoseMatInt(i,1:9),3,3)' * reshape(imuPoseMat(i,1:9),3,3))));
    transErr = norm(imuPoseMatInt(i,10:12) - imuPoseMat(i,10:12));
    ang_trans_error = [ang_trans_error; [angErr transErr]];
end

end