function plotPose()

inputFile = 'G:\matlab\data\vio_pose\imupose.txt';

% inputFile = 'G:\matlab\data\vio_pose\pose_old1.txt';
% inputFile = 'G:\matlab\data\vio_pose\pose_new1.txt';
Data = load(inputFile);

poseMat = [];
for i = 1 : size(Data,1)
    data = Data(i,:);
    xyzw = data(1:4);
    trans = data(5:7);
    
%     R2 = quatern2rotMat(xyzw([4 1 2 3]));
    R = quat2rotm(xyzw([4 1 2 3]));
    poseMat = [poseMat; [reshape(R,1,9), trans]];
%     rotm2quat(R)
%     rotMat2quatern(R2)
end











inputFile = 'G:\matlab\data\vio_pose\output.txt';
inputFile = 'G:\matlab\data\vio_pose\output_ref.txt';
% inputFile = 'G:\matlab\data\vio_pose\output_ref_cam0.txt';
% inputFile = 'G:\matlab\data\vio_pose\output_direct_cam0.txt';
Data = load(inputFile);

poseMat = [];
T = {};
time = [];
for i = 1 : size(Data,1)
    data = Data(i,:);
    time = [time; data(1)];
    trans = data(2:4);
    xyzw = data(5:8);
%     R = quatern2rotMat(xyzw([4 1 2 3]));
    R2 = quat2rotm(xyzw([4 1 2 3]));
    T{i,1} = [R2 trans';0 0 0 1];
    poseMat = [poseMat; [reshape(R2,1,9), trans]];
%     rotm2quat(R2)
%     rotMat2quatern(R)
end


end