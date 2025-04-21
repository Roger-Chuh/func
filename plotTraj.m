function [poseMat, time] = plotTraj(filename)
inputFile = filename;

% inputFile = 'G:\matlab\data\vio_pose\pose_old1.txt';
% inputFile = 'G:\matlab\data\vio_pose\pose_new1.txt';
Data = load(inputFile);

poseMat = [];
time = [];
for i = 1 : size(Data,1)
    data = Data(i,:);
    xyzw = data(5:8);
    trans = data(2:4);
    if length(data) == 11
        vel = data(9:11);
    else
        vel = [];
    end
%     R2 = quatern2rotMat(xyzw([4 1 2 3]));
    R = quat2rotm(xyzw([4 1 2 3]));
    poseMat = [poseMat; [reshape(R,1,9), trans vel]];
%     rotm2quat(R)
%     rotMat2quatern(R2)
   time = [time; data(1)];
end

figure,plotPath(poseMat);
end