function readSimData2()
% close all
inputDir = 'G:\matlab\data\sim\keyframe5';
inputDir = 'G:\matlab\data\sim\traj1';
inputDir = 'G:\matlab\data\gt_traj';
% dirInfo = dir(fullfile(inputDir, 'all_points_*.txt'));
K = eye(3);
K(1,1) = 236.8820281;
K(2,2) = 236.7437341;
K(1,3) = 325.0065396;
K(2,3) = 241.9557923;
dis = [0.2124214875, -0.1719834148, 0.06211290885, -0.009899091649]';


R =[ 0, 0, -1; -1, 0, 0; 0, 1, 0];
R_bc = R;
t_bc = [0.05, 0.04, 0.03];
% t_bc = [0.0, 0.0, 0.0];

R = rodrigues([0.1 0.2 0.3]);
q = rotMat2quatern(R);
R2 = quatern2rotMat(q);
err = R2-R;

% for i = 1 : length(dirInfo)
%     
%     data = load(fullfile(inputDir, dirInfo(i).name));
%     
%     
%     %figure(1);subplot(1,2,1);pcshow(data(:,1:3),'MarkerSize',100);
%     %     figure(1);imshow(zeros(480, 640));hold on;plot(data(:,5), data(:,6),'or');
%     %     drawnow;
%     
% end


% all_points = load(fullfile(inputDir,'all_points.txt')); % golden pose
% cam_pose = load(fullfile(inputDir,'cam_pose.txt')); % golden pose
imu_pose = load(fullfile(inputDir,'trajectory_cam0.txt')); % golden pose
% imu_int_pose = load(fullfile(inputDir,'imu_int_pose.txt')); % integrated pose
% camPoseMat = [];
% for i = 1 : size(imu_pose,1)
%     timestamp = cam_pose(i,1);
%     data = cam_pose(i,2:8);
%     % cam in world
%     R = quatern2rotMat(data(1:4));
%     t = data(5:7);
%     T_cw = inv([R' t';0 0 0 1]);
%     T_wc = inv(T_cw);
%    
% end

imuPoseMat = [];
for i = 1 : size(imu_pose,1)
    data = imu_pose(i,2:8);
    R = quatern2rotMat(data([7 4 5 6]));
%     R = quatern2rotMat(data([4 : 7]));
    t = data(1:3);
    T_wc = ([R t';0 0 0 1]);
%     T_wc = inv(T_wc);
if i == 1
    T_base = T_wc;
    T_wc = eye(4);
else
    T_wc = inv(T_base) * T_wc;
end
    imuPoseMat = [imuPoseMat; [reshape(T_wc(1:3,1:3),1,9) T_wc(1:3,4)']];
end
% figure,plotPath(camPoseMat(1:10:end,:))
figure,plotPath(imuPoseMat(1:20:end,:))
end
function [pt2d, inImageFlag] = projectKB8(pt3d, K, dis, R, t)
pt3d = ([R t; 0 0 0 1] * pt3d')';
pt3d = pt3d(:,1:3);

inImageFlag = pt3d(:,3) > 5;

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