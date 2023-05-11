function [poseMat, G, pose,vel1,dis1] = integrateImu(acc, gyro, g0, scaleAcc, scaleGyro, biasAcc, biasGyro,vel0,dis0, dt,R)
% vel0 = [0 0 0];
% dis0 = [0 0 0]';
% R = eye(3);


poseMat = [reshape(R,1,9) dis0'];

aRot = [];
aRot1 = [];
eulerWorldSingle = [];
AWorld = [];
Angle = [];
Ang = [];
G = [];
aOrig = [];
Vel = [];
DisPlace = [];
ROT = [];
eulerWorldAcumu = [];

for i = 2 : size(acc,1) - 0
    if i == 43
        avi = 1;
    end
    
    [vel1, rot1, t1,angle,ang,dis1,aWorld] = IntgAcc3(scaleAcc*acc(i-1,:)' + biasAcc - g0, [] ,scaleGyro*gyro(i-1,:)' + biasGyro, [], vel0, dis0, dt,R);
    sdbjk = 1;
    
%     [vel1, rot1, t1,angle,ang,dis1,aWorld] = IntgAcc3((a(i-1,:)').*sf + ba - g0, (a(i,:)').*sf + ba  - g0 ,w(i-1,:)' + bw, w(i,:)' + bw, vel0, dis0, dt,R);
    
    
    
    % [vel1, rot1, t1,angle,ang,dis1,aWorld] = IntgAcc3((a(i-1,:)').*sf + ba - g0, (a(i,:)').*sf + ba  - g0 ,[0;0;0], [0;0;0], vel0, dis0, dt,R);
    % [vel1, rot1, t1,angle,ang,dis1,aWorld] = IntgAcc3((a(startt2,:)').*sf + ba - ((a(startt2,:)').*sf + ba), (a(startt2+1,:)').*sf + ba  - g0 ,[0;0;0], [0;0;0], vel0, dis0, dt,R);
    % [vel1, rot1, t1,angle,ang,dis1,aWorld] = IntgAcc3((a(startt2,:)').*sf + ba - ((a(startt2,:)').*sf + ba), (a(startt2+1,:)').*sf + ba  - g0 ,w(i-1,:)' + bw, w(i,:)' + bw, vel0, dis0, dt,R);
    
    eulerWorldSingle = [eulerWorldSingle; rodrigues(rot1)'];
    
    
    AWorld = [AWorld; aWorld'];
    
    
    aRot = [aRot; (R*(scaleAcc*(acc(i-1,:)') + biasAcc - g0))'];
%     aRot1 = [aRot1; (R*(scaleAcc(acc(i-1,:)') + biasAcc) - (sf.*(mean(a(startt1:startt2,:),1)') + ba))'];
    
% %     aRot = [aRot; (R*((acc(i-1,:)').*sf + ba - g0))'];
% %     aRot1 = [aRot1; (R*((acc(i-1,:)').*sf + ba) - (sf.*(mean(a(startt1:startt2,:),1)') + ba))'];
    
    
    
    Angle = [Angle;angle];
    Ang = [Ang;ang];
    g0 = rot1'*g0;  % g from world to local
    G = [G; g0'];
    
    dis0 = dis1; % world position
        
    aOrig = [aOrig; (R*(scaleAcc*(acc(i,:)') + biasAcc))'];
    
    Vel = [Vel; vel1]; % world velocity
    
    
    DisPlace = [DisPlace; dis1']; % pose in world

    
    vel0 = vel1;

    
    R = R*rot1;
    % pose = [pose; Rt(:)'];
    ROT = [ROT; rodrigues(R)'];  % rotation from local to world
    eulerWorldAcumu = [eulerWorldAcumu; rodrigues(R)'];
    % % pose = [rodrigues(R)' dis1'];
    
     
    
    poseMat = [poseMat; [reshape(R,1,9) dis1']];
    
    
    
end


pose = [ROT DisPlace Vel];







end