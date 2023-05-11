close all

intrMat = [  5.19278503e+02, 0., 6.64875183e+02; 0., 5.19322693e+02,       3.93181030e+02; 0., 0., 1.];



a = load('D:\Temp\20201218\map\lab.txt');
b = load('D:\Temp\20201218\map\lab_ori.txt');

a = load('D:\Temp\20201218\map2\bk\bk\lab.txt');
b = load('D:\Temp\20201218\map2\bk\bk\lab2.txt');


a = load('D:\Temp\20201218\map2\bk\bk\lab2.txt');
b = load('D:\Temp\20201218\map2\bk\bk\lab3.txt');


baseId = [7 48];

use7 = 0;1;

if use7
    a = load('D:\Work\marker_mapper1.0.15\20201224_17959\result7\lab.yml');
    b = load('D:\Work\marker_mapper1.0.15\20201224_17959\result7\lab.yml');
    pose = load('D:\Work\marker_mapper1.0.15\20201224_17959\result7\lab.log');
    points = loadpcd('D:\Work\marker_mapper1.0.15\20201224_17959\result7\lab.pcd');
    
%     a = load('D:\Work\marker_mapper1.0.15\20201224_17959\result35\lab.yml');
%     b = load('D:\Work\marker_mapper1.0.15\20201224_17959\result35\lab.yml');
%     pose = load('D:\Work\marker_mapper1.0.15\20201224_17959\result35\lab.log');
%     points = loadpcd('D:\Work\marker_mapper1.0.15\20201224_17959\result35\lab.pcd');
    
else
    a = load('D:\Work\marker_mapper1.0.15\20201224_17959\result45\lab.yml');
    b = load('D:\Work\marker_mapper1.0.15\20201224_17959\result45\lab.yml');
    pose = load('D:\Work\marker_mapper1.0.15\20201224_17959\result45\lab.log');
    points = loadpcd('D:\Work\marker_mapper1.0.15\20201224_17959\result45\lab.pcd');
end

xyz1 = [];
validId = [];
idRange = [39:50];
planeNorm = [];
for i = 1 : size(a,1)
    
    marker3d = 1000.*reshape(a(i,2:13),3,[])';
    
    if ismember(a(i,1), idRange)
        validId = [validId; [(4*i-3) : 4*i]'];
        B = fitplane(marker3d');B = B./norm(B(1:3));
        if use7
            planeNorm = [planeNorm sign(B(1)).*B(1:3)];
        else
            planeNorm = [planeNorm sign(B(3)).*B(1:3)];
        end
    end
    
    if a(i,1) == baseId(1)
        
        oldBase = marker3d;
        
    end
    
    if a(i,1) == baseId(2)
        
        newBase = marker3d;
        
    end
    xyz1 = [xyz1;marker3d];
    
end
xyz1_1 = xyz1;
if 0
    rt1 = AlignCoord(oldBase,newBase);
    
    xyzExp = (rt1(1:3,1:3)*newBase' + repmat(rt1(1:3,4), 1, size(newBase,1)))';
    
    xyz1 = (rt1(1:3,1:3)*xyz1_1' + repmat(rt1(1:3,4), 1, size(xyz1_1,1)))';
    
    
    xyz1_2 = xyz1;
    
    coordShiftMat = eye(3);  roty(-90)*rotx(90);
    
    xyz1 = (coordShiftMat*xyz1_2')';
else
    coordShiftMat = eye(3);  
end

xyzPlane = xyz1(validId,:);

% figure,pcshow(xyzPlane, 'MarkerSize', 100)

ang2 = real(CalcDegree2(repmat(planeNorm(:,1)',size(planeNorm,2),1),planeNorm'));
figure,subplot(1,2,1);plotQuiver(planeNorm');subplot(1,2,2);hist(ang2,10);grid on;


poseMat = [];
for i = 1 : size(pose,1)
    t = 1000.*pose(i,2:4);
    r = quatern2rotMat(pose(i,[8 5 6 7]))';
    R(:,:,i) = r;
    RR(:,:,i) = r';
    shiftedPose0 = [coordShiftMat [0 0 0]';0 0 0 1]*[r t'; 0 0 0 1];
%     shiftedPose = [roty(90)*rotx(90) [0 0 0]';0 0 0 1]*shiftedPose0;
    shiftedPose = [eye(3) [0 0 0]';0 0 0 1]*shiftedPose0;
    poseMat = [poseMat; [reshape(shiftedPose(1:3,1:3),1,9) shiftedPose(1:3,4)']];
end
figure,plotPath(poseMat)


% BB = fitplane(poseMat(:,10:12)');BB = BB./norm(BB(1:3));
[~, ~, BB] = robustfitXYZ(poseMat(:,10:12),[]);
planeErr1 = dot(pextend(poseMat(:,10:12)'), repmat(BB, 1, size(poseMat,1)));

% CC = fitplane(xyzPlane');CC = CC./norm(CC(1:3));
[~, ~, CC] = robustfitXYZ(xyzPlane,[]);
planeErr2 = dot(pextend(xyzPlane'), repmat(CC, 1, size(xyzPlane,1)));

figure,subplot(1,2,1);hist(planeErr1, 100);grid on;title('pose plane');subplot(1,2,2),hist(planeErr2, 20);grid on;title('ground plane');


figure,plotPath(poseMat)


figure,subplot(1,3,1);pcshow(points(1:3,:)', 'MarkerSize',100);
subplot(1,3,2);pcshow(xyzPlane, 'MarkerSize',100);
subplot(1,3,3);pcshow(poseMat(:,10:12), 'MarkerSize',100);
% subplot(2,2,4);plotPath(poseMat)

planeNormAng = CalcDegree(BB(1:3)', CC(1:3)')
[slamPoseMat2, b2c] = rebaseStartRot(poseMat, rotz(90)*roty(90), xyzPlane, intrMat);


MarkerPoseData = load('D:\Work\marker_mapper1.0.15\MarkerPose\5\MarkerPose.txt');
c2b = load('D:\Work\marker_mapper1.0.15\c2b\1\c2b.txt');
plotMarkerPose(MarkerPoseData, inv(c2b));










frameId = 1;
pose1 = inv([reshape(poseMat(frameId,1:9),3,3) poseMat(frameId,10:12)';0 0 0 1]);
% transPt = EucildTransform(xyzPlane, pose1(1:3,1:3), pose1(1:3,4));
transPt = EucildTransform(xyz1, pose1(1:3,1:3), pose1(1:3,4));
figure,pcshow(transPt,'MarkerSize', 100)

slamPoseMat2 = rebaseStartRot(poseMat, rotz(90)*roty(90));

xyz2 = [];
for i = 1 : size(b,1)
    xyz2 = [xyz2;1000.*reshape(b(i,2:13),3,[])'];
    
end


rt = AlignCoord(xyz1,xyz2);
rt1 = AlignCoord(xyz2(21:24,:),xyz2(33:36,:));

xyz22 = (rt(1:3,1:3)*xyz2' + repmat(rt(1:3,4), 1, size(xyz2,1)))';




fid1 = fopen('base_10.txt','w');
fprintf(fid1,'%%');
fprintf(fid1,sprintf('YAML:1.0\n---\naruco_bc_dict: ARUCO_MIP_36h12\naruco_bc_nmarkers: %d\naruco_bc_mInfoType: 1\naruco_bc_markers:\n', size(b,1)));





for j = 1 : size(b,1)
    xyzMat = xyz22(4*j-3:4*j,:)'./1000;
    
    fprintf(fid1,'   - { id:%d, corners:[ [%0.8f, %0.8f, %0.8f],\n                       [%0.8f, %0.8f, %0.8f],\n                       [%0.8f, %0.8f, %0.8f],\n                       [%0.8f, %0.8f, %0.8f] ]}\n',b(j,1), xyzMat(:));
   
%     bb(j,:) = [b(j,1) xyzMat(:)'./1000];
    
end
fclose(fid1);