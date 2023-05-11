function [slamPoseMat22, b2c] = rebaseStartRot(slamPoseMat, extraRotMat, xyzPlane, intrMat)

baseNum = 1;
slamPoseMat2 = [];poseFidd = [];
rt1675 = [reshape(slamPoseMat(baseNum,1:9),3,3) slamPoseMat(baseNum,10:12)';0 0 0 1];
rt1675_ = [extraRotMat [0 0 0]';0 0 0 1]*rt1675;
for ii = baseNum:size(slamPoseMat,1)
    tmpRT =  [reshape(slamPoseMat(ii,1:9),3,3) slamPoseMat(ii,10:12)';0 0 0 1];
    tmpRt2 = inv(rt1675)*tmpRT;
    slamPoseMat2 = [slamPoseMat2;[reshape(tmpRt2(1:3,1:3),1,9) tmpRt2(1:3,4)']];
    poseFidd = [poseFidd;[tmpRt2(1:3,4)' tmpRt2(3,1:3)]];
end



% slamPoseMat22 = [];
% for i = 1 : size(slamPoseMat2,1)
%     tmpRT_ =  [reshape(slamPoseMat2(i,1:9),3,3) slamPoseMat2(i,10:12)';0 0 0 1];
%     
% end


if 0
    slamPoseMat22 = [reshape(eye(3),1,9) [0 0 0]];
    
    for i = 1 : size(slamPoseMat,1) - 1
        Tp =  [reshape(slamPoseMat(i,1:9),3,3) slamPoseMat(i,10:12)';0 0 0 1];
        Tc =  [reshape(slamPoseMat(i+1,1:9),3,3) slamPoseMat(i+1,10:12)';0 0 0 1];
        p2c = inv(Tc) * Tp;
        
        slamPoseMatLast =  [reshape(slamPoseMat22(end,1:9),3,3) slamPoseMat22(end,10:12)';0 0 0 1];
        slamPoseMatNew = slamPoseMatLast * inv(p2c);
        slamPoseMat22 = [slamPoseMat22; [reshape(slamPoseMatNew(1:3,1:3),1,9) slamPoseMatNew(1:3,4)']];
    end
else
    slamPoseMat22 = []; % [reshape(eye(3),1,9) [0 0 0]];
    Tbase =  [reshape(slamPoseMat(1,1:9),3,3) slamPoseMat(1,10:12)';0 0 0 1];
    for i = 1 : size(slamPoseMat,1) - 0
        %         Tp =  [reshape(slamPoseMat(i,1:9),3,3) slamPoseMat(i,10:12)';0 0 0 1];
        Tc =  [reshape(slamPoseMat(i+0,1:9),3,3) slamPoseMat(i+0,10:12)';0 0 0 1];
%         orig2c = inv(Tc) * Tbase;  % Tc \ Tp;  % inv(Tc) * Tp;
        c2orig = inv(Tbase)*Tc; %inv(orig2c);
        %         slamPoseMatLast =  [reshape(slamPoseMat22(end,1:9),3,3) slamPoseMat22(end,10:12)';0 0 0 1];
        %         slamPoseMatNew = slamPoseMatLast * inv(p2c);  % slamPoseMatLast / p2c; % slamPoseMatLast * inv(p2c);
        %         %     slamPoseMatBody = c2b*slamPoseMatNew;
        %         slamPoseMat22 = [slamPoseMat22; [reshape(slamPoseMatNew(1:3,1:3),1,9) slamPoseMatNew(1:3,4)']];
        slamPoseMat22 = [slamPoseMat22; [reshape(c2orig(1:3,1:3),1,9) c2orig(1:3,4)']];
    end
end
figure,plotPath(slamPoseMat22)

BB = fitplane(slamPoseMat22(:,10:12)');BB = BB./norm(BB(1:3));
[~, ~, BB] = robustfitXYZ(slamPoseMat22(:,10:12),[]);
planeErr1 = dot(pextend(slamPoseMat22(:,10:12)'), repmat(BB, 1, size(slamPoseMat22, 1)));

% CC =  fitplane(slamPoseMat(:,10:12)');CC = CC./norm(CC(1:3));
[~, ~, CC] = robustfitXYZ(slamPoseMat(:,10:12),[]);
pose_ = inv(Tbase)*pextend(slamPoseMat(:,10:12)');
pose_diff = pose_(1:3,:)' - slamPoseMat22(:,10:12);
tp = inv(Tbase(1:3,1:3))*CC(1:3);
normDiff = sign(tp(2)).*tp - BB(1:3);
[~, ~, groundPlane2] = robustfitXYZ(slamPoseMat(:,10:12),[]);
planeErr2 = dot(pextend(slamPoseMat(:,10:12)'), repmat(groundPlane2, 1, size(slamPoseMat, 1)));

gVec = BB(1:3);
gVec = sign(gVec(2)).*gVec;

rMatNew = roty(180)*CorrectRot([0;1;0], gVec);
poseNew = rMatNew*slamPoseMat22(:,10:12)';




c2b = [rMatNew [0 0 0]';0 0 0 1];

b2c = inv(c2b);


if 0
    img = imread('D:\Work\marker_mapper1.0.15\20201228_11633\video\L\img_L_0001.png');
    imggg = ImgTransform(img, intrMat, rMatNew);figure,imshow(imggg);
    testBirdView2(intrMat,imggg)
end

end