function plotMarkerPose(MarkerPoseData, b2c)
a = MarkerPoseData;
slamPoseMat = [];
for k = 1 : size(a,1)
    r = rodrigues(a(k,3:5));
    t = 1000.*a(k,6:8)';
    T = inv([r t; 0 0 0 1]);
    slamPoseMat(k,:) = [reshape(T(1:3,1:3),1,9) T(1:3,4)'];
end
% figure,plot3(poseMat(:,10), poseMat(:,11), poseMat(:,12),'.r');hold on;plot3(poseMat(1,10), poseMat(1,11), poseMat(1,12),'ob');axis equal
% % figure,plotPath(slamPoseMat)


c2b = inv(b2c);

if 0
    slamPoseMat22 = [reshape(eye(3),1,9) [0 0 0]];
    
    for i = 1 : size(slamPoseMat,1) - 1
        Tp =  [reshape(slamPoseMat(i,1:9),3,3) slamPoseMat(i,10:12)';0 0 0 1];
        Tc =  [reshape(slamPoseMat(i+1,1:9),3,3) slamPoseMat(i+1,10:12)';0 0 0 1];
        p2c = inv(Tc) * Tp;  % Tc \ Tp;  % inv(Tc) * Tp;
        
        slamPoseMatLast =  [reshape(slamPoseMat22(end,1:9),3,3) slamPoseMat22(end,10:12)';0 0 0 1];
        slamPoseMatNew = slamPoseMatLast * inv(p2c);  % slamPoseMatLast / p2c; % slamPoseMatLast * inv(p2c);
        %     slamPoseMatBody = c2b*slamPoseMatNew;
        slamPoseMat22 = [slamPoseMat22; [reshape(slamPoseMatNew(1:3,1:3),1,9) slamPoseMatNew(1:3,4)']];
    end
else
    slamPoseMat22 = []; % [reshape(eye(3),1,9) [0 0 0]];
    Tbase =  [reshape(slamPoseMat(1,1:9),3,3) slamPoseMat(1,10:12)';0 0 0 1];
    for i = 1 : size(slamPoseMat,1) - 0
%         Tp =  [reshape(slamPoseMat(i,1:9),3,3) slamPoseMat(i,10:12)';0 0 0 1];
        Tc =  [reshape(slamPoseMat(i+0,1:9),3,3) slamPoseMat(i+0,10:12)';0 0 0 1];
        orig2c = inv(Tc) * Tbase;  % Tc \ Tp;  % inv(Tc) * Tp;
        c2orig = inv(orig2c);
%         slamPoseMatLast =  [reshape(slamPoseMat22(end,1:9),3,3) slamPoseMat22(end,10:12)';0 0 0 1];
%         slamPoseMatNew = slamPoseMatLast * inv(p2c);  % slamPoseMatLast / p2c; % slamPoseMatLast * inv(p2c);
%         %     slamPoseMatBody = c2b*slamPoseMatNew;
%         slamPoseMat22 = [slamPoseMat22; [reshape(slamPoseMatNew(1:3,1:3),1,9) slamPoseMatNew(1:3,4)']];
        slamPoseMat22 = [slamPoseMat22; [reshape(c2orig(1:3,1:3),1,9) c2orig(1:3,4)']];
    end
    
end

if 0
    if 0
        BB = fitplane(slamPoseMat22(:,10:12)');BB = BB./norm(BB(1:3));
    else
        [~, ~, BB] = robustfitXYZ(slamPoseMat22(:,10:12),[]);
        planeErr1 = dot(pextend(slamPoseMat22(:,10:12)'), repmat(BB, 1, size(slamPoseMat22, 1)));
        figure,hist(planeErr1, 1000)
    end
    gVec = BB(1:3);
    gVec = sign(gVec(2)).*gVec;
    
    c2b_ = roty(180)*CorrectRot([0;1;0], gVec);
    
    c2b = [c2b_ [0 0 0]';0 0 0 1];
    
    c2b_plane = [c2b; BB'];
    
%     planeErr1 = dot(pextend(slamPoseMat22(:,10:12)'), repmat(BB, 1, size(slamPoseMat22,1)));

    
end

bodyPose = [];
for i = 1 : size(slamPoseMat22, 1)
    camPose =  [reshape(slamPoseMat22(i,1:9),3,3) slamPoseMat22(i,10:12)';0 0 0 1];
    rtBody = c2b*camPose;
    bodyPose = [bodyPose; [reshape(rtBody(1:3,1:3),1,9) rtBody(1:3,4)']];
end


figure,subplot(2,2,1);plotPath(slamPoseMat);title('cam in marker');
subplot(2,2,2);plotPath(slamPoseMat22);title('cam in cam');
subplot(2,2,3);plotPath(bodyPose);title('cam in body');
subplot(2,2,4);plot(bodyPose(:,11));grid on;title('cam height');
end