function pixDist = remapRectFisheyeRadTanThinPrism(pixRect, KRect, param, RL)

alpha = 0;

rays = inv(KRect)*pextend(pixRect);


% Rotation: (or affine transformation):

rays2 = RL'*rays;

x = [rays2(1,:)./rays2(3,:);rays2(2,:)./rays2(3,:)];


% Add distortion:
if 0
    xd = apply_distortion(x,distCoeff);
    
    
    % Reconvert in pixels:
    
    px2_ = KDistort(1,1)*(xd(1,:) + alpha*xd(2,:)) + KDistort(1,3);
    py2_ = KDistort(2,2)*xd(2,:) + KDistort(2,3);
    pixDist = [px2_;py2_]';
else
%     param = [KDistort(1,1); KDistort(2,2); KDistort(1,3); KDistort(2,3); distCoeff];
    [pixDist] = projectFisheyeRadTanThinPrism(pextend(pextend(x))', param, eye(3), [0 0 0]');
end
end