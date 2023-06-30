function pixDistorted = remapRect2(pixRect_, Knew,distCoeffnew, Kold,distCoeffold, RL)
% pixRect_ = pixRect_;
%% pixUndist表示不包含distCoeffRect畸变的bearing vector
[pixRect, pixUndist] = Orig2Rect2(pixRect_, Knew, Knew, eye(3), distCoeffnew);






alpha = 0;


rays = inv(Knew)*pextend(pixRect');

rays = pixUndist;
% Rotation: (or affine transformation):

rays2 = RL'*rays;

x = [rays2(1,:)./rays2(3,:);rays2(2,:)./rays2(3,:)];


% Add distortion:
xd = apply_distortion(x,distCoeffold);


% Reconvert in pixels:

px2_ = Kold(1,1)*(xd(1,:) + alpha*xd(2,:)) + Kold(1,3);
py2_ = Kold(2,2)*xd(2,:) + Kold(2,3);
pixDistorted = [px2_;py2_]';

end