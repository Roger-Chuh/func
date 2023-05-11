function [groundXYZ, indd, groundPlane] = robustfitXYZ(pt3dGround,ratio)

if ~exist('ratio')
    ratio = 0.5;
end


b = robustfit(pt3dGround(:,1:2),pt3dGround(:,3)); groundPlane = [b(2);b(3);-1;b(1)]; groundPlane = groundPlane./norm(groundPlane(1:3));
if groundPlane(4)>0
    groundPlane = -groundPlane;
end

% % pt3dErr = [groundPlane'*[pt3dGround ones(size(pt3dGround,1),1)]']';
pt3dErr = [pt3dGround ones(size(pt3dGround,1),1)]*groundPlane;
thr = 1;
% % % pt3dErr

% % planeErr = sqrt(sum(pt3dErr.^2)/length(pt3dErr));

while sum(abs(pt3dErr)<thr)<0.85*size(pt3dGround,1)
thr = thr+1;
end
indd = abs(pt3dErr)<thr;
groundXYZ =pt3dGround(indd,:);

end