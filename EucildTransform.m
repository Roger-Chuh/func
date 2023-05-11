function transPt = EucildTransform(pt, rotMat, transVec)

% EucildTransform transforms points from source coordinate system to target
% coordinate system in Euclidian space.
% pt: numPt by 3 or 2 matrix. Each row is the coordinates of a point in the
% form of [x,y,z] or [x,y] in source coordinate system
% rotMat: 3 by 3 or 2 by 2 rotation matrix that transforms point from
% source to target coordinate system
% transVec: 3 by 1 or 2 by 1 translation vector (in target coordinate
% system) that transforms point from source to target coordinate system.
% Returned transPt is transformed points of pt, a matrix with the same size
% as pt.
%
% By Ji Zhou

transPt = bsxfun(@plus, rotMat*pt', transVec)';

end