function [ptIcs, tgtPt3d] = TransformAndProject(pt3d, intrMat, rotMat2Ccs, transVec2Ccs)

% TransformAndProject transforms 3d points in a coordinate system to an
% camera coordinate system and project them to image plane of camera.
% pt3d: numPt by 3 matrix, list of 3d points. Each row is the coordinates
% of a point in the form of [x,y,z]
% intrMat: 3 by 3 intrinsic matrix
% rotMatToCcs: 3 by 3 rotation matrix transforming from current coordinate
% system to camera coordinate system.
% transVecToCcs: 3 by 1 translation vector (in ICS).
% Returned ptIcs is a numPt by 2 matrix with each row a projected image
% point corresponding to each 3d point in pt3d.
%
% By Ji Zhou

tgtPt3d = EucildTransform(pt3d, rotMat2Ccs, transVec2Ccs);
ptIcs = ProjectToImage(tgtPt3d, intrMat);

end