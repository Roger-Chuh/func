function ptIcs = ProjectToImage(ptCcs, intrMat)

% ProjectToImage projects 3d points in camera coordinate system to image
% coordinate system of that camera.
% ptCcs: numPt by 3 matrix with each row is the coordinates of a point in
% camera coordinate system
% intrMat: 3 by 3 camera intrinsic matrix
% Returned ptIcs is numPt by 2 matrix which is the image point coordinates
% projected by ptCcs.
%
% By Ji Zhou

ptIcs = ptCcs*intrMat';
ptIcs = bsxfun(@rdivide, ptIcs(:,1:2), ptIcs(:,3));

end