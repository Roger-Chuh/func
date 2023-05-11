function [rotMatLeft, rotMatRight, intrMatLeftNew, intrMatRightNew, kcLeftNew, kcRightNew] = ...
    GetRectifyPlane(stereoParam, imgSize, varargin)

global cfg



if (nargin <= 2)
    replaceF = [];
elseif (nargin == 3)
    replaceF = varargin{1};
else
    error('Too many input arguments');
end





focLeft = stereoParam.focLeft;
focRight = stereoParam.focRight;
cenLeft = stereoParam.cenLeft;
cenRight = stereoParam.cenRight;
alphaLeft = stereoParam.alphaLeft;
alphaRight = stereoParam.alphaRight;
kcLeft = stereoParam.kcLeft;
kcRight = stereoParam.kcRight;
rotVecRef = stereoParam.rotVecRef;
transVecRef = stereoParam.transVecRef;

% Bring the 2 cameras in the same orientation by rotating them "minimally":
if 1 % isempty(cfg.KK_newL)
    rotMatCFR = rodrigues(-rotVecRef/2); % rotation converting right camera frame to common frame
    rotMatCFL = rotMatCFR';            % rotation converting left camera frame to common frame
else
    rotMatCFR = eye(3); %rodrigues(-rotVecRef/2); % rotation converting right camera frame to common frame
    rotMatCFL = rodrigues(rotVecRef); % rotMatCFR';
end
%% 
eppVec = rotMatCFR * transVecRef; % vector of epipolar line in common rectangular frame
% % % % % eppVec = rotMatCFR * stereoParam.baselineVec;




% Rotate both cameras so as to bring the translation vector in alignment
% with the (1;0;0) axis:
if abs(eppVec(1)) > abs(eppVec(2))
    uu = [1;0;0]; % Horizontal epipolar lines
    stereoType = 0;
else
    uu = [0;1;0]; % Vertical epipolar lines
    stereoType = 1;
end
if (dot(eppVec, uu) < 0)
    uu = -uu;
end
ww = cross(eppVec,uu);
ww = ww/norm(ww); % rotation axis
rotVecEpp = acos(abs(dot(eppVec,uu))/(norm(eppVec)*norm(uu)))*ww; % acos is angle to rotate
rotMatEpp = rodrigues(rotVecEpp);



% Global rotations to be applied to both views:
rotMatLeft = rotMatEpp * rotMatCFL;
rotMatRight = rotMatEpp * rotMatCFR;

baseline_new = rotMatRight * transVecRef;

imgH = imgSize(1);
imgW = imgSize(2);
% Computation of the *new* intrinsic parameters for both left and right
% cameras:
% Vertical focal length *MUST* be the same for both images (here, we are
% trying to find a focal length that retains as much information contained
% in the original distorted images):
if (kcLeft(1) < 0)
    focYLeftNew = focLeft(2) * (1 + kcLeft(1)*(imgW^2 + imgH^2)/(4*focLeft(2)^2));
else
    focYLeftNew = focLeft(2);
end
if (kcRight(1) < 0)
    focYRightNew = focRight(2) * (1 + kcRight(1)*(imgW^2 + imgH^2)/(4*focRight(2)^2));
else
    focYRightNew = focRight(2);
end

%%  try the larger focal length, see what happens.
focYNew = min(focYLeftNew, focYRightNew);
focYNew0 = focYNew;
if ~isempty(replaceF)
    if focYNew < replaceF
        focYNew = replaceF + cfg.ext_focal;  % 1300;800;
    else
        focYNew = focYNew + cfg.ext_focal;
    end
else
    try
        focYNew = focYNew + cfg.ext_focal;
    catch
        focYNew = focYNew + cfg.ext_focal;
    end
end
% % % % % focYNew = max(focYLeftNew, focYRightNew);
% focYNew = mean([focYLeftNew, focYRightNew]);






% For simplicity, let's pick the same value for the horizontal focal length
% as the vertical focal length (resulting into square pixels):

if 1 % abs(focYNew0 - focLeft(1)) < abs(focYNew0 - focLeft(2))
    
    focLeftNew = round([focYNew; cfg.fxy_ratio * focYNew]);
    focRightNew = round([focYNew; cfg.fxy_ratio * focYNew]);
else
    
    
end
% Select the new principal points to maximize the visible area in the
% rectified images
pt = [0 imgW-1 imgW-1 0; 0 0 imgH-1 imgH-1];
pt = [imgW/2-1 imgW-1 imgW/2-1 0; 0 imgH/2-1 imgH-1 imgH/2-1];

intrMatLeft = [focLeft(1),0,cenLeft(1); 0,focLeft(2),cenLeft(2); 0,0,1];
% % % % %  ptUndistLeft = squeeze(cv.undistortPoints(permute(pt,[3,2,1]), intrMatLeft, kcLeft))';
if 1
    ptUndistLeft = normalize_pixel([0  imgW-1 imgW-1 0; 0 0 imgH-1 imgH-1],focLeft,cenLeft,kcLeft,0);
else
    ptUndistLeft = normalize_pixel(pt,focLeft,cenLeft,kcLeft,0);
end

cenLeftNew = [(imgW-1)/2;(imgH-1)/2] - ...
    mean(project_points2(...
                         [ptUndistLeft;[1 1 1 1]],...
                         rodrigues(rotMatLeft),zeros(3,1),focLeftNew,[0;0],zeros(5,1),0 ...
                         ),...
         2);
     
intrMatRight = [focRight(1),0,cenRight(1); 0,focRight(2),cenRight(2); 0,0,1];
% % % % %  ptUndistRight = squeeze(cv.undistortPoints(permute(pt,[3,2,1]), intrMatRight, kcRight))';
if 1
    ptUndistRight = normalize_pixel([0  imgW-1 imgW-1 0; 0 0 imgH-1 imgH-1],focRight,cenRight,kcRight,0);
else
    ptUndistRight = normalize_pixel(pt,focRight,cenRight,kcRight,0);
end
cenRightNew = [(imgW-1)/2;(imgH-1)/2] - ...
    mean(project_points2(...
                         [ptUndistRight;[1 1 1 1]],...
                         rodrigues(rotMatRight),zeros(3,1),focRightNew,[0;0],zeros(5,1),0 ...
                         ),...
         2);

% For simplivity, set the principal points for both cameras to be the
% average of the two principal points.
if ~stereoType
    % Horizontal stereo
    cenYNew = (cenLeftNew(2) + cenRightNew(2))/2;
    cenLeftNew = [cenLeftNew(1); cenYNew];
    cenRightNew = [cenRightNew(1); cenYNew];
    
    if 1 % to enforce cenXL = cenXY
        cenXNew = (cenLeftNew(1) + cenRightNew(1))/2;
        cenLeftNew = [cenXNew; cenYNew];
        cenRightNew = [cenXNew; cenYNew];
    end


else
    % Vertical stereo
    cenXNew = (cenLeftNew(1) + cenRightNew(1))/2;
    cenLeftNew = [cenXNew; cenLeftNew(2)];
    cenRightNew = [cenXNew; cenRightNew(2)];
    
    if 1 % to enforce cenXL = cenXY
        cenYNew = (cenLeftNew(2) + cenRightNew(2))/2;
        cenLeftNew = [cenXNew; cenYNew];
        cenRightNew = [cenXNew; cenYNew];
    end
end

% Assume after rectification, there is no skew and distortion
alphaLeftNew = 0;
alphaRightNew = 0;
kcLeftNew = [0;0;0;0;0];
kcRightNew = [0;0;0;0;0];

intrMatLeftNew = [focLeftNew(1) focLeftNew(1)*alphaLeftNew cenLeftNew(1);0 focLeftNew(2) cenLeftNew(2); 0 0 1];
intrMatRightNew = [focRightNew(1) focRightNew(1)*alphaRightNew cenRightNew(1);0 focRightNew(2) cenRightNew(2); 0 0 1];

end