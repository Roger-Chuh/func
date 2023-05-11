function [Y,dYdom,dYdT] = rigid_motion_cube(X,om,T,cubeIdx,config);

%rigid_motion.m
%
%[Y,dYdom,dYdT] = rigid_motion(X,om,T)
%
%Computes the rigid motion transformation Y = R*X+T, where R = rodrigues(om).
%
%INPUT: X: 3D structure in the world coordinate frame (3xN matrix for N points)
%       (om,T): Rigid motion parameters between world coordinate frame and camera reference frame
%               om: rotation vector (3x1 vector); T: translation vector (3x1 vector)
%
%OUTPUT: Y: 3D coordinates of the structure points in the camera reference frame (3xN matrix for N points)
%        dYdom: Derivative of Y with respect to om ((3N)x3 matrix)
%        dYdT: Derivative of Y with respect to T ((3N)x3 matrix)
%
%Definitions:
%Let P be a point in 3D of coordinates X in the world reference frame (stored in the matrix X)
%The coordinate vector of P in the camera reference frame is: Y = R*X + T
%where R is the rotation matrix corresponding to the rotation vector om: R = rodrigues(om);
%
%Important function called within that program:
%
%rodrigues.m: Computes the rotation matrix corresponding to a rotation vector



if nargin < 3+2,
   T = zeros(3,1);
   if nargin < 2+2,
      om = zeros(3,1);
      if nargin < 1+2,
         error('Need at least a 3D structure as input (in rigid_motion.m)');
         return;
      end;
   end;
end;


[R,dRdom] = rodrigues(om);

[m,n] = size(X);
% Y = zeros(3, length(cell2mat(cubeIdx')));
T1 = [R T; 0 0 0 1];
T2 = T1*config.rt2;
T3 = T1*config.rt3;

Y1 =           R*X(:,cubeIdx{1}) + repmat(T,        [1 length(cubeIdx{1})]);
Y2 = T2(1:3,1:3)*X(:,cubeIdx{2}) + repmat(T2(1:3,4),[1 length(cubeIdx{2})]);
Y3 = T3(1:3,1:3)*X(:,cubeIdx{3}) + repmat(T3(1:3,4),[1 length(cubeIdx{3})]);
Y = [Y1 Y2 Y3];

if nargout > 1,
   

dYdR = zeros(3*n,9);
dYdT = zeros(3*n,3);

dYdR(1:3:end,1:3:end) =  X';
dYdR(2:3:end,2:3:end) =  X';
dYdR(3:3:end,3:3:end) =  X';

dYdT(1:3:end,1) =  ones(n,1);
dYdT(2:3:end,2) =  ones(n,1);
dYdT(3:3:end,3) =  ones(n,1);

dYdom = dYdR * dRdom;

end
end




