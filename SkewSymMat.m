function s = SkewSymMat(v)
% Form skew symmetrical matrix of 3D vector v.
% By Ji Zhou

s = [0, -v(3), v(2); v(3), 0, -v(1); -v(2), v(1), 0];


end

