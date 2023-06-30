function res = LeftMultiSkew(n, R)
res = eye(3);
for i = 1:3
    res(:,i) = cross(n, R(:,i));
end

end

% R  = rodrigues([0.1 0.2 0.3]);
% SkewSymMat([2.1 2.3 -4.6])*R - LeftMultiSkew([2.1 2.3 -4.6]', R)