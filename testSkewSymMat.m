function testSkewSymMat()



a1 = [0.1 0.2 0.3]';
a2 = [0.1 0.2 -0.22]';
a3 = [0.1 -0.12 -0.22]';
X = SkewSymMat(a1)+SkewSymMat(a2) + SkewSymMat(a3);
Y = SkewSymMat(a1 + a2 + a3);

x = SkewSymMat(2.5.*a1);
y = 2.5.*SkewSymMat(a1);
err_ = x-y;
err = X-Y;

n = [1 2 3]';
% n = n./norm(n);
R = rodrigues([0.3 0.2 -0.1]);
res1 = LeftMultiSkew( n, R);
res1 - SkewSymMat(n)*R
res1 - SkewSymMat(R*n)

res2 =  RightMultiSkew1( R,  n);
res2' - R*SkewSymMat(n)
res2 - SkewSymMat(R*n)


res2 =  RightMultiSkew( R,  n);
res2 - R*SkewSymMat(n)

end
function res = LeftMultiSkew( n, R)

for i =  1 : 3
    res(:,i) = cross(n, R(:,i));
end

end
function res =  RightMultiSkew( R,  n)
%  R = R';
for i = 1 : 3
    res(i,:) = cross(R(i,:),n);
end
end
function res =  RightMultiSkew1( R,  n)
 R = R';
for i = 1 : 3
    res(:,i) = cross(R(:,i),n);
end
end