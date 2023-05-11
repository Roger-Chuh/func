% VI_Verification
R = rodrigues([0.21 0.2 0.3]);
c = [0.1 0.2 0.3]';
SkewSymMat((R-eye(3))*c);
(R-eye(3))*SkewSymMat(c);
a = SkewSymMat(R*c) - SkewSymMat(c);
b = SkewSymMat((R-eye(3))*c);
a-b

R = rodrigues([0.21 0.2 0.3]);
c = [0.1 0.2 0.3]';
SkewSymMat(R*c);
SkewSymMat(R*c) - R*SkewSymMat(c)*R'














D = sym('D',[3,3])
d = sym('d',[3,1])
lambda = sym('lambda');
g = sym('g');

assume(D,'real')
assume(d,'real')
assume(lambda,'real')
assume(g,'real')


expression = det((D - lambda*eye(3,3))^2 - (1/g^2)*(d*d'));
collected = collect(expression, lambda)