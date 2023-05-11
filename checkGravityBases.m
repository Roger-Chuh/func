function checkGravityBases()

g = [0.486262; -9.7774;0.368008];
Ag = [  -9.7774  0.178948;
-0.486262  -3.59816;
0   -95.834];


gg = g./norm(g);
Ag2 = [1-gg(1)^2/(1+gg(3)) -gg(1)*gg(2)/(1+gg(3));-gg(1)*gg(2)/(1+gg(3)) 1-gg(2)^2/(1+gg(3)); -gg(1) -gg(2)];
axis = cross(Ag2(:,1), Ag2(:,2));

R_ = [Ag2 axis];
det(R_)

A0 = [ cross(g,[0 0 -1]') cross(g,cross(g,[0 0 -1]'))];
A = [ cross(g,[0 0 -1]')./norm(cross(g,[0 0 -1]')) cross(g,cross(g,[0 0 -1]'))./norm(cross(g,cross(g,[0 0 -1]')))];
R = [-g'./norm(g); A' ];
RR = [-g'; Ag'];
RRR = [-g'./norm(g); Ag'];
R*g
RR*g
RRR*g


dx = [0.002; -0.003];

g_new = g + A(:,1)*dx(1) + A(:,2)*dx(2);
end