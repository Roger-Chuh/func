function [img_undist, remapX, remapY] = GenUndistortMap(img, K_pinhole, width, height, param)
[xMat, yMat] = meshgrid(1 : width, 1 : height);
pix = [xMat(:) yMat(:)];
rays = (inv(K_pinhole) * pextend(pix'))';
[pt2d] = ProjectKB8(rays(:,1), rays(:,2), rays(:,3), param);

remapX = reshape(pt2d(:,1), height, width);
remapY = reshape(pt2d(:,2), height, width);

img_undist = uint8(interp2(xMat, yMat, double(img(:,:,1)), remapX, remapY));
end
function y = pextend(x)
y = [x; ones(1,size(x,2))];
end
function [pt2d] = ProjectKB8(X, Y, Z, param)
fx = param(1);
fy = param(2);
cx = param(3);
cy = param(4);
k1 = param(5);
k2 = param(6);
k3 = param(7);
k4 = param(8);

x = X./Z;
y = Y./Z;
z = 1;%Z;


r2 = x .* x + y .* y;
r = sqrt(r2);

R2 = X .* X + Y .* Y;
R = sqrt(R2);

theta = atan(r);
theta2 = theta .* theta;
theta4 = theta2 .* theta2;
theta6 = theta4 .* theta2;
theta8 = theta6 .* theta2;

r_theta = theta .* (1 + k1 .* theta2 + k2 .* theta4 + k3 .* theta6 + k4 .* theta8);

if r > 1e-18
    norm_inv = 1./r;
else
    norm_inv = ones(length(r), 1);
end
% norm_inv = r > 1e-8 ? double(1.0) / r : 1;

mx = r_theta .* x .* norm_inv;
my = r_theta .* y .* norm_inv;


pt2d(:,1) = fx .* mx + cx;
pt2d(:,2) = fy .* my + cy;


end