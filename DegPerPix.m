function DegPerPix()

close all;

width = 640;
height = 480;
f = 235; % 100
radius = 1;

K = [f 0 (width+1)/2-1; 0 f (height+1)/2-1;0 0 1];

k1 = 0.21841;
k2 = -0.174732;
k3 = 0.0657894;
k4 = -0.0105908;


[xMat, yMat] = meshgrid(0 : width-1, 0 : height-1);
pix = [xMat(:) yMat(:)];
center = [(width+1)/2-1 (height+1)/2-1];

pix_flow = pix - center;
[dir, ~] = NormalizeVector(pix_flow);

if 1
    pix1 = pix + radius.*dir;
    pix2 = pix - radius.*dir;
else
    pix1 = pix - [1 0];
    pix2 = pix + [1 0];
end
[dir1, dist1] = NormalizeVector(pix1 - pix);
[dir2, dist2] = NormalizeVector(pix - pix2);

figure,plot(dir1 - dir);
figure,plot(dir2 - dir);
figure,plot(dist1 - radius);
figure,plot(dist2 - radius);
if 0
    backProj1 = (inv(K)*(pextend(pix1')))';
    backProj2 = (inv(K)*(pextend(pix2')))';
else
    backProj1 = unprojectKB8( K(1,1),  K(2,2),  K(1,3),  K(2,3),  k1,  k2,  k3,  k4, pix1);
    backProj2 = unprojectKB8( K(1,1),  K(2,2),  K(1,3),  K(2,3),  k1,  k2,  k3,  k4, pix2);
    [~, dis1] = NormalizeVector(backProj1);
    [~, dis2] = NormalizeVector(backProj2);
end
backProj11 =  normr(backProj1);
backProj22 =  normr(backProj2);
validFlag = abs(1-dis1) < 0.00001 &  abs(1-dis2) < 0.00001;
if 0
    cov0 = reshape((acos(dot(backProj11', backProj22'))), height, width);
    cov = reshape(rad2deg(acos(dot(backProj11', backProj22'))), height, width);
    %     figure,surface(cov);
    figure,surface(cov);%axis equal;
else
    [~, err] = NormalizeVector(backProj11 - backProj22);
    cov = reshape(err, height, width);
    noise = 1./(f*cov);
%     noise(noise>0.1) = 0;
    noise(noise>0.9) = 0;
    noise(~validFlag) = 0;
    figure,surface(noise);
    figure,contour(noise, 100);colorbar; axis equal;
    
end

end
function [theta,d_func_d_theta1] = solveTheta( k1,  k2,  k3,  k4,   r_theta ,d_func_d_theta, ITER)


theta = r_theta;
for  i = ITER:-1:1
    theta2 = theta .* theta;
    
    func = k4 .* theta2;
    func = func + k3;
    func = func .* theta2;
    func = func + k2;
    func = func .* theta2;
    func = func + k1;
    func = func .* theta2;
    func = func + 1.0;
    func = func .* theta;
    
    d_func_d_theta = (9) .* k4 .* theta2;
    d_func_d_theta = d_func_d_theta + (7) * k3;
    d_func_d_theta = d_func_d_theta .* theta2;
    d_func_d_theta = d_func_d_theta + (5) * k2;
    d_func_d_theta = d_func_d_theta .* theta2;
    d_func_d_theta = d_func_d_theta + (3) * k1;
    d_func_d_theta = d_func_d_theta .* theta2;
    d_func_d_theta = d_func_d_theta + (1);
    
    
    %theta = theta + (r_theta - func) / d_func_d_theta;
    theta_fix = (r_theta - func) ./ d_func_d_theta;
    theta = theta+ theta_fix;
%     if (abs(theta_fix) < 1e-8)
%         break;
%     end
end
d_func_d_theta1 = d_func_d_theta;
end
function p3d = unprojectKB8( fx,  fy,  cx,  cy,  k1,  k2,  k3,  k4, proj)


mx = (proj(:,1) - cx) ./ fx;
my = (proj(:,2) - cy) ./ fy;




theta = 0.*ones(size(proj,1),1);
sin_theta = 0.*ones(size(proj,1),1);
cos_theta = 1.*ones(size(proj,1),1);
thetad = sqrt(mx .* mx + my .* my);

thetad = min([[[max([thetad,-3.141592653/2*ones(length(thetad),1)]')]' 3.141592653/2*ones(length(thetad),1)]]')';
scaling = 1.0.*ones(size(proj,1),1);
d_func_d_theta = 0.*ones(size(proj,1),1);
flag = thetad > 1e-8;
if (1)
    [theta(flag),d_func_d_theta1(flag)] = solveTheta(k1, k2, k3, k4, thetad(flag), d_func_d_theta(flag),20);
    d_func_d_theta(flag) = d_func_d_theta1(flag);
    sin_theta(flag) = sin(theta(flag));
    cos_theta(flag) = cos(theta(flag));
    scaling(flag) = sin_theta(flag) ./ thetad(flag);
end

p3d = [mx .* scaling my .* scaling cos_theta];


end