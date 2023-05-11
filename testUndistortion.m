function testUndistortion()

close all
J = checkerboard(40,4,7);
J = J > 0.5;
J = uint8(255.*J);
if 0
    width_virtual = size(J,2);
    height_virtual = size(J,1);
    width_gt = width_virtual+100;
    height_gt = height_virtual+100;
else
    if 0
        width_gt = size(J,2);
        height_gt = size(J,1);
        width_virtual = width_gt+100;
        height_virtual = height_gt+100;
    else
        scale = [1.5 2];
        width_gt = size(J,2);
        height_gt = size(J,1);
        width_virtual = width_gt*scale(2);
        height_virtual = height_gt*scale(1);
    end
    
end
% intrMat = [499 0 319; 0 501 241;0 0 1];
% intrMat_virtual = [510 0 (width_virtual+1)/2; 0 520 (height_virtual+1)/2;0 0 1];
intrMat_gt = [400 0 (width_gt+1)/2+10; 0 450 (height_gt+1)/2-10;0 0 1];
% intrMat_gt = [400 0 (width_gt+1)/2+0; 0 450 (height_gt+1)/2-0;0 0 1];
intrMat_virtual = [510 0 (width_virtual+1)/2; 0 520 (height_virtual+1)/2;0 0 1];
intrMat_virtual = [intrMat_gt(1,1)*1.1*scale(1) 0 (width_virtual+1)/2+30; 0 intrMat_gt(2,2)*0.7*scale(2) (height_virtual+1)/2-50;0 0 1];
% k = -[0.3 0.2 0.00001 -0.00003 0.003]';
k = -[0.02 0.008 0.00001 -0.00003 0]';
% k = [0.178 0.1 0.00001 -0.00003 0]';
% k = [1.5 0.9 0.00001 -0.00003 0]';
k = [1.0 0.7 0.00001 -0.00003 0]';
k = [0.7 0.5 0.00001 -0.00003 0]';
R = eye(3);
R = rodrigues([0.03,-0.04 -0.05]);
% R = eye(3);


[xMat_virtual, yMat_virtual] = meshgrid(1:width_virtual, 1:height_virtual);
pix_virtual = [xMat_virtual(:) yMat_virtual(:)];
[xMat_gt, yMat_gt] = meshgrid((intrMat_virtual(1,3)-intrMat_gt(1,3)+1):(intrMat_virtual(1,3)-intrMat_gt(1,3)+width_gt), (intrMat_virtual(2,3)-intrMat_gt(2,3)+1):(intrMat_virtual(2,3)-intrMat_gt(2,3)+height_gt));
[xMat_gt, yMat_gt] = meshgrid(1 : width_gt, 1 : height_gt);
pix_gt = [xMat_gt(:) yMat_gt(:)];

%target图（virtual，distorted）的坐标是grid，需要把它映射到gt图像


if 0
    pixUndist_gt = remapRect(pix_virtual, intrMat_virtual, intrMat_gt, R, k);
    pixUndist_gt_no_rot = remapRect(pix_virtual, intrMat_virtual, intrMat_gt, eye(3), k);
    pixOrig_virtual = Orig2Rect(pixUndist_gt, intrMat_gt, intrMat_virtual, R, k);
    % pixOrig = Orig2Rect(pixUndist, intrMat_virtual, intrMat_gt, R, k);
    % pixOrig2 = Orig2Rect(pix, intrMat_virtual, intrMat_gt, R, k);
    % pixOrig_gt2 = Orig2Rect(pix_virtual, intrMat_gt, intrMat_virtual, R, k);
    if 0
        pixOrig_gt2 = Orig2Rect(pix_gt, intrMat_gt, intrMat_virtual, R, k);
    else
        pixOrig_gt2_ = Orig2Rect(pix_gt, intrMat_gt, intrMat_virtual, R, k);
        pixOrig_gt2 = Orig2Rect(pix_gt, intrMat_gt, intrMat_gt, R, k);
    end
else
    pixUndist_gt = Orig2Rect(pix_virtual, intrMat_virtual, intrMat_gt, R, k);
    pixUndist_gt_no_rot = Orig2Rect(pix_virtual, intrMat_virtual, intrMat_gt, eye(3), k);
    pixOrig_virtual = remapRect(pixUndist_gt, intrMat_gt, intrMat_virtual, R, k);
    % pixOrig = Orig2Rect(pixUndist, intrMat_virtual, intrMat_gt, R, k);
    % pixOrig2 = Orig2Rect(pix, intrMat_virtual, intrMat_gt, R, k);
    % pixOrig_gt2 = Orig2Rect(pix_virtual, intrMat_gt, intrMat_virtual, R, k);
    if 0
        pixOrig_gt2 = remapRect(pix_gt, intrMat_gt, intrMat_virtual, R, k);
    else
        pixOrig_gt2_ = remapRect(pix_gt, intrMat_gt, intrMat_virtual, R, k);
        if 0
            intrMat_gt_mod = intrMat_virtual./scale;
            intrMat_gt_mod(3,3) = 1;
        else
            intrMat_gt_mod = eye(3);
            intrMat_gt_mod(1,1) = intrMat_virtual(1,1)/scale(2);
            intrMat_gt_mod(1,3) = intrMat_virtual(1,3)/scale(2);
            intrMat_gt_mod(2,2) = intrMat_virtual(2,2)/scale(1);
            intrMat_gt_mod(2,3) = intrMat_virtual(2,3)/scale(1);
        end
        if 0
            pixOrig_gt2 = remapRect(pix_gt, intrMat_gt, intrMat_gt, R, k);
        else
            pixOrig_gt2 = remapRect(pix_gt, intrMat_gt, intrMat_gt_mod, R, k);
        end
    end
    
    
end

% 对【无畸变】施加【枕形畸变】变成【枕形畸变】（投影仪的成像过程）
J_pincushion = interp2(xMat_gt,yMat_gt,double(J),reshape(pixUndist_gt(:,1),height_virtual,width_virtual),reshape(pixUndist_gt(:,2),height_virtual,width_virtual));

%对【枕形畸变】施加【桶形畸变】变成【无畸变】
J0 = interp2(xMat_virtual,yMat_virtual,double(J_pincushion),reshape(pixOrig_gt2_(:,1),height_gt,width_gt),reshape(pixOrig_gt2_(:,2),height_gt,width_gt));

%% 对【无畸变】施加【桶形畸变】变成【桶形畸变】
J_barrel = interp2(xMat_gt,yMat_gt,double(J),reshape(pixOrig_gt2(:,1),height_gt,width_gt),reshape(pixOrig_gt2(:,2),height_gt,width_gt));

%对【桶形畸变】施加【枕形畸变】变成【无畸变】
J00 = interp2(xMat_gt,yMat_gt,double(J_barrel),reshape(pixUndist_gt(:,1),height_virtual,width_virtual),reshape(pixUndist_gt(:,2),height_virtual,width_virtual));
J000 = interp2(xMat_gt,yMat_gt,double(J_barrel),reshape(pixUndist_gt_no_rot(:,1),height_virtual,width_virtual),reshape(pixUndist_gt_no_rot(:,2),height_virtual,width_virtual));
if 0
    figure,subplot(2,2,1);imshow(J,[]);title('原始图像','Color','r','FontSize',20);
    subplot(2,2,2);imshow(J_pincushion,[]);title(sprintf('光机的输出\n包含枕形畸变与旋转（畸变与旋转矩阵由标定获得）'),'Color','r','FontSize',20);
    subplot(2,2,3);imshow(J_barrel,[]);title(sprintf('对原始图像施加桶形畸变以及反方向的旋转\n（由枕形畸变参数和旋转矩阵求反函数得到）'),'Color','b','FontSize',20);
    subplot(2,2,4);imshow(J00,[]);title(sprintf('桶形畸变的图像经过光机作用后\n恢复成不含畸变且旋转补偿后的原始图像'),'Color','b','FontSize',20);
else
    figure,subplot(2,2,1);imshow(J,[]);title('\bf 原始图像','Color','r','FontSize',20);
    subplot(2,2,2);imshow(J_pincushion,[]);title({'\bf 光机的输出';'包含枕形畸变与旋转（畸变与旋转矩阵由标定获得）'},'Color','r','FontSize',20);
    subplot(2,2,3);imshow(J_barrel,[]);title({'\bf 对原始图像施加桶形畸变以及反方向的旋转';'（由枕形畸变参数和旋转矩阵求反函数得到）'},'Color','b','FontSize',20);
    subplot(2,2,4);imshow(J00,[]);title({'\bf 桶形畸变的图像经过光机作用后';'恢复成不含畸变且旋转补偿后的原始图像'},'Color','b','FontSize',20);
end

%        figure,imshowpair(J,J00);
if 0
    figure,imshowpair(imresize(J,scale),J00);
else
    figure,imshowpair(imresize(J,[height_virtual, width_virtual]),J00);
end
       if 0
           num = 200000;
           ind = randperm(size(pixUndist,1));
           
           pixUndist_ = pixUndist(1:ind(1:num),:);
           pixOrig_ = pixOrig(1:ind(1:num),:);
       end

figure,imshow(zeros(height_virtual,width_virtual));hold on;plot(pixUndist_gt(:,1),pixUndist_gt(:,2),'.r');
% % figure,subplot(1,2,1),imshow(zeros(480,640));hold on;plot(pixUndist_(:,1),pixUndist_(:,2),'.r');
% %        subplot(1,2,2),imshow(zeros(480,640));hold on;plot(pixOrig_(:,1),pixOrig_(:,2),'.g');

err = pixOrig_virtual - pix_virtual;

if 1
    figure,plot(err);
end


end
function pixDist = remapRect(pixRect, KRect, KDistort,  R, distCoeff)

alpha = 0;
pixRectHomo = [pixRect'; ones(1,size(pixRect,1))];
rays = inv(KRect)*pixRectHomo;


% Rotation: (or affine transformation):

rays2 = R'*rays;

x = [rays2(1,:)./rays2(3,:);rays2(2,:)./rays2(3,:)];


% Add distortion:
xd = apply_distortion(x,distCoeff);


% Reconvert in pixels:

px2_ = KDistort(1,1)*(xd(1,:) + alpha*xd(2,:)) + KDistort(1,3);
py2_ = KDistort(2,2)*xd(2,:) + KDistort(2,3);
pixDist = [px2_;py2_]';

end
function [xd,dxddk] = apply_distortion(x,k)


% Complete the distortion vector if you are using the simple distortion model:
length_k = length(k);
if length_k <5 ,
    k = [k ; zeros(5-length_k,1)];
end;


[m,n] = size(x);

% Add distortion:

r2 = x(1,:).^2 + x(2,:).^2;

r4 = r2.^2;

r6 = r2.^3;


% Radial distortion:

cdist = 1 + k(1) * r2 + k(2) * r4 + k(5) * r6;

if nargout > 1,
    dcdistdk = [ r2' r4' zeros(n,2) r6'];
end;


xd1 = x .* (ones(2,1)*cdist);



if nargout > 1,
    dxd1dk = zeros(2*n,5);
    dxd1dk(1:2:end,:) = (x(1,:)'*ones(1,5)) .* dcdistdk;
    dxd1dk(2:2:end,:) = (x(2,:)'*ones(1,5)) .* dcdistdk;
end;


% tangential distortion:

a1 = 2.*x(1,:).*x(2,:);
a2 = r2 + 2*x(1,:).^2;
a3 = r2 + 2*x(2,:).^2;

delta_x = [k(3)*a1 + k(4)*a2 ;
    k(3) * a3 + k(4)*a1];



if nargout > 1,
    ddelta_xdk = zeros(2*n,5);
    ddelta_xdk(1:2:end,3) = a1';
    ddelta_xdk(1:2:end,4) = a2';
    ddelta_xdk(2:2:end,3) = a3';
    ddelta_xdk(2:2:end,4) = a1';
end;

xd = xd1 + delta_x;

if nargout > 1,
    dxddk = dxd1dk + ddelta_xdk ;
    if length_k < 5,
        dxddk = dxddk(:,1:length_k);
    end;
end;


end

function pixRect = Orig2Rect(pix, intrMatOld, intrMatNew, R, kc)

[pixUndist] = normalize_pixel(pix',[intrMatOld(1,1);intrMatOld(2,2)],[intrMatOld(1,3);intrMatOld(2,3)],kc,0);
pixUndistHomo = [pixUndist; ones(1, size(pixUndist,2))];
pixUndistR = R*pixUndistHomo;
pixRect = intrMatNew*pixUndistR;
pixRect = [pixRect(1,:)./pixRect(3,:); pixRect(2,:)./pixRect(3,:)];
pixRect = pixRect(1:2,:)';

end
function [xn] = normalize_pixel(x_kk,fc,cc,kc,alpha_c)

%normalize
%
%[xn] = normalize_pixel(x_kk,fc,cc,kc,alpha_c)
%
%Computes the normalized coordinates xn given the pixel coordinates x_kk
%and the intrinsic camera parameters fc, cc and kc.
%
%INPUT: x_kk: Feature locations on the images
%       fc: Camera focal length
%       cc: Principal point coordinates
%       kc: Distortion coefficients
%       alpha_c: Skew coefficient
%
%OUTPUT: xn: Normalized feature locations on the image plane (a 2XN matrix)
%
%Important functions called within that program:
%
%comp_distortion_oulu: undistort pixel coordinates.

if nargin < 5,
    alpha_c = 0;
    if nargin < 4;
        kc = [0;0;0;0;0];
        if nargin < 3;
            cc = [0;0];
            if nargin < 2,
                fc = [1;1];
            end;
        end;
    end;
end;


% First: Subtract principal point, and divide by the focal length:
x_distort = [(x_kk(1,:) - cc(1))/fc(1);(x_kk(2,:) - cc(2))/fc(2)];

% Second: undo skew
x_distort(1,:) = x_distort(1,:) - alpha_c * x_distort(2,:);

if norm(kc) ~= 0,
    % Third: Compensate for lens distortion:
    xn = comp_distortion_oulu(x_distort,kc);
else
    xn = x_distort;
end;

end

function [x] = comp_distortion_oulu(xd,k);

%comp_distortion_oulu.m
%
%[x] = comp_distortion_oulu(xd,k)
%
%Compensates for radial and tangential distortion. Model From Oulu university.
%For more informatino about the distortion model, check the forward projection mapping function:
%project_points.m
%
%INPUT: xd: distorted (normalized) point coordinates in the image plane (2xN matrix)
%       k: Distortion coefficients (radial and tangential) (4x1 vector)
%
%OUTPUT: x: undistorted (normalized) point coordinates in the image plane (2xN matrix)
%
%Method: Iterative method for compensation.
%
%NOTE: This compensation has to be done after the subtraction
%      of the principal point, and division by the focal length.


if length(k) == 1,
    
    [x] = comp_distortion(xd,k);
    
else
    
    k1 = k(1);
    k2 = k(2);
    k3 = k(5);
    p1 = k(3);
    p2 = k(4);
    
    x = xd; 				% initial guess
    
    for kk=1:20,
        
        r_2 = sum(x.^2);
        k_radial =  1 + k1 * r_2 + k2 * r_2.^2 + k3 * r_2.^3;
        delta_x = [2*p1*x(1,:).*x(2,:) + p2*(r_2 + 2*x(1,:).^2);
            p1 * (r_2 + 2*x(2,:).^2)+2*p2*x(1,:).*x(2,:)];
        x = (xd - delta_x)./(ones(2,1)*k_radial);
        
    end;
    
end;

end