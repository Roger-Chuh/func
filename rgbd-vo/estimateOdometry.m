function T = estimateOdometry(K, imgRGB_a, imgRGB_b, imgZ_a, imgZ_b, T0)

a_k = transformToParameters(T0);
a_k0 = a_k;
%             imgZ_a = imageA.getDepthImage();
if ndims(imgRGB_a) == 3
    imgRGB_a = double(rgb2gray(imgRGB_a));
    imgRGB_b = double(rgb2gray(imgRGB_b));
    imgRGB_a = imgRGB_a/255.0;
    imgRGB_b = imgRGB_b/255.0;
end

[rows, cols] = size(imgZ_a);
numPts = rows*cols;

%imgRGB_a = reshape(imgRGB_a,numPts,3);
imgI_a = reshape(imgRGB_a,1,numPts);

%             imgZ_b = double(imageB.getDepthImage());
[gradZ_x, gradZ_y] = gradient(imgZ_b);
% imgRGB_b = imgRGB_b/255.0;
[gradI_x, gradI_y] = gradient(imgRGB_b);
imgZ_b = reshape(imgZ_b,1,numPts);
%imgRGB_b = reshape(imgRGB_b,numPts,3);
imgI_b = reshape(imgRGB_b,1,numPts);


pts = getPointCloud(imgZ_a, cols, rows,K);
T = parametersToTransform(a_k);
prev_sqError = Inf;
MAXITER=200;
for iter=1:MAXITER
    %warpA = imageA.warpImage(T);
    %warpA_Z = warpA.getDepthImage();
    %imgB_Z = imageB.getDepthImage();
    %figure(1), subplot(1,3,1), imshow(imgB_Z,[0 8]);
    %figure(1), subplot(1,3,2), imshow(warpA_Z,[0 8]);
    %figure(1), subplot(1,3,3), imshow(abs(warpA_Z - imgB_Z),[ ]);
    %error = abs(warpA_Z - imgB_Z);
    %warpA_Z = reshape(warpA_Z,1,numPts);
    %error = reshape(error,1,numPts);
    %error_mask=zeros(1,numPts);
    %error_mask(mask==1)=1;
    %error_mask(isnan(warpA_Z))=1;
    %error(error_mask==1)=[];
    %errorEst = sum(error.*error)
    %numValid = length(error)
    
    acc_hessian = zeros(6,6);
    depth_grad = zeros(6,1);
    intensity_grad = zeros(6,1);
    numEqs = 0;
    sqError = 0;
    PtIdx = [];
    if 0
        for ptIdx=1:numPts
            if (~isnan(pts(3,ptIdx)))
                pt = T*pts(1:4,ptIdx);
                intensity = imgI_a(ptIdx);
                pt2d = imageB.project(pt);
                if (pt2d(1)>= 1 && pt2d(1) < cols && ...
                        pt2d(2)>= 1 && pt2d(2) < rows)
                    x0y0 = floor(pt2d(1:2));
                    x1y1 = floor(pt2d(1:2)+1);
                    x0y1 = [x0y0(1); x1y1(2)];
                    x1y0 = [x1y1(1); x0y0(2)];
                    corners = [x0y0 x1y0 x0y1 x1y1];
                    alpha_xy = [pt2d(1)-x0y0(1); pt2d(2)-x0y0(2)];
                    idxs = sub2ind(size(imgZ_a), corners(2,:), corners(1,:));
                    fc = [imgZ_b(idxs)' gradZ_x(idxs)' gradZ_y(idxs)' imgI_b(idxs)' gradI_x(idxs)' gradI_y(idxs)'];
                    fc_interp = (fc(1,:)*(1-alpha_xy(1)) + fc(2,:)*alpha_xy(1)) * (1-alpha_xy(2)) + (fc(3,:)*(1-alpha_xy(1)) + fc(4,:)*alpha_xy(1)) * alpha_xy(2);
                    if (~isnan(sum(fc_interp)))
                        %ipt3 = imageB.getPoint3d(pt2d, fc_interp(1));
                        %[Jx, Jy, Jz] = imageB.getJacobian(ipt3);
                        [Jx, Jy, Jz] = imageB.getJacobian(pt);
                        residual_Z = fc_interp(1) - pt(3);
                        residual_I = fc_interp(4) - intensity;
                        gradZ_interp = fc_interp(2:3)*[Jx; Jy] - Jz;
                        gradI_interp = fc_interp(5:6)*[Jx; Jy];
                        
                        depth_grad = depth_grad + residual_Z*gradZ_interp';
                        intensity_grad = intensity_grad + residual_I*gradI_interp';
                        acc_hessian = acc_hessian + gradZ_interp'*gradZ_interp;
                        acc_hessian = acc_hessian + gradI_interp'*gradI_interp;
                        sqError = sqError + residual_Z*residual_Z;
                        sqError = sqError + residual_I*residual_I;
                        
                        numEqs = numEqs + 1;
                        PtIdx = [PtIdx; ptIdx];
                    end
                else
                    %disp('error');
                end
            end
        end
    end
    intensityMat = reshape(imgI_a, rows, cols);
    ptsMat = reshape(pts(3,:), rows, cols);
    ptList = T*pts(1:4,:);
    pt2dList = project(K, ptList);
    
    [xu, yu] = meshgrid(1:cols ,1:rows);
    Z_b = interp2(xu,yu,double(reshape(imgZ_b, rows, cols)),reshape(pt2dList(1,:), rows, cols),reshape(pt2dList(2,:), rows, cols));
    Z_x = interp2(xu,yu,double(gradZ_x),reshape(pt2dList(1,:), rows, cols),reshape(pt2dList(2,:), rows, cols));
    Z_y = interp2(xu,yu,double(gradZ_y),reshape(pt2dList(1,:), rows, cols),reshape(pt2dList(2,:), rows, cols));
    I_b = interp2(xu,yu,double(reshape(imgI_b, rows, cols)),reshape(pt2dList(1,:), rows, cols),reshape(pt2dList(2,:), rows, cols));
    I_x = interp2(xu,yu,double((gradI_x)),reshape(pt2dList(1,:), rows, cols),reshape(pt2dList(2,:), rows, cols));
    I_y = interp2(xu,yu,double((gradI_y)),reshape(pt2dList(1,:), rows, cols),reshape(pt2dList(2,:), rows, cols));
    
    fc_interp_mat = [Z_b(:), Z_x(:), Z_y(:), I_b(:), I_x(:), I_y(:)];
    
    fc_interp_mat_sum = sum(fc_interp_mat');
    
    validFlag = (~isnan(pts(3,:)) & ~isnan(fc_interp_mat_sum) & (pt2dList(1,:)>= 1 & pt2dList(1,:) < cols & pt2dList(2,:)>= 1 & pt2dList(2,:) < rows));
    
    unValidFlag = ~validFlag;
    fc_interp_mat(unValidFlag,:) = nan;
    %                 fc_interp_mat_sum = sum(fc_interp_mat');
    
    
    [JxMat, JyMat, JzMat] = getJacobian(K, ptList);
    
    
    
    
    residual_Z_mat = fc_interp_mat(:,1) - ptList(3,:)';
    residual_I_mat = fc_interp_mat(:,4) - imgI_a';
    
    gradZ_interp_mat = [dot(fc_interp_mat(:,2:3)', [JxMat(:,1) JyMat(:,1)]');...
        dot(fc_interp_mat(:,2:3)', [JxMat(:,2) JyMat(:,2)]');...
        dot(fc_interp_mat(:,2:3)', [JxMat(:,3) JyMat(:,3)]');...
        dot(fc_interp_mat(:,2:3)', [JxMat(:,4) JyMat(:,4)]');...
        dot(fc_interp_mat(:,2:3)', [JxMat(:,5) JyMat(:,5)]');...
        dot(fc_interp_mat(:,2:3)', [JxMat(:,6) JyMat(:,6)]')]' - JzMat;
    
    gradI_interp_mat = [dot(fc_interp_mat(:,5:6)', [JxMat(:,1) JyMat(:,1)]');...
        dot(fc_interp_mat(:,5:6)', [JxMat(:,2) JyMat(:,2)]');...
        dot(fc_interp_mat(:,5:6)', [JxMat(:,3) JyMat(:,3)]');...
        dot(fc_interp_mat(:,5:6)', [JxMat(:,4) JyMat(:,4)]');...
        dot(fc_interp_mat(:,5:6)', [JxMat(:,5) JyMat(:,5)]');...
        dot(fc_interp_mat(:,5:6)', [JxMat(:,6) JyMat(:,6)]')]';
    
    
    
    depth_grad_mat = (repmat(residual_Z_mat, 1, 6).*gradZ_interp_mat);
    depth_grad_mat(isnan(depth_grad_mat)) = 0;
    depth_grad_mat(unValidFlag,:) = 0;
    depth_grad_mat_sum = sum(depth_grad_mat)';
    
    
    intensity_grad_mat = repmat(residual_I_mat, 1, 6).*gradI_interp_mat;
    intensity_grad_mat(isnan(intensity_grad_mat)) = 0;
    intensity_grad_mat(unValidFlag,:) = 0;
    intensity_grad_mat_sum = sum(intensity_grad_mat)';
    
    
    
    
    gradZ_interp_mat1 = repmat(gradZ_interp_mat,6,1);
    gradZ_interp_mat1 = reshape(gradZ_interp_mat1, size(gradZ_interp_mat,1),6,6);
    gradZ_interp_mat1 = permute(gradZ_interp_mat1, [3 2 1]);
    
    gradZ_interp_mat2 = permute(gradZ_interp_mat1, [2 1 3]);
    
    
    
    gradI_interp_mat1 = repmat(gradI_interp_mat, 6, 1);
    gradI_interp_mat1 = reshape(gradI_interp_mat1, size(gradI_interp_mat,1),6,6);
    gradI_interp_mat1 = permute(gradI_interp_mat1, [3 2 1]);
    
    gradI_interp_mat2 = permute(gradI_interp_mat1, [2 1 3]);
    
    
    
    acc_hessian_mat = gradZ_interp_mat1.*gradZ_interp_mat2 + gradI_interp_mat1.*gradI_interp_mat2;
    acc_hessian_mat(isnan(acc_hessian_mat)) = 0;
    acc_hessian_mat(:,:,unValidFlag) = 0;
    
    acc_hessian_mat_sum = sum(acc_hessian_mat,3);
    
    
    sqError_mat = residual_Z_mat.^2 + residual_I_mat.^2;
    
    numEqs_mat = sum(~isnan(sqError_mat));
    
    sqError_mat(isnan(sqError_mat)) = 0;
    sqError_mat(unValidFlag) = 0;
    sqError_mat_sum = sum(sqError_mat);
    
    
    if 0
        acc_hessian_mat_sum - acc_hessian
        sqError_mat_sum - sqError
        intensity_grad_mat_sum - intensity_grad
        depth_grad_mat_sum - depth_grad
    end
    
    acc_hessian = acc_hessian_mat_sum;
    sqError = sqError_mat_sum;
    intensity_grad = intensity_grad_mat_sum;
    depth_grad = depth_grad_mat_sum;
    numEqs = sum(validFlag);
    
    %numEqs
    sqError
    %sqError/numEqs
    if (sqError > prev_sqError || numEqs < 6)
        break;
    end
    prev_sqError = sqError;
    acc_hessian = acc_hessian;% + numEqs*numEqs*(0.5)*eye(6);
    %acc_hessian = acc_hessian - 0*numEqs*numEqs*(0.1^iter)*eye(6);
    delta_a_k = inv(acc_hessian)*depth_grad;
    %delta_a_k
    
    delta_T = rgbd_odom.parametersToTransform(delta_a_k');
    % store previous R,t in case error increases
    prev_T = T;
    T=inv(delta_T)*T;
end
iter;
% if error went up get minimum as previous R and t
if (iter < MAXITER)
    T = prev_T;
end
% a_k = transformToParameters(T);

if 0
    warpA = imageA.warpImage(inv(T));
    figure,imshowpair(warpA.imgRGB, imageB.imgRGB);
end

end
function p3d = getPointCloud(imgZ, width, height,K)
P2_to_P3 = inv([K(1,1) 0 0 K(1,3); 0 K(2,2) 0 K(2,3);0 0 1 0;0 0 0 1]);
numPixels = width*height;
[xi, yi]=meshgrid( 1:width, 1:height);
p3d = [reshape( xi, 1, numPixels);
    reshape( yi, 1, numPixels);
    reshape( imgZ, 1, numPixels);
    ones( 1, numPixels)];
p3d = P2_to_P3*p3d;
p3d(1:2,:) = p3d(1:2,:).*p3d(3,:);
end
function T = parametersToTransform(a_k)
theta_k = sqrt(a_k(4:6)*a_k(4:6)');
if (theta_k > 1e-10)
    n_k = a_k(4:6)/theta_k;
    n_cross_k = [0 -n_k(3) n_k(2); n_k(3) 0 -n_k(1); -n_k(2) n_k(1) 0];
    R = cos(theta_k)*eye(3) + (1-cos(theta_k))*(n_k'*n_k) + sin(theta_k)*n_cross_k;
    V = eye(3) + ((1-cos(theta_k))/theta_k)*n_cross_k + (1 - sin(theta_k)/theta_k)*(n_k'*n_k);
    t = V*a_k(1:3)';
else
    R = eye(3);
    t = a_k(1:3)';
end

T=[R t; 0 0 0 1];
end
function pts2d = project(K, pts3d)
P3_to_P2 = ([K(1,1) 0 0 K(1,3); 0 K(2,2) 0 K(2,3);0 0 1 0;0 0 0 1]);
pts3d(1:2,:) = pts3d(1:2,:)./pts3d(3,:);
pts2d = round(P3_to_P2*pts3d, 10);
end
function [Jx_w, Jy_w, Jz_w] = getJacobian(K, pt)

f.x = K(1,1);
f.y = K(2,2);

ONEVEC = ones(1, size(pt,2));
ZEROVEC = zeros(1, size(pt,2));

Jx_w = (f.x./pt(3,:) .* [ONEVEC; ...
    ZEROVEC; ...
    -pt(1,:)./pt(3,:); ...
    -pt(1,:).*pt(2,:)./pt(3,:); ...
    pt(3,:) + (pt(1,:).*pt(1,:)./pt(3,:)); ...
    -pt(2,:)])';
Jy_w = (f.y./pt(3,:) .* [ZEROVEC; ...
    ONEVEC; ...
    -pt(2,:)./pt(3,:); ...
    -pt(3,:) + (pt(2,:).*pt(2,:)./pt(3,:)); ...
    pt(1,:).*pt(2,:)./pt(3,:); ...
    pt(1,:)])';
Jz_w = [ZEROVEC; ZEROVEC; ONEVEC; pt(2,:); -pt(1,:); ZEROVEC]';
end
function a_k = transformToParameters(T)
R = T(1:3,1:3);
t = T(1:3,4);
theta = acos(0.5*(trace(R)-1));

if (theta > 1e-10)
    a_k(4:6) = theta*(1/(2*sin(theta)))*[R(3,2)-R(2,3), R(1,3)-R(3,1), R(2,1)-R(1,2)];
    w = a_k(4:6);
    w_cross = [0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0];
    V_inv = eye(3) - 0.5*w_cross + ((1 - theta*sin(theta)/(2*(1-cos(theta))))/theta^2)*w'*w;
    a_k(1:3) = (V_inv*t)';
else
    a_k(4:6) = [0 0 0];
    a_k(1:3) = t;
end

end