function SymSphere()
  close all



colorStack = {'k';'m';'c'};
colorStack = repmat(colorStack,100,1);

L2R = [rodrigues([0.02; 0.01; -0.003]) [-50.1 0.01 -0.1]';0 0 0 1];

markerSize = 800; 150; 100; 200;
markerRow = 3; 15; 4; 7; 15;20; 7;
markerCol = 3; 15; 4; 7; 15;20; 7;

[xMat0, yMat0] = meshgrid(markerSize : markerSize : (markerCol - 0) * markerSize, markerSize : markerSize : (markerRow - 0) * markerSize);
xMat0 = xMat0';
yMat0 = yMat0';
marker1 = [xMat0(:) yMat0(:) zeros(markerRow*markerCol,1)];

Cam2MarkerRot = rodrigues([-0.01 1.0 -0.02] + 0.05*(rand(1,3)-0.5));
Cam2MarkerRot = rodrigues([-0.01 1.0 -0.02] + 0.00*(rand(1,3)-0.5));
Cam2MarkerTrans = [-100 1100 -1000]';
Cam2Marker = [Cam2MarkerRot Cam2MarkerTrans; 0 0 0 1];
Marker2CamL = inv(Cam2Marker);
Marker2CamR = L2R*Marker2CamL;


width = 640; height = 480;
K = eye(3);
K(1,1) = 236.8820281-0;
K(2,2) = 236.7437341-0;
K(1,3) = 325.0065396;
K(2,3) = 241.9557923;
dist = [0.2124214875, -0.1719834148, 0.06211290885, -0.009899091649]';
% dist = [0.06314965524 -0.07006363688 0.01760113423 0.0001556135281]';
% dist = [0.46314965524 -0.47006363688 0.1760113423 0.001556135281]';


T00 = eye(4);
T01 = [ -7.4656856642149472e-01, -7.0503228692542969e-02, -6.6156229515842346e-01, -3.5247847596384795e-02;
    -6.2170843696788058e-01, -2.8011666159265608e-01,7.3144601666913212e-01, 4.1311291893071443e-02;
    -2.3688392734482649e-01, 9.5737346455921379e-01,1.6529384236437439e-01, -3.4150037619299412e-02;
    0., 0., 0., 1. ];

T02 = [ -4.3634354777464124e-01, 4.7520621068004289e-02, -8.9852440083102414e-01, -6.8945219190976056e-02;
    8.9764389729039329e-01, 9.1761591914208740e-02, -4.3106292337277669e-01, -1.1696462645325466e-01;
    6.1965651555951469e-02, -9.9464647027106146e-01, -8.2696174062438135e-02, -6.6681057254321219e-02;
    0., 0., 0., 1. ];

T03 = [ 9.2370408189253994e-01, 3.7671281461209621e-01, -6.9700964140337329e-02, -2.1413263321725769e-02;
    -3.7601452313871109e-01, 8.5662781853225878e-01, -3.5327306563822758e-01, -1.0905148943810981e-01;
    -7.3374706022088163e-02, 3.5252834754625934e-01, 9.3292010198755815e-01, -2.2677364423121626e-02;
    0., 0., 0., 1. ];

T01(1:3,4) = 1000.*T01(1:3,4);
T02(1:3,4) = 1000.*T02(1:3,4);
T03(1:3,4) = 1000.*T03(1:3,4);

[xMat, yMat] = meshgrid(1 : width, 1 : height);
pix = [xMat(:) yMat(:)];

pt3d = unprojectKB8( K(1,1),  K(2,2),  K(1,3),  K(2,3),  dist(1),  dist(2),  dist(3),  dist(4), pix);
[~,norm_3d] = NormalizeVector(pt3d);
[pt2d, inImageFlag] = projectKB8(pextend(pt3d')', K, dist, eye(3), [0;0;0]);


[pt2dL00, pt3dL] = TransformAndProject(marker1, K, Marker2CamL(1:3,1:3), Marker2CamL(1:3,4));

[pt2d_show, inImageFlag_show] = projectKB8(pextend(marker1')', K, dist, Marker2CamL(1:3,1:3), Marker2CamL(1:3,4));
[pt2d_showR, inImageFlag_showR] = projectKB8(pextend(marker1')', K, dist, Marker2CamR(1:3,1:3), Marker2CamR(1:3,4));
figure(20),imshow(zeros(480, 640));hold on; plot(pt2d_show(:,1), pt2d_show(:,2), 'or','MarkerSize',5,'LineWidth',5);
plot(pt2d_showR(:,1), pt2d_showR(:,2), 'ob','MarkerSize',5,'LineWidth',5);
plot(pt2d_show(1,1), pt2d_show(1,2), 'xg','MarkerSize',10,'LineWidth',10);

% bearing vector in left camera
pt3d_show = unprojectKB8( K(1,1),  K(2,2),  K(1,3),  K(2,3),  dist(1),  dist(2),  dist(3),  dist(4), pt2d_show);
pt3d_showR = unprojectKB8( K(1,1),  K(2,2),  K(1,3),  K(2,3),  dist(1),  dist(2),  dist(3),  dist(4), pt2d_showR);

bearing_L_in_R = (L2R * pextend(pt3d_show'))';
oc_L_in_R = [L2R(1:3,4)' 1];



MarkerId = 1 : size(marker1,1);
MarkerId = reshape(MarkerId, markerRow, markerRow);

plane_norm_err = [];
for i = 1 : markerRow
     if 1
        markerId = MarkerId(:,i)';
    else
        markerId = MarkerId(i,:);
     end
    pt3d_show_left = pt3d_show(markerId,:);
    pt3d_show_list = bearing_L_in_R(markerId,:);
    pt3d_showR_list = pt3d_showR(markerId,:);
    for j = 1 : size(pt3d_show_list,1)
        [~,~,plane_norm_] = svd([pt3d_show_list(j,:); oc_L_in_R; [0 0 0 1]],0);
        plane_norm_ = plane_norm_(:,4);
        plane_norm_ = plane_norm_./norm(plane_norm_(1:3));
        plane_norm_ = sign(plane_norm_(2)).*plane_norm_;
        
        plane_norm_check = cross(oc_L_in_R(1:3),pt3d_show_list(j,1:3));
        plane_norm_check_check = SkewSymMat(oc_L_in_R(1:3))*pt3d_show_list(j,1:3)';
        
        plane_norm_check = plane_norm_check./norm(plane_norm_check);       
        plane_norm_check_check = plane_norm_check_check./norm(plane_norm_check_check);
        
        pt3dR = pt3d_showR_list(j,:);
        pt_in_plane_err = dot(plane_norm_, [pt3dR 1]);
        
        % formulation
        plane_norm_compute = (SkewSymMat(L2R(1:3,4)) * L2R(1:3,1:3)*pt3d_show_left(j,:)')';
        plane_norm_compute = plane_norm_compute./norm(plane_norm_compute);
        
        plane_norm_err = [plane_norm_err; (plane_norm_(1:3)' - plane_norm_check) pt_in_plane_err (plane_norm_check - plane_norm_check_check') (plane_norm_compute -  plane_norm_check_check') (plane_norm_compute - plane_norm_(1:3)')];
    end

end
    
[~,err] = NormalizeVector(pt2d - pix);
validId = abs(norm_3d-1) < 0.0000001;
mask = zeros(height, width);
mask(validId) = 1;
figure(30),imshow(mask)
% figure,plot(err(validId))
figure(10),hold on;grid on;pcshow(pt3d(validId,:));hold on;plot3(pt3d_show(:,1), pt3d_show(:,2), pt3d_show(:,3), 'xr','MarkerSize',10,'LineWidth',10);
plot3(pt3dL(:,1)./5000, pt3dL(:,2)./5000, pt3dL(:,3)./5000, 'xb','MarkerSize',10,'LineWidth',10);
plot3(0,0,0, 'xk','MarkerSize',15,'LineWidth',15);





% figure(40);hold on;pcshow([xyz0],[1 0 0]);pcshow([xyz1],[0 1 0]);pcshow([xyz2],[0 0 1]);pcshow([xyz3],[1 1 0]);
% plot3(0,0,0, 'xk','MarkerSize',15,'LineWidth',15);


distance = 3000; 500; 300; 1000;



% idVec = [0 0 0];

%% base 0
idVec = [1 2 3];
[mat_0_in_1, mat_0_in_2, mat_0_in_3, matAll_to_0] = CalcRegion(distance, K, dist, pt3d, validId,T01,T02,T03, idVec);
% figure,imshow([mat_0_in_1, mat_0_in_2, mat_0_in_3], []);

%% base 1
idVec = [0 2 3];
T10 = inv(T01);
T12 = inv(T01) * T02;
T13 = inv(T01) * T03;
[mat_1_in_0, mat_1_in_2, mat_1_in_3, matAll_to_1] = CalcRegion(distance,K, dist, pt3d, validId,T10,T12,T13, idVec);
% figure,imshow([mat_1_in_0 mat_1_in_2 mat_1_in_3], []);

%% base 2
idVec = [0 1 3];
T20 = inv(T02);
T21 = inv(T02) * T01;
T23 = inv(T02) * T03;
[mat_2_in_0, mat_2_in_1, mat_2_in_3, matAll_to_2] = CalcRegion(distance,K, dist, pt3d, validId,T20,T21,T23, idVec);
% figure,imshow([mat_2_in_0, mat_2_in_1, mat_2_in_3], []);

%% base 3
idVec = [0 1 2];
T30 = inv(T03);
T31 = inv(T03) * T01;
T32 = inv(T03) * T02;
[mat_3_in_0, mat_3_in_1, mat_3_in_2, matAll_to_3] = CalcRegion(distance,K, dist, pt3d, validId,T30,T31,T32, idVec);
% figure,imshow([mat_3_in_0, mat_3_in_1, mat_3_in_2], []);

comb = [[mat_0_in_1, mat_0_in_2, mat_0_in_3];[mat_1_in_0 mat_1_in_2 mat_1_in_3];[mat_2_in_0, mat_2_in_1, mat_2_in_3];[mat_3_in_0, mat_3_in_1, mat_3_in_2]];

if 1
    figure,imagesc(comb);axis equal;colorbar;
    
    figure,subplot(4,4,1);imshow(mat_0_in_1, []);title('1 to 0');
    subplot(4,4,2);imshow(mat_0_in_2, []);title('2 to 0');
    subplot(4,4,3);imshow(mat_0_in_3, []);title('3 to 0');
    subplot(4,4,4);imshow(matAll_to_0, []);title('all cams seen by cam0');
    
    subplot(4,4,5);imshow(mat_1_in_0, []);title('0 to 1');
    subplot(4,4,6);imshow(mat_1_in_2, []);title('2 to 1');
    subplot(4,4,7);imshow(mat_1_in_3, []);title('3 to 1');
    subplot(4,4,8);imshow(matAll_to_1, []);title('all cams seen by cam1');
    
    subplot(4,4,9);imshow(mat_2_in_0, []);title('0 to 2');
    subplot(4,4,10);imshow(mat_2_in_1, []);title('1 to 2');
    subplot(4,4,11);imshow(mat_2_in_3, []);title('3 to 2');
    subplot(4,4,12);imshow(matAll_to_2, []);title('all cams seen by cam2');
    
    subplot(4,4,13);imshow(mat_3_in_0, []);title('0 to 3');
    subplot(4,4,14);imshow(mat_3_in_1, []);title('1 to 3');
    subplot(4,4,15);imshow(mat_3_in_2, []);title('2 to 3');
    subplot(4,4,16);imshow(matAll_to_3, []);title('all cams seen by cam3');
end


for i = 1 :1: markerRow
    
    %     markerId = (i*markerRow - markerRow+1): i*markerRow;
    if 1
        markerId = MarkerId(:,i)';
    else
        markerId = MarkerId(i,:);
    end
    pt3d_show_list = pt3d_show(markerId,:);
    
    line_stack = [];
    for j = 1 : 1: size(pt3d_show_list,1)
        line_stack = [line_stack; [0 0 0];pt3d_show_list(j,:) ];
    end
    line_stack = [line_stack; [0 0 0]];
    figure(10), plot3(pt3d_show_list(:,1), pt3d_show_list(:,2), pt3d_show_list(:,3), '-sg','MarkerSize',5,'LineWidth',5);
    plot3(pt3dL(markerId,1)./5000, pt3dL(markerId,2)./5000, pt3dL(markerId,3)./5000, '-sb','MarkerSize',5,'LineWidth',5);
    line(line_stack(:,1),line_stack(:,2),line_stack(:,3),'Color',colorStack{i},'LineWidth', 2,'MarkerSize', 2);
    
    [~,~,plane_norm] = svd(pextend(pt3d_show_list(:,1:3)')',0);
    plane_norm = plane_norm(:,4);
    plane_norm = plane_norm./norm(plane_norm(1:3));
    plane_norm = sign(plane_norm(2)).*plane_norm;
    
    plotQuiver(plane_norm(1:3)',colorStack{i});
    %     plotQuiver(plane_norm(1:3)');
    hori_van_dir_err(i,1) = sum(abs(dot(repmat(plane_norm,1,size(pt3d_show_list,1)),pextend(pt3d_show_list(:,1:3)'))));
    plane_norms(i,:) = plane_norm(1:3)';
end


[~,~,plane_norm_norm] = svd(pextend(plane_norms(:,1:3)')',0);
plane_norm_norm = plane_norm_norm(:,4);
plane_norm_norm = plane_norm_norm./norm(plane_norm_norm(1:3));
plane_norms_dir_err = abs(dot(repmat(plane_norm_norm,1,size(plane_norms,1)),pextend(plane_norms(:,1:3)')));

figure(10),plotQuiver(plane_norm_norm(1:3)','r');

Marker2CamL(1:3,1:3)
plane_norm_norm


% plane_norm1 = [0.11 0.9, 0.1];
% plane_norm1 = plane_norm1./norm(plane_norm1);
% plane1 = [plane_norm1 0];


end
function [pt2d, inImageFlag] = projectKB8(pt3d, K, dis, R, t)
pt3d = ([R t; 0 0 0 1] * pt3d')';
pt3d = pt3d(:,1:3);

% inImageFlag = pt3d(:,3) > 0.005;

fx = K(1,1);
fy = K(2,2);
cx = K(1,3);
cy = K(2,3);
k1 = dis(1);
k2 = dis(2);
k3 = dis(3);
k4 = dis(4);

x = pt3d(:,1);
y = pt3d(:,2);
z = pt3d(:,3);

r2 = x .* x + y .* y;
r = sqrt(r2);

theta = atan2(r, z);
theta2 = theta .* theta;
theta4 = theta2 .* theta2;
theta6 = theta4 .* theta2;
theta8 = theta6 .* theta2;

r_theta = theta .* (1 + k1 .* theta2 + k2 .* theta4 + k3 .* theta6 + k4 .* theta8);
if r > 1e-8
    norm_inv = 1./r;
else
    norm_inv = ones(length(r), 1);
end
% norm_inv = r > 1e-8 ? double(1.0) / r : 1;

mx = r_theta .* x .* norm_inv;
my = r_theta .* y .* norm_inv;

pt2d(:,1) = fx .* mx + cx;
pt2d(:,2) = fy .* my + cy;



inImageFlag = find(pt2d(:,1) > 1 & pt2d(:,1) < 640-1 & pt2d(:,2) > 1 & pt2d(:,2) < 480-1 & pt3d(:,3)> 10);

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
    d_func_d_theta = d_func_d_theta + (7) .* k3;
    d_func_d_theta = d_func_d_theta .* theta2;
    d_func_d_theta = d_func_d_theta + (5) .* k2;
    d_func_d_theta = d_func_d_theta .* theta2;
    d_func_d_theta = d_func_d_theta + (3) .* k1;
    d_func_d_theta = d_func_d_theta .* theta2;
    d_func_d_theta = d_func_d_theta + (1);
    
    
    %theta = theta + (r_theta - func) / d_func_d_theta;
    theta_fix = (r_theta - func) ./ d_func_d_theta;
    theta = theta+ theta_fix;
    if 0 % (abs(theta_fix) < 1e-8)
        break;
    end
end
d_func_d_theta1 = d_func_d_theta;
end
function p3d = unprojectKB8( fx,  fy,  cx,  cy,  k1,  k2,  k3,  k4, proj)


mx = (proj(:,1) - cx) ./ fx;
my = (proj(:,2) - cy) ./ fy;




theta = 0;
sin_theta = 0;
cos_theta = 1;
thetad = sqrt(mx .* mx + my .* my);

thetad = min(max(thetad,-3.141592653/2),3.141592653/2);
scaling = 1.0;
d_func_d_theta = 0;
if (thetad > 1e-8)
    [theta,d_func_d_theta1] = solveTheta(k1, k2, k3, k4, thetad, d_func_d_theta,3);
    d_func_d_theta = d_func_d_theta1;
    sin_theta = sin(theta);
    cos_theta = cos(theta);
    scaling = sin_theta ./ thetad;
end

p3d = [mx .* scaling my .* scaling cos_theta];


end
function [BW10, BW20, BW30, matAll] = CalcRegion(distance, K, dist, pt3d, validId,T01,T02,T03, idVec)
% function [mat1, mat2, mat3] = CalcRegion(K, dist, pt3d, validId,T01,T02,T03, idVec)
% 10 12 13
validId = find(validId);
[xMat, yMat] = meshgrid(1:640, 1:480);
pix = [xMat(:) yMat(:)];
mat1 = -1.*ones(480, 640);
mat2 = mat1;
mat3 = mat1;
xyz0 = distance.*pt3d(validId,:);
xyz1 = ((T01)*pextend(xyz0'))'; xyz1 = xyz1(:,1:3);
xyz2 = ((T02)*pextend(xyz0'))'; xyz2 = xyz2(:,1:3);
xyz3 = ((T03)*pextend(xyz0'))'; xyz3 = xyz3(:,1:3);

[pt2d_1_in_0, validFlag10] = projectKB8(pextend(xyz1')', K, dist, eye(3), [0 0 0]');
BW10 = CalcValidPix(pt2d_1_in_0, validFlag10);
mat1(validId(validFlag10)) = idVec(1);
BW10(BW10 > 0) = idVec(1);
% figure(51),imshow(zeros(480,640));hold on;plot(pt2d_1_in_0(validFlag10,1), pt2d_1_in_0(validFlag10,2),'.r');title('1 in 0');

[pt2d_2_in_0, validFlag20] = projectKB8(pextend(xyz2')', K, dist, eye(3), [0 0 0]');
BW20 = CalcValidPix(pt2d_2_in_0, validFlag20);
mat2(validId(validFlag20)) = idVec(2);
BW20(BW20 > 0) = idVec(2);
% figure(52),imshow(zeros(480,640));hold on;plot(pt2d_2_in_0(validFlag20,1), pt2d_2_in_0(validFlag20,2),'.r');title('2 in 0');

[pt2d_3_in_0, validFlag30] = projectKB8(pextend(xyz3')', K, dist, eye(3), [0 0 0]');
BW30 = CalcValidPix(pt2d_3_in_0, validFlag30);
mat3(validId(validFlag30)) = idVec(3);
BW30(BW30 > 0) = idVec(3);
% figure(50),imshow(zeros(480,640));hold on;plot(pt2d_3_in_0(validFlag30,1), pt2d_3_in_0(validFlag30,2),'.r');title('3 in 0');

matAll = BW10 + BW20 + BW30;
end
function BW2 = CalcValidPix(pt2d_2_in_0, validFlag20)
pixels_hit = unique(round(pt2d_2_in_0(validFlag20,:)),'rows');
ind = sub2ind([480 640], pixels_hit(:,2),pixels_hit(:,1));
validMask = zeros(480, 640);
validMask(ind) = 1;
se = strel('square',3);
BW2 = imdilate(validMask,se);
BW2(BW2==0) = -1;
end