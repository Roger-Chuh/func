function SimPixelTrace()


% T = [rodrigues([0.9 0.2 0.3]) [100 200 200]';0 0 0 1];
angle = 100;
trans_scale = 0.3;
T = [roty(angle) trans_scale.*[1000 -1000 1000]';0 0 0 1];
K = [200 0 320;0 200 240;0 0 1];

width = 640;
height = 480;

depth_scale = 3;
xyz = depth_scale*[-1500 200 1000];

idepth = 1/xyz(3);

RKi = T(1:3,1:3)*inv(K);
t = T(1:3,4);

host_pix_ = pflat(K*xyz');
host_pix = host_pix_(1:2)';


target_pix_check = RKi * [host_pix 1]' + t*idepth;
target_pix_check = target_pix_check./target_pix_check(3);
target_pix = K*target_pix_check;

F = inv(K')*SkewSymMat(T(1:3,4)) * T(1:3,1:3)*inv(K);
err = target_pix'*F*host_pix_;
pole1 = null(F);
pole_in_1 = pole1./pole1(3);
pole2 = null(F');
pole_in_2 = pole2./pole2(3);

pole_in_2_check = pflat(K*T(1:3,4));

epiline_host_in_target = F*host_pix_; %这是target帧上极线的方向，但是极限的方向是用host的pixel坐标计算得到
epiline_host_in_target = epiline_host_in_target./norm(epiline_host_in_target(1:2));
slop_target = -epiline_host_in_target(1)/epiline_host_in_target(2);


%正的纯平移
target_pix_t_positive = inv(K) * [host_pix 1]' + t*idepth;
target_pix_t_positive = pflat(K*target_pix_t_positive);
%负的纯平移
target_pix_t_negative = inv(K) * [host_pix 1]' - t*idepth;
target_pix_t_negative = pflat(K*target_pix_t_negative);

%正的旋转+平移
target_pix_rt_positive = RKi * [host_pix 1]' + t*idepth;
target_pix_rt_positive = pflat(K*target_pix_rt_positive);
%负的旋转+平移
target_pix_rt_negative = RKi * [host_pix 1]' - t*idepth;
target_pix_rt_negative = pflat(K*target_pix_rt_negative);


[ptIcs, tgtPt3d] = TransformAndProject(xyz, K, T(1:3,1:3), T(1:3,4));
van_pt = pflat(K*T(1:3,4));

len1 = norm(target_pix_rt_positive - target_pix_rt_negative);
len2 = norm(target_pix_t_positive - target_pix_t_negative);

dist1 = norm(target_pix_rt_positive(1:2)' - host_pix)^2 + norm(target_pix_rt_negative(1:2)' - host_pix)^2;
dist2 = norm(target_pix_t_positive(1:2)' - host_pix)^2 + norm(target_pix_t_negative(1:2)' - host_pix)^2;

figure,imshow(ones(480, 640));hold on;plot(host_pix(:,1), host_pix(:,2),'.r');plot(pole_in_2_check(1),pole_in_2_check(2),'or');
plot(target_pix_rt_positive(1), target_pix_rt_positive(2),'.g');
plot(van_pt(1), van_pt(2),'xy');
plot(target_pix_rt_negative(1), target_pix_rt_negative(2),'.b');
plot(target_pix_t_positive(1), target_pix_t_positive(2),'.m');
plot(target_pix_t_negative(1), target_pix_t_negative(2),'.c');
title(sprintf('[len1 / len2]: [%f %f], ratio = %f, depth-scale: %f, angle: %f, trans-scale: %f\n[dist1 + dist2]: [%f + %f] = %f, value: %f', len1, len2, len1/len2, depth_scale, angle, trans_scale, dist1, dist2,dist1+dist2, sqrt((dist1+dist2)/2)/(width+height)));

dir1 = target_pix_rt_positive - target_pix_rt_negative;
dir1 = dir1(1:2);
dir1 = dir1./norm(dir1);

dir2 = van_pt - target_pix_rt_positive;
dir2 = dir2(1:2);
dir2 = dir2./norm(dir2);
[dir1' - dir2' dir1' (slop_target-dir1(2)/dir1(1))]
end