function testEpipolarFactorForCalibration()

pose{1,1} = [rodrigues([0 0 0]) [0 0 0]';0 0 0 1]; % host frame c0
% pose{2,1} = [rodrigues([0.1 0.02 0.03]) [10 20 30]';0 0 0 1]; % target frame c0
%% delta pose is identity, we have covisibility between
pose{2,1} = [rodrigues([0 0 0]) [0 0 0]';0 0 0 1]; % target frame c0

pose{3,1} = [rodrigues([0.1 -0.02 0.03]) [10 20 -30]';0 0 0 1]; % taregt frame c1
pose{4,1} = [rodrigues([-0.1 0.02 0.03]) [-10 20 30]';0 0 0 1]; % taregt frame c2
pose{5,1} = [rodrigues([0.1 -0.02 -0.03]) [-10 -20 30]';0 0 0 1]; % taregt frame c3

T_host = pose{1,1};
T_target = pose{2,1};
T_bc{1,1} = eye(4);
T_bc{2,1} = inv(pose{2}) * pose{3,1};
T_bc{3,1} = inv(pose{2}) * pose{4,1};
T_bc{4,1} = inv(pose{2}) * pose{5,1};

XYZ = [10 20 1300; 12 22 1302;13 23 1303;14 24 1304;15 25 1305;16 26 1306;17 27 1307;18 28 1308;19 29 1309;161 261 13106;126 226 12306;136 236 13306];
% XYZ2 = XYZ + 2;
% XYZ = [XYZ; XYZ2];
% XYZ2 = XYZ + 2;
% XYZ = [XYZ; XYZ2];
% XYZ2 = XYZ + 2;
% XYZ = [XYZ; XYZ2];
% XYZ2 = XYZ + 2;
% XYZ = [XYZ; XYZ2];
% XYZ2 = XYZ + 2;
% XYZ = [XYZ; XYZ2];
% XYZ = XYZ(1:1,:);
[host_bearings, dist] = NormalizeVector(XYZ);
lambdas = 1 ./ dist;


K = [400 0 320; 0 400 240;0 0 1];

poses = {};


Tbw{1, 1} = [rodrigues([0.001 0.002 0.003]) [10 20 30]';0 0 0 1];
Tbw{2, 1} = [rodrigues([0.001 -0.002 0.003]) [10 20 2]';0 0 0 1];
Tbw{3, 1} = [rodrigues([0.001 -0.002 -0.003]) [10 2 30]';0 0 0 1];
Tbw{4, 1} = [rodrigues([-0.001 0.002 -0.003]) [10 11 30]';0 0 0 1];
Tbw{5, 1} = [rodrigues([-0.001 0.002 0.003]) [10 -20 30]';0 0 0 1];

VM = {};
for k = 1 : length(Tbw)
    VM{k,1} = {};
    for i = 1 : length(T_bc)
        Tcw = inv(T_bc{i}) * Tbw{k};
        [ptIcs, tgtPt3d] = TransformAndProject(XYZ, K, Tcw(1:3, 1:3), Tcw(1:3, 4));
        [bearings, ~] = NormalizeVector((inv(K) * pextend(ptIcs'))');
        VM{k,1}{i,1} = bearings;
    end
end
% 这里一定要固定住2个pose，这样尺度才会被固定住，就像vio里可以用纯极线约束的原因是，Tbc是固定住的，所以不担心尺度有问题。
fix_pose_num = 2;

pose_num = length(T_bc) - 2;
pose_offset = 6;
iter_num = 10;

T_bc_gt = T_bc;
for i = fix_pose_num+1 : length(T_bc)
    T_bc{i} = [rodrigues(0.01*(rand(3,1)-0.5)) 100*(rand(3,1)-0.5); 0 0 0 1] * T_bc{i};
    
end

for k = 1 : iter_num
    
    H = zeros(pose_offset * pose_num, pose_offset * pose_num);
    b = zeros(pose_offset * pose_num,1);
    energy = 0;
    for i = 1 : length(VM) % loop timestamps
        for cam_i = 1 : length(VM{i})
            vm_is = VM{i}{cam_i};
            Tbci = T_bc{cam_i};
            for cam_j = 1 : length(VM{i})
                if cam_i >= cam_j
                    continue;
                end
                if cam_i <= 2 || cam_j <= 2
%                     continue;
                end
                vm_js = VM{i}{cam_j};
                Tbcj = T_bc{cam_j};
                assert(size(vm_is,1) == size(vm_js,1));
                for id = 1 : size(vm_is,1)
                    vm_i = vm_is(id,:)';
                    vm_j = vm_js(id,:)';
                    [res, d_err_d_T_cami, d_err_d_T_camj] = computeEpipolarFactor(vm_i, vm_j, Tbci, Tbcj);
                    energy = energy + abs(res);
                    if cam_i == 1
                        d_err_d_T_cami = zeros(size(d_err_d_T_cami));
                    end
                    if cam_j == 1
                        d_err_d_T_camj = zeros(size(d_err_d_T_camj));
                    end
                    try
                        H(pose_offset * (cam_i-fix_pose_num-1) + 1 : pose_offset * (cam_i-fix_pose_num), pose_offset * (cam_i-fix_pose_num-1) + 1 : pose_offset * (cam_i-fix_pose_num)) = H(pose_offset * (cam_i-fix_pose_num-1) + 1 : pose_offset * (cam_i-fix_pose_num), pose_offset * (cam_i-fix_pose_num-1) + 1 : pose_offset * (cam_i-fix_pose_num)) + d_err_d_T_cami' * d_err_d_T_cami;
                    catch
                    end
                    try
                        H(pose_offset * (cam_j-fix_pose_num-1) + 1 : pose_offset * (cam_j-fix_pose_num), pose_offset * (cam_j-fix_pose_num-1) + 1 : pose_offset * (cam_j-fix_pose_num)) = H(pose_offset * (cam_j-fix_pose_num-1) + 1 : pose_offset * (cam_j-fix_pose_num), pose_offset * (cam_j-fix_pose_num-1) + 1 : pose_offset * (cam_j-fix_pose_num)) + d_err_d_T_camj' * d_err_d_T_camj;
                    catch
                    end
                    try
                        H(pose_offset * (cam_i-fix_pose_num-1) + 1 : pose_offset * (cam_i-fix_pose_num), pose_offset * (cam_j-fix_pose_num-1) + 1 : pose_offset * (cam_j-fix_pose_num)) = H(pose_offset * (cam_i-fix_pose_num-1) + 1 : pose_offset * (cam_i-fix_pose_num), pose_offset * (cam_j-fix_pose_num-1) + 1 : pose_offset * (cam_j-fix_pose_num)) + d_err_d_T_cami' * d_err_d_T_camj;
                    catch
                    end
                    try
                        H(pose_offset * (cam_j-fix_pose_num-1) + 1 : pose_offset * (cam_j-fix_pose_num), pose_offset * (cam_i-fix_pose_num-1) + 1 : pose_offset * (cam_i-fix_pose_num)) = H(pose_offset * (cam_j-fix_pose_num-1) + 1 : pose_offset * (cam_j-fix_pose_num), pose_offset * (cam_i-fix_pose_num-1) + 1 : pose_offset * (cam_i-fix_pose_num)) + d_err_d_T_camj' * d_err_d_T_cami;
                    catch
                    end
                    try
                        b(pose_offset * (cam_j-fix_pose_num-1) + 1 : pose_offset * (cam_j-fix_pose_num), 1) = b(pose_offset * (cam_j-fix_pose_num-1) + 1 : pose_offset * (cam_j-fix_pose_num), 1) + d_err_d_T_camj' * res;
                    catch
                    end
                    try
                        b(pose_offset * (cam_i-fix_pose_num-1) + 1 : pose_offset * (cam_i-fix_pose_num), 1) = b(pose_offset * (cam_i-fix_pose_num-1) + 1 : pose_offset * (cam_i-fix_pose_num), 1) + d_err_d_T_cami' * res;
                    catch
                    end
                    
                end
                
            end
        end
    end
    fprintf(sprintf('iter: %d, fx: %f\n', k, energy));
    inc = -inv(H) * b;
    update = reshape(inc, 6, []);
    for kk = fix_pose_num+1 : length(T_bc)
        T_bc{kk} = [rodrigues(update(1:3,kk-fix_pose_num)) update(4:6,kk-fix_pose_num); 0 0 0 1] * T_bc{kk};
    end
end


end
function [res, d_err_d_T_cami, d_err_d_T_camj] = computeEpipolarFactor(vm_i, vm_j, Tbci, Tbcj)
R_c0_ci = Tbci(1:3,1:3);
t_c0_ci = Tbci(1:3,4);
R_c0_cj = Tbcj(1:3,1:3);
t_c0_cj = Tbcj(1:3,4);

plane = SkewSymMat(R_c0_ci' * (t_c0_cj - t_c0_ci)) * R_c0_ci' * R_c0_cj * vm_j;


plane_normalized = plane./norm(plane);
res = dot(vm_i, plane_normalized);
if 0
    d_plane_d_R_cami = -R_c0_ci' * SkewSymMat(R_c0_cj * vm_j) * SkewSymMat(t_c0_ci) - R_c0_ci' * SkewSymMat(SkewSymMat(R_c0_cj * vm_j) * (t_c0_cj - t_c0_ci));
    d_plane_d_R_camj = -R_c0_ci' * SkewSymMat(SkewSymMat(t_c0_cj) * R_c0_cj * vm_j) + R_c0_ci' * SkewSymMat(t_c0_ci) * SkewSymMat(R_c0_cj * vm_j);
else
    d_plane_d_R_cami = SkewSymMat(R_c0_ci' * (t_c0_cj - t_c0_ci)) * R_c0_ci' * SkewSymMat(R_c0_cj * vm_j) - SkewSymMat(R_c0_ci' * R_c0_cj * vm_j) * R_c0_ci' * SkewSymMat(t_c0_cj);
    d_plane_d_R_camj = -d_plane_d_R_cami;
end
check = abs(d_plane_d_R_cami + d_plane_d_R_camj);
try
    assert(max(check(:)) < 1e-13);
catch
    sdfhk = 1;
end
d_plane_d_t_cami = SkewSymMat(R_c0_ci' * R_c0_cj * vm_j) * R_c0_ci';
d_plane_d_t_camj = -d_plane_d_t_cami;

J_err_plane = zeros(1,3);

J_err_plane(1) = vm_i(1)/norm(plane) - res * plane(1) / norm(plane)^3;
J_err_plane(2) = vm_i(2)/norm(plane) - res * plane(2) / norm(plane)^3;
J_err_plane(3) = vm_i(3)/norm(plane) - res * plane(3) / norm(plane)^3;

d_err_d_T_cami = J_err_plane * [d_plane_d_R_cami d_plane_d_t_cami];
d_err_d_T_camj = J_err_plane * [d_plane_d_R_camj d_plane_d_t_camj];

end