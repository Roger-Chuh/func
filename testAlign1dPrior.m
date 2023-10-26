function testAlign1dPrior()
global use_dist_error use_5dof pose_size opt_3dof add_noise right_perturbation include_first_pose offset fixed_pose_num use_reproj_factor use_epiplane_factor...
    fix_trans_norm use_xyz fix_rot
add_noise = 0;
use_5dof = 0;
fix_rot = 0;
right_perturbation = 0;

pose_num = 1;
cam_num = 4;

%% T_wc
pose{1,1} = [rodrigues([0 0 0]) [0 0 0]';0 0 0 1]; % host frame c0
pose{2,1} = [rodrigues([0.1 0.02 0.03]) [10 20 30]';0 0 0 1]; % target frame c0
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
[host_bearings, dist] = NormalizeVector(XYZ);
lambdas = 1 ./ dist;


K = [400 0 320; 0 400 240;0 0 1];

poses = {};

for j = 1 : size(XYZ, 1)
    target_vm_count = 1;
    target_bearing{j,1} = [];
    for i = 1 : length(pose)
        [ptIcs(i,:), tgtPt3d] = TransformAndProject(XYZ(j,:), K, pose{i,1}(1:3, 1:3)', -pose{i,1}(1:3, 1:3)'*pose{i,1}(1:3, 4));
        bearing = inv(K) * [ptIcs(i,:) 1]';
        bearings{j,1}(i,:) = bearing'./norm(bearing);
        if i == 1
            host_bearing(j,:) = bearings{j,1}(i,:);
        else
            target_bearing{j,1}(target_vm_count,:) = bearings{j,1}(i,:);
            target_vm_count = target_vm_count + 1;
        end
    end
end
% pt3d = triangulate(poses, bearings);

iter_num = 10;


lambdas_use = lambdas + 0.001;

T_target_use = T_target * [rodrigues([1.1 2.2 -0.3]) [0.1 10 -20]'; 0 0 0 1];
pose_offset = 6;
point_num = size(XYZ,1);
for k = 1 : iter_num
    
    H = zeros(pose_offset + point_num, pose_offset + point_num);
    b = zeros(pose_offset + point_num,1);
    fx = 0;
    
    Twb_host = T_host;
    Twb_target = T_target_use;
    Twc0_host = Twb_host * T_bc{1,1};
    
    
    [err2, d_err_d_target_pose2] = computeRelativePoseFactor(T_target_use, T_target);
    fx = fx + norm(err2);
    H(1:pose_offset,1:pose_offset) = H(1:pose_offset,1:pose_offset) + d_err_d_target_pose2' * d_err_d_target_pose2;
    b(1:pose_offset,1) = b(1:pose_offset,1) + (-d_err_d_target_pose2' * err2);
    
    for ii = 1 : point_num
        lambda_use = lambdas_use(ii);
        point_id = ii;
        
        for i = 1 : pose_num
            %             Twb_host = T_host;
            %             Twb_target = T_target_use;
            %             Twc0_host = Twb_host * T_bc{1,1};
            
            %             [err2, d_err_d_target_pose2] = computeRelativePoseFactor(T_target_use, T_target);
            %             fx = fx + norm(err2);
            %             H(1:pose_offset,1:pose_offset) = H(1:pose_offset,1:pose_offset) + d_err_d_target_pose2' * d_err_d_target_pose2;
            %             b(1:pose_offset,1) = b(1:pose_offset,1) + (-d_err_d_target_pose2' * err2);
            
            for j = 1 : cam_num
                Twci_target = Twb_target * T_bc{j,1};
                [err_, d_err_d_host_pose, d_err_d_target_pose, JacPoint_true] = computeReprojFactor(Twc0_host, Twci_target, lambda_use, host_bearing(ii,:)', target_bearing{ii,1}(j,:)');
                fx = fx + norm(err_);
                H(1:pose_offset,1:pose_offset) = H(1:pose_offset,1:pose_offset) + d_err_d_target_pose' * d_err_d_target_pose;
                H(pose_offset + point_id, pose_offset + point_id) = H(pose_offset + point_id, pose_offset + point_id) + JacPoint_true' * JacPoint_true;
                H(1:pose_offset,pose_offset + point_id) = H(1:pose_offset,pose_offset + point_id) + d_err_d_target_pose' * JacPoint_true;
                H(pose_offset + point_id,1:pose_offset) = H(pose_offset + point_id,1:pose_offset) + JacPoint_true' * d_err_d_target_pose;
                
                b(1:pose_offset,1) = b(1:pose_offset,1) + (-d_err_d_target_pose' * err_);
                b(pose_offset + point_id,1) = b(pose_offset + point_id,1) + (-JacPoint_true' * err_);
                
                %                 [err2, d_err_d_target_pose2] = computeRelativePoseFactor(T_target_use, T_target);
                %                 fx = fx + norm(err2);
                %                 H(1:pose_offset,1:pose_offset) = H(1:pose_offset,1:pose_offset) + d_err_d_target_pose2' * d_err_d_target_pose2;
                %                 b(1:pose_offset,1) = b(1:pose_offset,1) + (-d_err_d_target_pose2' * err2);
            end
        end
    end
    fprintf(sprintf('iter: %d, fx: %f\n', k, fx));
%         fx
    inc = inv(H) * b;
    T_target_use = [rodrigues(inc(1:3)) inc(4:6);0 0 0 1] * T_target_use;
    lambdas_use = lambdas_use + inc(pose_offset+1:end);
    
end

end

function point_3d = triangulate(poses, points)

design_matrix = zeros(length(poses)*2, 4);
for (i = 1 : length(poses))
    p0x = points(i,1);
    p0y = points(i,2);
    p0z = points(i,3);
    design_matrix(i*2-1,:) = p0x * poses{i,1}(3,:) - p0z*poses{i,1}(1,:);
    design_matrix(i*2,:) = p0y * poses{i,1}(3,:) - p0z*poses{i,1}(2,:);
end


[A, B, C] = svd(design_matrix,0);
[x, y, z] = svd(design_matrix' * design_matrix,0);
triangulated_point = C (:, end);
point_3d = triangulated_point(1:3) ./ triangulated_point(4);


errs = design_matrix*[point_3d; 1];
error =  norm(errs)/ length(errs);

end
function [A] = ProduceOtherOthogonalBasis(n)
N = n;
if (N(0+1) < 0)
    N(0+1) = -N(0+1);
end
if (N(1+1) < 0)
    N(1+1) = -N(1+1);
end
if (N(2+1) < 0)
    N(2+1) = -N(2+1);
end
minIdx = 0+1;
if (N(0+1) <= N(1+1))
    if (N(0+1) <= N(2+1))
        minIdx = 0+1;
    else
        minIdx = 2+1;
    end
else
    if (N(1+1) <= N(2+1))
        minIdx = 1+1;
    else
        minIdx = 2+1;
    end
end
A = zeros(3,2);
switch (minIdx)
    case 1
        A(:,1) = [0, -n(2+1), n(1+1)]';
        
    case 2
        A(:,1) = [n(2+1), 0, -n(0+1)]';
        
    case 3
        A(:,1) = -[-n(1+1), n(0+1), 0]';
end
A(:,2) = cross(n,A(:,1));


% A(:,1) = A(:,1)./norm(A(:,1));
% A(:,2) = A(:,2)./norm(A(:,2));


end
function [err_, d_err_d_host_pose, d_err_d_target_pose, JacPoint_true] = computeReprojFactor(Twc_host, Twc_target, lambda, n,ob)
global use_5dof right_perturbation fix_rot
Twc_h = Twc_host;
Twc_t = Twc_target;
Tcw_t = inv(Twc_t);
Tct_ch_ = Tcw_t * Twc_h;
R_th = Tct_ch_(1:3,1:3);
t_th = Tct_ch_(1:3,4);

project_vec_ = R_th * n + lambda * t_th;
project_vec_ = project_vec_ / lambda;
d = norm(project_vec_);
d2 = d * d;
project_vec_ = project_vec_ / d;
err_ = project_vec_ - ob;

dH_ = eye(3) - project_vec_ * project_vec_';
dH_ = dH_ * 1 / d2;
dH_true = d * dH_;
JacPoint_true = dH_true * (-R_th * n / lambda / lambda);

dH_point_ = JacPoint_true' * JacPoint_true;
db_point_ = JacPoint_true' * err_;

Rwc_h = Twc_h(1:3,1:3);
twc_h = Twc_h(1:3,4);
pro_vec_1_ = Rwc_h * n + twc_h * lambda;
%��������ϵ�µ�3d�㣨ע�ⲻ�Ǳ�����host����ϵ�µģ�
pro_vec_1_ = pro_vec_1_/ lambda;
pro_vec_2_ = R_th * n + R_th * lambda;
%target��ϵ�µ�3d�㣨ע�ⲻ�Ǳ�����host����ϵ�µģ�
pro_vec_2_ = pro_vec_2_/ lambda;
Rcw_t = Tcw_t(1:3,1:3);
%       Mat3RightMultiplySkewM3V3(Rcw_t, pro_vec_1_, Asb_);
Asb_ = Rcw_t * SkewSymMat(pro_vec_1_);
JacPose_(1:3,1:3) = -Asb_;
JacPose_(1:3,4:6) = Rcw_t;
JacPose_true = dH_true * JacPose_;
%       JacPose_ = dH_ * JacPose_;
%       dH_pose_h_ = JacPose_true' * JacPose_true;
%       db_pose_h_ = JacPose_true' * err_;
%       dH_pose_h_point_ = JacPose_true' * JacPoint_true;


if right_perturbation
    d_err_d_R_host = dH_true * (-R_th * SkewSymMat(n./lambda));
    d_err_d_R_target = dH_true * SkewSymMat(pro_vec_2_);
    
    d_err_d_t_host = dH_true * R_th;
    d_err_d_t_target = -dH_true;
    
end


if ~use_5dof
    %     assert(right_perturbation == 0)
    if ~right_perturbation
        d_err_d_host_pose = JacPose_true;
        d_err_d_target_pose = -JacPose_true;
    else
        d_err_d_host_pose = [d_err_d_R_host d_err_d_t_host];
        d_err_d_target_pose = [d_err_d_R_target d_err_d_t_target];
    end
else
    d_err_d_host_pose = zeros(3,5);
    d_err_d_target_pose = zeros(3,5);
    
    twc_i = Twc_host(1:3,4);
    twc_j = Twc_target(1:3,4);
    
    Ag_i= ProduceOtherOthogonalBasis(twc_i);
    Ag_j = ProduceOtherOthogonalBasis(twc_j);
    if ~right_perturbation
        d_err_d_host_pose(:,1:3) = JacPose_true(:,1:3);
        d_err_d_target_pose(:,1:3) = -JacPose_true(:,1:3);
    else
        d_err_d_host_pose(:,1:3) = d_err_d_R_host;
        d_err_d_target_pose(:,1:3) = d_err_d_R_target;
    end
    
    d_err_d_host_pose(:,4:5) = dH_true * (-Twc_target(1:3,1:3)' * SkewSymMat(Twc_host(1:3,4))*Ag_i);
    d_err_d_target_pose(:,4:5) = dH_true * (Twc_target(1:3,1:3)' * SkewSymMat(Twc_target(1:3,4))*Ag_j);
    
    if (fix_rot)
        d_err_d_host_pose(:,1:2) = dH_true * (-Twc_target(1:3,1:3)' * SkewSymMat(Twc_host(1:3,4))*Ag_i);
        d_err_d_target_pose(:,1:2) = dH_true * (Twc_target(1:3,1:3)' * SkewSymMat(Twc_target(1:3,4))*Ag_j);
        
    end
    
end


if ~right_perturbation && use_5dof
    
    %     d_err_d_R_host = dH_true * Twc_target(1:3,1:3)' * SkewSymMat(Twc_host(1:3,1:3) * n ./ lambda);
    %     d_err_d_R_target = dH_true * Twc_target(1:3,1:3)' * SkewSymMat(Twc_host(1:3,1:3) * n ./ lambda + Twc_host(1:3,4) - Twc_target(1:3,4));
    %
    %     d_err_d_host_pose = zeros(3,5);
    %     d_err_d_target_pose = zeros(3,5);
    %
    %     d_err_d_host_pose(:,1:3) = d_err_d_R_host;
    %     d_err_d_target_pose(:,1:3) = d_err_d_R_target;
    %
    %     d_err_d_host_pose(:,4:5) = -Twc_target(1:3,1:3)' * SkewSymMat(Twc_host(1:3,4))*Ag_i;
    %     d_err_d_target_pose(:,4:5) = Twc_target(1:3,1:3)' * SkewSymMat(Twc_target(1:3,4))*Ag_j;
end



end
function [err, d_err_d_target_pose] = computeRelativePoseFactor(T, T_prior)

deltaPose = inv(T) * T_prior;

err = [rodrigues(deltaPose(1:3,1:3)); deltaPose(1:3,4)];

d_err_d_target_pose = zeros(6,6);
d_err_d_target_pose(1:3,1:3) = -JrInv(err(1:3)) * T_prior(1:3,1:3)';  % d_R_err_d_R
d_err_d_target_pose(4:6,1:3) = T(1:3,1:3)' * SkewSymMat(T_prior(1:3,4));  % d_t_err_d_R
d_err_d_target_pose(4:6,4:6) = -T(1:3,1:3)'; % d_t_err_d_t

end
function J =  JrInv(phi)
EPSILON = 1e-16;
EPSILONSQRT = sqrt(EPSILON);
kOur_PI = 3.141592653;
J = eye(3);

phi_norm2 = norm(phi)^2;
phi_hat = SkewSymMat(phi);
phi_hat2 = phi_hat * phi_hat;

J = J + phi_hat / 2;
if (phi_norm2 > EPSILON)
    phi_norm = sqrt(phi_norm2);
    
    %     assert(phi_norm <= kOur_PI + EPSILON && "We require that the angle is in range [0, pi].");
    
    if (phi_norm < kOur_PI - EPSILONSQRT)
        
        J = J + phi_hat2 * (1 / phi_norm2 - (1 + cos(phi_norm)) / (2 * phi_norm * sin(phi_norm)));
    else
        J = J + phi_hat2 / (kOur_PI * kOur_PI);
    end
else
    J = J + phi_hat2 / 12;
end

end