function testEpiplaneFactor()
global use_dist_error use_5dof pose_size opt_3dof add_noise right_perturbation include_first_pose offset fixed_pose_num use_reproj_factor use_epiplane_factor...
        fix_trans_norm
add_noise = 1;
opt_3dof = 0;
use_dist_error = 0;1;
use_5dof = 1;
right_perturbation = 0;
use_reproj_factor = 0;
use_epiplane_factor = 1;
pose_num = 5;


if ~use_5dof
%     right_perturbation = 0;
end

fixed_pose_num = 1;
include_first_pose = 0;1;

if include_first_pose
    offset = 1;
else
    offset = 0;
end

if use_5dof
    pose_size = 5;
else
    pose_size = 6;
end

if(use_5dof)
   fix_trans_norm = true; 
end


if opt_3dof
    pose_size = 3;
end

n = rand(3,1);
n = n./norm(n);
pt = 0.5*n;
n1 = rand(3,1);
n1 = n1./norm(n1);
n2 = cross(n,n1);
dir = n1;

pt=dir-dot(dir,n) * n;
a = (dir-pt)./norm(dir-pt);
[a-n]
dot(pt, n)



K = [300 0 320; 0 300 240; 0 0 1];

xyz = [10 200 1000;200 -20 1002;-10 500 1003;-100 -50 990;300 -50 990;-100 -250 990;-100 -150 990];

xyz = [xyz; xyz + [-10 20 3]; xyz + [10 50 3]; xyz + [-70 -20 3]; xyz + [-10 20 30]];

[~,depths] = NormalizeVector(xyz);

idps = 1./depths;

pose{1,1} = eye(4);
for i = 1 : pose_num %size(xyz,1)
    if 0 %i == 1
        pose{i,1} = eye(4);
    else
        pose{i,1} = [rodrigues(0.2*(rand(3,1)-0.5)) 0.5*(rand(3,1)-0.5); 0 0 0 1];
    end
    
    bearing = pose{i,1}(1:3,1:3)*xyz' + repmat(pose{i,1}(1:3,4), 1, size(xyz,1));
    [bearing_] = NormalizeVector(bearing');
    bearings{i,1} = bearing_;
    
end


opt_fids = fixed_pose_num+1:pose_num;%size(xyz,1);

pose_Twc = pose;
for k = 1 : fixed_pose_num
    pose_err{k,1} = eye(4);
    pose_Twc{k,1} = inv(pose{k,1});
end
for i = opt_fids
    pose_err{i,1} = [rodrigues(0.1*(rand(3,1)-0.5)) 0.1*(rand(3,1)-0.5); 0 0 0 1];
    if add_noise
        %        pose_Twc{i,1} = pose_err{i,1} * inv(pose{i,1});
        pose_Twc{i,1} = inv(pose{i,1}) * pose_err{i,1};
        Twc_gt_ = inv(pose{i,1});
        pose_Twc{i,1}(1:3,4) = pose_Twc{i,1}(1:3,4)./norm(pose_Twc{i,1}(1:3,4)).*norm(pose{i,1}(1:3,4));
        %        pose_err{i,1} = pose_Twc{i,1} * pose{i,1};
        pose_err{i,1} = pose{i,1} * pose_Twc{i,1};
        if opt_3dof
            pose_Twc{i,1}(1:3,4) = Twc_gt_(1:3,4);
            %            pose_err{i,1} = pose_Twc{i,1} * pose{i,1};
            pose_err{i,1} = pose{i,1} * pose_Twc{i,1};
        end
    else
        pose_Twc{i,1} = inv(pose{i,1});
    end
end


host_fid = 1;

[xMat, yMat] = meshgrid(1 : pose_num, 1 :pose_num);
pairs = [xMat(:) yMat(:)];

loss_vec = [];


for iter = 1 : 30
    
    H = zeros(pose_size*(pose_num-fixed_pose_num+1-1+offset), pose_size*(pose_num-fixed_pose_num+1-1+offset));
    b = zeros(pose_size*(pose_num-fixed_pose_num+1-1+offset), 1);
    err_sum = 0;
    err_count = 0;
    for pid = 1 : size(xyz,1)
        pt = xyz(pid,:);
        host_bearing = bearings{1}(pid,:);
%         assert(norm(pt./norm(pt) - host_bearing) < 0.000001);
        for pair_id = 1 : size(pairs, 1)
            host_fid = pairs(pair_id, 1);
            target_fid = pairs(pair_id, 2);
            if host_fid >= target_fid
                continue;
            end
            if host_fid > fixed_pose_num
                continue;
            end
            if target_fid <= fixed_pose_num
                continue;
            end
%             [host_fid target_fid]
            host_bearing = bearings{host_fid}(pid,:);
            Twc_host = pose_Twc{host_fid,1};
            Twc_host_gt = inv(pose{host_fid,1});
            
            target_bearing = bearings{target_fid}(pid,:);
            reproj =  pose{target_fid,1}(1:3,1:3)*pt' + pose{target_fid,1}(1:3,4);
            assert(norm(target_bearing' - reproj./norm(reproj)) < 0.000001);
            
            pt_at_host = pose{host_fid,1}(1:3,1:3)*pt' + pose{host_fid,1}(1:3,4);
            
            Twc_target = pose_Twc{target_fid,1};
            Twc_target_gt = inv(pose{target_fid,1});
            [err, d_err_d_host_pose, d_err_d_target_pose] = computeEpiplaneFacor(Twc_host, Twc_target, host_bearing, target_bearing);
            if use_epiplane_factor
                err_sum = err_sum + 235*norm(err);
                err_count = err_count+1;
            end
            if(host_fid <= fixed_pose_num)
                d_err_d_host_pose = zeros(size(d_err_d_host_pose));
            end
            if(target_fid <= fixed_pose_num)
                d_err_d_target_pose = zeros(size(d_err_d_target_pose));
            end
%             delta_target_pose = pose_err{target_fid,1};
%             delta_target_pose_vec = [rodrigues(delta_target_pose(1:3,1:3));delta_target_pose(1:3,4)];
%             delta_host_pose = pose_err{host_fid,1};
%             delta_host_pose_vec = [rodrigues(delta_host_pose(1:3,1:3));delta_host_pose(1:3,4)];
%             if ~use_5dof
%                 err_comp = d_err_d_target_pose * delta_target_pose_vec + d_err_d_host_pose * delta_host_pose_vec;
%             else
%                 
%             end
            
            d_err_d_host_pose = d_err_d_host_pose(:,1:pose_size);
            d_err_d_target_pose = d_err_d_target_pose(:,1:pose_size);
            %       [err err_comp (err - err_comp)]
            
            if use_epiplane_factor
                H = FillMatrix(H, (host_fid-2+offset-fixed_pose_num+1)*pose_size+1, (host_fid-2+offset-fixed_pose_num+1)*pose_size+1, pose_size, pose_size, d_err_d_host_pose' * d_err_d_host_pose);
                H = FillMatrix(H, (target_fid-2+offset-fixed_pose_num+1)*pose_size+1, (target_fid-2+offset-fixed_pose_num+1)*pose_size+1, pose_size, pose_size, d_err_d_target_pose' * d_err_d_target_pose);
                H = FillMatrix(H, (host_fid-2+offset-fixed_pose_num+1)*pose_size+1, (target_fid-2+offset-fixed_pose_num+1)*pose_size+1, pose_size, pose_size, d_err_d_host_pose' * d_err_d_target_pose);
                H = FillMatrix(H, (target_fid-2+offset-fixed_pose_num+1)*pose_size+1, (host_fid-2+offset-fixed_pose_num+1)*pose_size+1, pose_size, pose_size, d_err_d_target_pose' * d_err_d_host_pose);
                
                b = FillMatrix(b, (host_fid-2+offset-fixed_pose_num+1)*pose_size+1, 1, pose_size, 1, -d_err_d_host_pose' * err);
                b = FillMatrix(b, (target_fid-2+offset-fixed_pose_num+1)*pose_size+1, 1, pose_size, 1, -d_err_d_target_pose' * err);
            end
            %       d_err_d_host_pose *
            
            
            
            if use_reproj_factor
                [err1, d_err_d_host_pose1, d_err_d_target_pose1] = computeReprojFactor(Twc_host, Twc_target, 1./norm(pt_at_host), host_bearing', target_bearing');
                err_sum = err_sum + 235*norm(err1);
                err_count = err_count+1;
                if(host_fid <= fixed_pose_num)
                    d_err_d_host_pose1 = zeros(size(d_err_d_host_pose1));
                end
                if(target_fid <= fixed_pose_num)
                    d_err_d_target_pose1 = zeros(size(d_err_d_target_pose1));
                end
                d_err_d_host_pose1 = d_err_d_host_pose1(:,1:pose_size);
                d_err_d_target_pose1 = d_err_d_target_pose1(:,1:pose_size);
                H = FillMatrix(H, (host_fid-2+offset-fixed_pose_num+1)*pose_size+1, (host_fid-2+offset-fixed_pose_num+1)*pose_size+1, pose_size, pose_size, d_err_d_host_pose1' * d_err_d_host_pose1);
                H = FillMatrix(H, (target_fid-2+offset-fixed_pose_num+1)*pose_size+1, (target_fid-2+offset-fixed_pose_num+1)*pose_size+1, pose_size, pose_size, d_err_d_target_pose1' * d_err_d_target_pose1);
                H = FillMatrix(H, (host_fid-2+offset-fixed_pose_num+1)*pose_size+1, (target_fid-2+offset-fixed_pose_num+1)*pose_size+1, pose_size, pose_size, d_err_d_host_pose1' * d_err_d_target_pose1);
                H = FillMatrix(H, (target_fid-2+offset-fixed_pose_num+1)*pose_size+1, (host_fid-2+offset-fixed_pose_num+1)*pose_size+1, pose_size, pose_size, d_err_d_target_pose1' * d_err_d_host_pose1);
                
                b = FillMatrix(b, (host_fid-2+offset-fixed_pose_num+1)*pose_size+1, 1, pose_size, 1, -d_err_d_host_pose1' * err1);
                b = FillMatrix(b, (target_fid-2+offset-fixed_pose_num+1)*pose_size+1, 1, pose_size, 1, -d_err_d_target_pose1' * err1);
            end
            
        end
        
        
        
    end
    
    if(include_first_pose)
        
        H(1:pose_size, 1:pose_size) = H(1:pose_size, 1:pose_size) + 1e10.*eye(pose_size);
    end
    dx = inv(H) * b;
    fprintf(sprintf('iter: %d, err_sum: %.20f, err_count: %d\n', iter, err_sum, err_count));
    loss_vec = [loss_vec; err_sum];
    err_vec = [];
    for ind = 1 : pose_num-1+offset-fixed_pose_num+1
        update_vec = dx((ind-1)*pose_size+1:ind*pose_size);
        if ~use_5dof
%             assert(right_perturbation == 0);
            if(~opt_3dof)
                if ~right_perturbation
                    pose_Twc{ind+1-offset+fixed_pose_num-1} = [rodrigues(update_vec(1:3)) update_vec(4:pose_size);0 0 0 1] * pose_Twc{ind+1-offset+fixed_pose_num-1};
                else
                    pose_Twc{ind+1-offset+fixed_pose_num-1} = pose_Twc{ind+1-offset+fixed_pose_num-1} * [rodrigues(update_vec(1:3)) update_vec(4:pose_size);0 0 0 1];
                end
            else
                pose_Twc{ind+1-offset+fixed_pose_num-1}(1:3,1:3) = [rodrigues(update_vec(1:3))] * pose_Twc{ind+1-offset+fixed_pose_num-1}(1:3,1:3);
            end
        else
            Twc_old = pose_Twc{ind+1-offset+fixed_pose_num-1};
            if ~opt_3dof
                delta_r = ProduceOtherOthogonalBasis(Twc_old(1:3,4)) * update_vec(4:pose_size);
                trans_new_2dof = rodrigues(delta_r) * Twc_old(1:3,4);
                pose_Twc{ind+1-offset+fixed_pose_num-1}(1:3,4) = trans_new_2dof;
            end
            if ~right_perturbation
                pose_Twc{ind+1-offset+fixed_pose_num-1}(1:3,1:3) = rodrigues(update_vec(1:3)) * Twc_old(1:3,1:3);
            else
                pose_Twc{ind+1-offset+fixed_pose_num-1}(1:3,1:3) = Twc_old(1:3,1:3) * rodrigues(update_vec(1:3));
            end
            
        end
        pose_gt = inv(pose{ind+1-offset+fixed_pose_num-1,1});
        pose_gt_est = pose_Twc{ind+1-offset+fixed_pose_num-1};
        delta_pose = inv(pose_gt_est) * pose_gt;
        rot_err = rad2deg(norm(rodrigues(pose_gt(1:3,1:3)' * pose_Twc{ind+1-offset+fixed_pose_num-1}(1:3,1:3))));
        norm_a = pose_gt(1:3,4)./norm(pose_gt(1:3,4));
        norm_b = pose_Twc{ind+1-offset+fixed_pose_num-1}(1:3,4)./norm(pose_Twc{ind+1-offset+fixed_pose_num-1}(1:3,4));
        trans_norm_err = norm(norm_a.*sign(norm_a(3)) - norm_b.*sign(norm_b(3)));
        err_vec = [err_vec; [rot_err trans_norm_err norm(delta_pose(1:3,4))]];
    end
    err_vec
end


figure(100),plot(loss_vec);title('loss_vec');

end
function H = FillMatrix(H, start_row, start_col, row_size, col_size, value)

if (start_row < 1 || start_col < 1)
    return
end

H(start_row:start_row+row_size-1, start_col:start_col+col_size-1) = H(start_row:start_row+row_size-1, start_col:start_col+col_size-1) + value;
end
function [err, d_err_d_host_pose, d_err_d_target_pose] = computeEpiplaneFacor(Twc_host, Twc_target, host_bearing, target_bearing)
global use_dist_error use_5dof pose_size
T_th = inv(Twc_target) * Twc_host;
plane = SkewSymMat(T_th(1:3,4)) * T_th(1:3,1:3) * host_bearing';
plane_norm = plane./norm(plane);

pt_in_plane = target_bearing' - (target_bearing * plane_norm) .*plane_norm;
plane_norm_check = (pt_in_plane - target_bearing')./norm(pt_in_plane - target_bearing');

if 0
    min([norm(plane_norm_check + plane_norm) norm(plane_norm_check - plane_norm)])
    dot(pt_in_plane,plane_norm)
end
if ~use_dist_error
    err = (pt_in_plane./norm(pt_in_plane) - target_bearing');
else
    res = dot(target_bearing', plane);
    err = dot(target_bearing', plane_norm);
end




d_err_d_pt = compute_d_bearing_d_pt_jac(pt_in_plane);

d_pt_d_plane_norm = compute_d_pt_d_plane_norm_jac(target_bearing', plane_norm);

d_plane_norm_d_plane = compute_d_bearing_d_pt_jac(plane);


[d_plane_d_R_host, d_plane_d_t_host, d_plane_d_R_target, d_plane_d_t_target] = compute_d_plane_d_pose(Twc_host, Twc_target, host_bearing');


if ~use_dist_error
    d_err_d_plane = d_err_d_pt * d_pt_d_plane_norm * d_plane_norm_d_plane;
    
    d_err_d_host_pose = d_err_d_plane * [d_plane_d_R_host d_plane_d_t_host];
    
    d_err_d_target_pose = d_err_d_plane * [d_plane_d_R_target d_plane_d_t_target];
else
    J_err_plane = compute_d_dist_d_plane(res,  target_bearing, plane);
    d_err_d_host_pose = J_err_plane * [d_plane_d_R_host d_plane_d_t_host];
    d_err_d_target_pose = J_err_plane * [d_plane_d_R_target d_plane_d_t_target];
end
end
function J_err_plane = compute_d_dist_d_plane(res, p3d_target, plane)
plane_norm = norm(plane);
J_err_plane(0+1) = p3d_target(0+1)/plane_norm - res * plane(0+1) / plane_norm^3;
J_err_plane(1+1) = p3d_target(1+1)/plane_norm - res * plane(1+1) / plane_norm^3;
J_err_plane(2+1) = p3d_target(2+1)/plane_norm - res * plane(2+1) / plane_norm^3;
end
function d_bearing_d_pt = compute_d_bearing_d_pt_jac(pt)

d_bearing_d_pt = (norm(pt).*eye(3) - pt * pt'./norm(pt))./(norm(pt)^2);

if 0
    
    ptNorm = norm(pt);
    ptNorm_3_2 = -0.5 * (1 / (ptNorm * ptNorm * ptNorm));
    d_err_d_pt3d = zeros(3,3);
    d_err_d_pt3d(0+1, 0+1) = 1.0 / ptNorm + pt(0+1) * (ptNorm_3_2 * 2 * pt(0+1));
    d_err_d_pt3d(0+1, 1+1) = pt(0+1) * (ptNorm_3_2 * 2 * pt(1+1));
    d_err_d_pt3d(0+1, 2+1) = pt(0+1) * (ptNorm_3_2 * 2 * pt(2+1));
    
    d_err_d_pt3d(1+1, 0+1) = pt(1+1) * (ptNorm_3_2 * 2 * pt(0+1));
    d_err_d_pt3d(1+1, 1+1) = 1.0 / ptNorm + pt(1+1) * (ptNorm_3_2 * 2 * pt(1+1));
    d_err_d_pt3d(1+1, 2+1) = pt(1+1) * (ptNorm_3_2 * 2 * pt(2+1));
    
    d_err_d_pt3d(2+1, 0+1) = pt(2+1) * (ptNorm_3_2 * 2 * pt(0+1));
    d_err_d_pt3d(2+1, 1+1) = pt(2+1) * (ptNorm_3_2 * 2 * pt(1+1));
    d_err_d_pt3d(2+1, 2+1) = 1.0 / ptNorm + pt(2+1) * (ptNorm_3_2 * 2 * pt(2+1));
end

end
function d_bearing_d_pt = compute_d_pt_d_plane_norm_jac(target_bearing, plane_norm)

d_bearing_d_pt = zeros(3,3);

d_bearing_d_pt(1,1) = -2*target_bearing(1)*plane_norm(1) - target_bearing(2)*plane_norm(2) - target_bearing(3)*plane_norm(3);
d_bearing_d_pt(1,2) = -target_bearing(2)*plane_norm(1);
d_bearing_d_pt(1,3) = -target_bearing(3)*plane_norm(1);

d_bearing_d_pt(2,1) = -target_bearing(1)*plane_norm(2);
d_bearing_d_pt(2,2) = -target_bearing(1)*plane_norm(1) - 2*target_bearing(2)*plane_norm(2) - target_bearing(3)*plane_norm(3);
d_bearing_d_pt(2,3) = -target_bearing(3)*plane_norm(2);

d_bearing_d_pt(3,1) = -target_bearing(1)*plane_norm(3);
d_bearing_d_pt(3,2) = -target_bearing(2)*plane_norm(3);
d_bearing_d_pt(3,3) = -target_bearing(1)*plane_norm(1) - target_bearing(2)*plane_norm(2) - 2*target_bearing(3)*plane_norm(3);

end
function [d_plane_d_R_host, d_plane_d_t_host, d_plane_d_R_target, d_plane_d_t_target] = compute_d_plane_d_pose(Twc_host, Twc_target, host_bearing)
global use_5dof opt_3dof right_perturbation
Rwc_j = Twc_target(1:3,1:3);
Rwc_i = Twc_host(1:3,1:3);
twc_j = Twc_target(1:3,4);
twc_i = Twc_host(1:3,4);
R_cj_ci = Twc_target(1:3,1:3)' * Twc_host(1:3,1:3)';


Ag_i= ProduceOtherOthogonalBasis(twc_i);
Ag_j = ProduceOtherOthogonalBasis(twc_j);

if ~right_perturbation
    d_plane_d_R_host = SkewSymMat(R_cj_ci * host_bearing) * Rwc_j' * SkewSymMat(twc_i) - SkewSymMat(Rwc_j' * (twc_i - twc_j)) * Rwc_j' * SkewSymMat(Rwc_i * host_bearing);
    d_plane_d_R_target = -d_plane_d_R_host;
else
    d_plane_d_R_host = -SkewSymMat(Rwc_j' * (twc_i - twc_j)) * R_cj_ci * SkewSymMat(host_bearing);
    d_plane_d_R_target = SkewSymMat(Rwc_j' * (twc_i - twc_j)) * SkewSymMat(R_cj_ci * host_bearing) - SkewSymMat(R_cj_ci * host_bearing) * SkewSymMat(Rwc_j' * (twc_i - twc_j));
end


if ~use_5dof
    d_plane_d_t_host = -SkewSymMat(R_cj_ci * host_bearing) * Rwc_j';
    d_plane_d_t_target = -d_plane_d_t_host;
else
    d_plane_d_t_host = SkewSymMat(R_cj_ci * host_bearing) * Rwc_j' * SkewSymMat(twc_i) * Ag_i;
    d_plane_d_t_target = -SkewSymMat(R_cj_ci * host_bearing) * Rwc_j' * SkewSymMat(twc_j) * Ag_j;
end

if opt_3dof
    d_plane_d_t_target = zeros(size(d_plane_d_t_target));
    d_plane_d_t_host = zeros(size(d_plane_d_t_target));
end
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
        A(:,1) = [-n(1+1), n(0+1), 0]';
end
A(:,2) = cross(n,A(:,1));


% A(:,1) = A(:,1)./norm(A(:,1));
% A(:,2) = A(:,2)./norm(A(:,2));


end
function [err_, d_err_d_host_pose, d_err_d_target_pose, JacPoint_true] = computeReprojFactor(Twc_host, Twc_target, lambda, n,ob)
global use_5dof right_perturbation
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
pro_vec_1_ = pro_vec_1_/ lambda;
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
    d_err_d_R_target = dH_true * SkewSymMat(pro_vec_1_);
    
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
    
    d_err_d_host_pose(:,4:5) = -Twc_target(1:3,1:3)' * SkewSymMat(Twc_host(1:3,4))*Ag_i;
    d_err_d_target_pose(:,4:5) = Twc_target(1:3,1:3)' * SkewSymMat(Twc_target(1:3,4))*Ag_j;
    
end

end