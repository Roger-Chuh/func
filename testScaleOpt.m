function testScaleOpt()
global use_dist_error use_5dof pose_size opt_3dof add_pose_noise add_point_noise right_perturbation include_first_pose offset fixed_pose_num use_reproj_factor use_epiplane_factor...
    fix_trans_norm use_xyz fix_rot opt_Tbc

close all;

add_pose_noise = true;
add_point_noise = true;
opt_Tbc = true;
if opt_Tbc
  use_5dof = true;
else
  use_5dof = false;
end
fix_Tbc0 = false;
opt_scale = true;
if opt_Tbc
    opt_scale = false;
end
if add_point_noise
%    fix_Tbc0 = true; 
end
% Tbcx = [rodrigues(rand(3, 1)) rand(3, 1); 0 0 0 1];
% Tbcy = [rodrigues(rand(3, 1)) rand(3, 1); 0 0 0 1];
% Twb1 = [rodrigues(rand(3, 1)) rand(3, 1); 0 0 0 1];
% Twb2 = [rodrigues(rand(3, 1)) rand(3, 1); 0 0 0 1];
%
% n = rand(3,1);
% n = n./norm(n);
% n = sign(n(3)).*n;
% lambda = 0.6;
% scale = 1.3;
% Tb2b1 = inv(Twb2) * Twb1;
%
% Tb2b1_scaled = Tb2b1;
% Tb2b1_scaled(1:3,4) = scale * Tb2b1(1:3,4);
%
%
% Tth = inv(Tbcy) * Tb2b1_scaled * Tbcx;
%
%
% X1 = Tth(1:3,1:3) * n + lambda * Tth(1:3,4);
%
% R_th = Tbcy(1:3,1:3)' * Twb2(1:3,1:3)' * Twb1(1:3,1:3) * Tbcx(1:3,1:3);
% t_th = Tbcy(1:3,1:3)' * (Twb2(1:3,1:3)' * (Twb1(1:3,1:3) * Tbcx(1:3,4) + scale * (Twb1(1:3,4) - Twb2(1:3,4))) - Tbcy(1:3,4));
% X2 = R_th * n + lambda * t_th;
%
%
% X1 - X2
K = [300 0 320; 0 300 240; 0 0 1];

xyz = [10 200 1000;200 -20 1002;-10 500 1003;-100 -50 990;300 -50 990;-100 -250 990;-100 -150 990];

xyz = [xyz; xyz + [-10 20 3]; xyz + [10 50 3]; xyz + [-70 -20 3]; xyz + [-10 20 30]; xyz + [-50 15 30]; xyz + [-10 27 -25]; xyz + [75 -20 -30];xyz + [14 -9 31]];

xyz = xyz(1:61,:);

fixed_pose_num = 5;
pose_num = 10;

each_host = ceil(size(xyz,1) / pose_num);

host_ids = [];
target_ids = [];
for(i = 1: size(xyz,1))
    
    host_ids = [host_ids; mod(i, pose_num-1)+1];
    target_ids = [target_ids; mod(i, pose_num-1)+2];
    
end


cam_num = 3;

Tbc = {};
for i = 1 : cam_num
    if 0%i == 1
        Tbc{i, 1} = eye(4);
    else
        Tbc{i, 1} = [rodrigues(0.3 * (rand(3,1) - 0.5)) 0.3 * (rand(3,1) - 0.5); 0 0 0 1];
    end
end


host_pids_each_fid = round(size(xyz,1) / pose_num);

pid_range = 1 : host_pids_each_fid : size(xyz,1);
pid_range = [pid_range size(xyz,1)];

pid_range = unique(pid_range);

fid_to_pid_range = {};
pid_to_host_fid = [];
for i = 1 : length(pid_range)-1
    if i ~= length(pid_range)-1
        fid_to_pid_range{i,1} = [pid_range(i) : pid_range(i+1)-1];
    else
        fid_to_pid_range{i,1} = [pid_range(i) : pid_range(i+1)];
    end
    pid_to_host_fid = [pid_to_host_fid;i * ones(length(fid_to_pid_range{i,1}),1)];
end

[~,depths] = NormalizeVector(xyz);

% idps = 1./depths;

idps_gt = {};
idps = {};
pose_wb{1,1} = eye(4);
trans_ang_err{1,1} = [0 0];
for i = 1 : pose_num %size(xyz,1)
    if i == 1
        pose_wb{i,1} = eye(4);
    else
        pose_wb{i,1} = [rodrigues(0.2*(rand(3,1)-0.5)) 0.5*(rand(3,1)-0.5); 0 0 0 1];
    end
    
    for j = 1 : cam_num
        pose_cw = inv(pose_wb{i,1} * Tbc{j});
        bearing = pose_cw(1:3,1:3)*xyz' + repmat(pose_cw(1:3,4), 1, size(xyz,1));
        [bearing_, dist] = NormalizeVector(bearing');
        bearings{j,i} = bearing_;
        if j == 1
            idps{1, i} = 1./dist + 0.0001;
            idps_gt{1, i} = 1./dist;
            if ~add_point_noise
                idps{1, i} = 1./dist;
                bearings{j,i} .* repmat(dist,1,3);
            end
        end
    end
end


opt_fids = fixed_pose_num+1:pose_num;%size(xyz,1);

pose_Twb_gt = pose_wb;
pose_Twb = pose_wb;
for k = 1 : fixed_pose_num
    pose_err{k,1} = eye(4);
    pose_Twb{k,1} = (pose_wb{k,1});
end
for i = opt_fids
    pose_err{i,1} = [rodrigues(0.1*(rand(3,1)-0.5)) 0.05*(rand(3,1)-0.5); 0 0 0 1];
    if add_pose_noise
        pose_Twb{i,1} = (pose_wb{i,1}) * pose_err{i,1};
    else
        pose_Twb{i,1} = (pose_wb{i,1});
    end
end


[xMat, yMat] = meshgrid(1 : pose_num, 1 :pose_num);
pairs = [xMat(:) yMat(:)];

if opt_scale
    scale = 2;
else
    scale = 1;
end
pose_size = 6;
if use_5dof
    TbcDim = 5;
else
    TbcDim = 6;
end
if opt_Tbc
    extra_dim = 1 + cam_num * TbcDim;
else
    extra_dim = 1;
end
start_idp_idx = extra_dim + pose_size * (pose_num - fixed_pose_num);
start_scale_idx = 0; % pose_size*(pose_num-fixed_pose_num) + size(xyz,1);
loss_vec = [];
Tbc_gt = Tbc;
for iter = 1 : 30
    
    H = zeros(pose_size*(pose_num-fixed_pose_num) + size(xyz,1) + extra_dim, pose_size*(pose_num-fixed_pose_num) + size(xyz,1) + extra_dim);
    b = zeros(pose_size*(pose_num-fixed_pose_num) + size(xyz,1) + extra_dim, 1);
    err_sum = 0;
    err_count = 0;
    for pid = 1 : size(xyz,1)
        pt = xyz(pid,:);
        host_fid = pid_to_host_fid(pid);
        host_idp = idps{1, host_fid}(pid);
        host_bearing = bearings{1,host_fid}(pid,:);
        pose_wc_host = pose_Twb{host_fid,1} * Tbc{1};
        for target_fid = 1 : pose_num
            if target_fid == host_fid
                continue;
            end
            for cid = 1 : cam_num
                if cid == 1
                    %                     continue;
                end
                target_bearing = bearings{cid, target_fid}(pid,:);
                pose_cw_target = inv(pose_Twb{target_fid,1} * Tbc{cid});
                Tth = pose_cw_target * pose_wc_host;
                reproj = Tth(1:3,1:3)*(host_bearing'/host_idp) +Tth(1:3,4);
                if ~add_pose_noise && ~add_point_noise
                    assert(norm(target_bearing' - reproj./norm(reproj)) < 0.000001);
                end
                host_cid = 1 - 1;
                target_cid = cid - 1;
                if host_cid == target_cid
                    is_same_cid = true;
                else
                   is_same_cid = false; 
                end
                [err, d_err_d_Twb1, d_err_d_Twb2, d_err_d_scale, d_err_d_rho, d_err_d_Tbcx, d_err_d_Tbcy] = computeScaledReprojFactor(is_same_cid, pose_Twb{host_fid,1}, pose_Twb{target_fid,1}, Tbc{1}, Tbc{cid}, host_idp, scale, host_bearing', target_bearing');
                if ~opt_scale
                    d_err_d_scale = zeros(size(d_err_d_scale));
                end
                
                if fix_Tbc0 && host_cid == 0
                    d_err_d_Tbcx = zeros(size(d_err_d_Tbcx));
                end
                if fix_Tbc0 && target_cid == 0
                    d_err_d_Tbcy = zeros(size(d_err_d_Tbcy));
                end
                
                
                err_sum = err_sum + 235 * norm(err);
                err_count = err_count + 1;
                
                host_fid_valid = host_fid > fixed_pose_num;
                target_fid_valid = target_fid > fixed_pose_num;
                
                
                if opt_Tbc
                    H = FillMatrix(H, 1 + host_cid * TbcDim + 1, 1 + host_cid * TbcDim + 1, TbcDim, TbcDim, d_err_d_Tbcx' * d_err_d_Tbcx);
                    b = FillMatrix(b, 1 + host_cid * TbcDim + 1, 1, TbcDim, 1, d_err_d_Tbcx' * err);
                    H = FillMatrix(H, 1 + target_cid * TbcDim + 1, 1 + target_cid * TbcDim + 1, TbcDim, TbcDim, d_err_d_Tbcy' * d_err_d_Tbcy);
                    b = FillMatrix(b, 1 + target_cid * TbcDim + 1, 1, TbcDim, 1, d_err_d_Tbcy' * err);
                    
                    H = FillMatrix(H, 1 + host_cid * TbcDim + 1, 1 + target_cid * TbcDim + 1, TbcDim, TbcDim, d_err_d_Tbcx' * d_err_d_Tbcy);
                    H = FillMatrix(H, 1 + target_cid * TbcDim + 1, 1 + host_cid * TbcDim + 1, TbcDim, TbcDim, d_err_d_Tbcy' * d_err_d_Tbcx);
                    
                    H = FillMatrix(H, 1 + host_cid * TbcDim + 1, start_idp_idx + pid, TbcDim, 1, d_err_d_Tbcx' * d_err_d_rho);
                    H = FillMatrix(H, start_idp_idx + pid, 1 + host_cid * TbcDim + 1, 1, TbcDim, d_err_d_rho' * d_err_d_Tbcx);
                    H = FillMatrix(H, 1 + target_cid * TbcDim + 1, start_idp_idx + pid, TbcDim, 1, d_err_d_Tbcy' * d_err_d_rho);
                    H = FillMatrix(H, start_idp_idx + pid, 1 + target_cid * TbcDim + 1, 1, TbcDim, d_err_d_rho' * d_err_d_Tbcy);
                    
                    H = FillMatrix(H, 1 + host_cid * TbcDim + 1, 1, TbcDim, 1, d_err_d_Tbcx' * d_err_d_scale);
                    H = FillMatrix(H, 1, 1 + host_cid * TbcDim + 1, 1, TbcDim, d_err_d_scale' * d_err_d_Tbcx);
                    H = FillMatrix(H, 1 + target_cid * TbcDim + 1, 1, TbcDim, 1, d_err_d_Tbcy' * d_err_d_scale);
                    H = FillMatrix(H, 1, 1 + target_cid * TbcDim + 1, 1, TbcDim, d_err_d_scale' * d_err_d_Tbcy);
                    
                    if host_fid_valid
                        H = FillMatrix(H, 1 + host_cid * TbcDim + 1, (host_fid - fixed_pose_num - 1) * pose_size + 1 + extra_dim, TbcDim, pose_size, d_err_d_Tbcx' * d_err_d_Twb1);
                        H = FillMatrix(H, (host_fid - fixed_pose_num - 1) * pose_size + 1 + extra_dim, 1 + host_cid * TbcDim + 1, pose_size, TbcDim, d_err_d_Twb1' * d_err_d_Tbcx);
                        
                        H = FillMatrix(H, 1 + target_cid * TbcDim + 1, (host_fid - fixed_pose_num - 1) * pose_size + 1 + extra_dim, TbcDim, pose_size, d_err_d_Tbcy' * d_err_d_Twb1);
                        H = FillMatrix(H, (host_fid - fixed_pose_num - 1) * pose_size + 1 + extra_dim, 1 + target_cid * TbcDim + 1, pose_size, TbcDim, d_err_d_Twb1' * d_err_d_Tbcy);
                    end
                    
                    if target_fid_valid
                        H = FillMatrix(H, 1 + host_cid * TbcDim + 1, (target_fid - fixed_pose_num - 1) * pose_size + 1 + extra_dim, TbcDim, pose_size, d_err_d_Tbcx' * d_err_d_Twb2);
                        H = FillMatrix(H, (target_fid - fixed_pose_num - 1) * pose_size + 1 + extra_dim, 1 + host_cid * TbcDim + 1, pose_size, TbcDim, d_err_d_Twb2' * d_err_d_Tbcx);
                        
                        H = FillMatrix(H, 1 + target_cid * TbcDim + 1, (target_fid - fixed_pose_num - 1) * pose_size + 1 + extra_dim, TbcDim, pose_size, d_err_d_Tbcy' * d_err_d_Twb2);
                        H = FillMatrix(H, (target_fid - fixed_pose_num - 1) * pose_size + 1 + extra_dim, 1 + target_cid * TbcDim + 1, pose_size, TbcDim, d_err_d_Twb2' * d_err_d_Tbcy);
                    end
                end
                
                
                if host_fid_valid
                    H = FillMatrix(H, (host_fid - fixed_pose_num - 1) * pose_size + 1 + extra_dim, (host_fid - fixed_pose_num - 1) * pose_size + 1 + extra_dim, pose_size, pose_size, d_err_d_Twb1' * d_err_d_Twb1);
                    b = FillMatrix(b, (host_fid - fixed_pose_num - 1) * pose_size + 1 + extra_dim, 1, pose_size, 1, d_err_d_Twb1' * err);
                    
                    H = FillMatrix(H, (host_fid - fixed_pose_num - 1) * pose_size + 1 + extra_dim, start_idp_idx + pid, pose_size, 1, d_err_d_Twb1' * d_err_d_rho);
                    H = FillMatrix(H, start_idp_idx + pid, (host_fid - fixed_pose_num - 1) * pose_size + 1 + extra_dim, 1, pose_size, d_err_d_rho' * d_err_d_Twb1);
                    
                    H = FillMatrix(H, (host_fid - fixed_pose_num - 1) * pose_size + 1 + extra_dim, start_scale_idx + 1, pose_size, 1, d_err_d_Twb1' * d_err_d_scale);
                    H = FillMatrix(H, start_scale_idx + 1, (host_fid - fixed_pose_num - 1) * pose_size + 1 + extra_dim, 1, pose_size, d_err_d_scale' * d_err_d_Twb1);
                end
                if target_fid_valid
                    H = FillMatrix(H, (target_fid - fixed_pose_num - 1) * pose_size + 1 + extra_dim, (target_fid - fixed_pose_num - 1) * pose_size + 1 + extra_dim, pose_size, pose_size, d_err_d_Twb2' * d_err_d_Twb2);
                    b = FillMatrix(b, (target_fid - fixed_pose_num - 1) * pose_size + 1 + extra_dim, 1, pose_size, 1, d_err_d_Twb2' * err);
                    
                    H = FillMatrix(H, (target_fid - fixed_pose_num - 1) * pose_size + 1 + extra_dim, start_idp_idx + pid, pose_size, 1, d_err_d_Twb2' * d_err_d_rho);
                    H = FillMatrix(H, start_idp_idx + pid, (target_fid - fixed_pose_num - 1) * pose_size + 1 + extra_dim, 1, pose_size, d_err_d_rho' * d_err_d_Twb2);
                    
                    H = FillMatrix(H, (target_fid - fixed_pose_num - 1) * pose_size + 1 + extra_dim, start_scale_idx + 1, pose_size, 1, d_err_d_Twb2' * d_err_d_scale);
                    H = FillMatrix(H, start_scale_idx + 1, (target_fid - fixed_pose_num - 1) * pose_size + 1 + extra_dim, 1, pose_size, d_err_d_scale' * d_err_d_Twb2);
                end
                H = FillMatrix(H, start_idp_idx + pid, start_idp_idx + pid, 1, 1, d_err_d_rho' * d_err_d_rho);
                b = FillMatrix(b, start_idp_idx + pid, 1, 1, 1, d_err_d_rho' * err);
                H = FillMatrix(H, start_scale_idx + 1, start_scale_idx + 1, 1, 1, d_err_d_scale' * d_err_d_scale);
                b = FillMatrix(b, start_scale_idx + 1, 1, 1, 1, d_err_d_scale' * err);
                
                H = FillMatrix(H, start_idp_idx + pid, start_scale_idx + 1, 1, 1, d_err_d_rho' * d_err_d_scale);
                H = FillMatrix(H, start_scale_idx + 1, start_idp_idx + pid, 1, 1, d_err_d_scale' * d_err_d_rho);
                if host_fid_valid && target_fid_valid
                    H = FillMatrix(H, (host_fid - fixed_pose_num - 1) * pose_size + 1 + extra_dim, (target_fid - fixed_pose_num - 1) * pose_size + 1 + extra_dim, pose_size, pose_size, d_err_d_Twb1' * d_err_d_Twb2);
                    H = FillMatrix(H, (target_fid - fixed_pose_num - 1) * pose_size + 1 + extra_dim, (host_fid - fixed_pose_num - 1) * pose_size + 1 + extra_dim, pose_size, pose_size, d_err_d_Twb2' * d_err_d_Twb1);
                end
            end
        end
    end
    H11 = H(1:start_idp_idx,1:start_idp_idx);
    H12 = H(1:start_idp_idx,start_idp_idx+1:end);
    H22 = H(start_idp_idx+1:end,start_idp_idx+1:end);
    b1 = b(1:start_idp_idx);
    b2 = b(start_idp_idx+1:end);
    
    loss_vec = [loss_vec; err_sum];
    
    fprintf(sprintf('iter: %d, err_sum: %.15f, err_count: %d, mean_err: %.15f, scale: %.15f\n', iter, err_sum, err_count, err_sum / err_count, scale));
    if 1
        if add_point_noise
            if opt_scale
                dx = -inv(H) * b;
            else
                dx = -inv(H(2:end, 2:end)) * b(2:end);
                dx = [0;dx];
            end
        else
            if opt_scale
                dx = -inv(H(1:start_idp_idx,1:start_idp_idx)) * b(1:start_idp_idx);
            else
                dx = -inv(H(2:start_idp_idx,2:start_idp_idx)) * b(2:start_idp_idx);
                dx = [0;dx];
            end
        end
        dp = reshape(dx(extra_dim + 1:pose_size * (pose_num - fixed_pose_num) + extra_dim), 6, []);
        for id = 1 : size(dp,2)
            pose_id = id + fixed_pose_num;
            pose_Twb{pose_id,1} = Exp(dp(1:3,id), dp(4:6,id)) * pose_Twb{pose_id,1};
        end
        if add_point_noise
            for point_id = 1 : size(xyz,1)
                host_frame_id = pid_to_host_fid(point_id);
                idps{1, host_frame_id}(point_id) = idps{1, host_frame_id}(point_id) + dx(start_idp_idx + point_id);
            end
        end
        scale = scale + dx(start_scale_idx + 1);
        if opt_Tbc
            d_bc = reshape(dx(2:extra_dim), TbcDim, []);
            if use_5dof
                for idx = 1 : cam_num
                    dT = Exp(d_bc(1:3,idx), zeros(3, 1));
                    Tbc{idx,1}(1:3,1:3) = dT(1:3,1:3) * Tbc{idx,1}(1:3,1:3);
                    trans = Tbc{idx,1}(1:3,4);
                    Tbc{idx,1}(1:3,4) = rodrigues(ProduceOtherOthogonalBasis(trans) * d_bc(4:5,idx)) * trans;
                end
            else
                for idx = 1 : cam_num
                    Tbc{idx,1} = Exp(d_bc(1:3,idx), d_bc(4:6,idx)) * Tbc{idx,1};
                end
            end
            
        end
    else
        
    end
end

figure,plot(diff(loss_vec));

end
function result =  Exp(w, v)


omega = w;
theta_sq = sum(omega.^2);

if (theta_sq < 1e-10)
    theta = 0;
else
    theta = sqrt(theta_sq);
end
so3 = rodrigues(omega);
Omega = SkewSymMat(omega);
Omega_sq = Omega * Omega;
V = zeros(3, 3);

if (theta < 1e-5)
    V = so3;
    
else
    theta_sq = theta * theta;
    V = (eye(3) + ((1) - cos(theta)) / (theta_sq)*Omega + (theta - sin(theta)) / (theta_sq * theta) * Omega_sq);
end

tran = V * v;
result = eye(4);

result(1:3,1:3) = so3;
result(1:3,4) = tran;

end
function H = FillMatrix(H, start_row, start_col, row_size, col_size, value)

if (start_row < 1 || start_col < 1)
    return
end

H(start_row:start_row+row_size-1, start_col:start_col+col_size-1) = H(start_row:start_row+row_size-1, start_col:start_col+col_size-1) + value;
end
function [err, d_err_d_Twb1, d_err_d_Twb2, d_err_d_scale, d_err_d_rho, d_err_d_Tbcx, d_err_d_Tbcy] = computeScaledReprojFactor(is_same_cid, Twb1, Twb2, Tbcx, Tbcy, rho, scale, host_bearing, target_bearing)
global use_5dof
err = []; d_err_d_Twb1 = []; d_err_d_Twb2 = []; d_err_d_scale = []; d_err_d_rho = [];
d_err_d_Tbcx = [];  d_err_d_Tbcy = [];

Rwb1 = Twb1(1:3,1:3);
twb1 = Twb1(1:3,4);
Rwb2 = Twb2(1:3,1:3);
twb2 = Twb2(1:3,4);
Rbcx = Tbcx(1:3,1:3);
tbcx = Tbcx(1:3,4);
Rbcy = Tbcy(1:3,1:3);
tbcy = Tbcy(1:3,4);

Tb2b1 = inv(Twb2) * Twb1;

Tb2b1_scaled = Tb2b1;
Tb2b1_scaled(1:3,4) = scale * Tb2b1(1:3,4);

Tth = inv(Tbcy) * Tb2b1_scaled * Tbcx;
X = Tth(1:3,1:3) * host_bearing + rho * Tth(1:3,4);
err = X./norm(X) - target_bearing;

d_bearing_d_X = compute_d_bearing_d_pt_jac(X);

d_X_d_Rwb1 = -Rbcy' * Rwb2' * SkewSymMat(Rwb1 * Rbcx * host_bearing) - rho * Rbcy' * (Rwb2' * SkewSymMat(Rwb1 * tbcx + scale * twb1));
d_X_d_Rwb2 =  Rbcy' * Rwb2' * SkewSymMat(Rwb1 * Rbcx * host_bearing) + rho * Rbcy' * (Rwb2' * SkewSymMat(Rwb1 * tbcx + scale * twb1));

d_X_d_twb1 = scale * rho * Rbcy' * Rwb2';
d_X_d_twb2 = -scale * rho * Rbcy' * Rwb2';

d_X_d_rho = Rbcy' * (Rwb2' * (Rwb1 * tbcx + scale * (twb1 - twb2)) - tbcy);

d_X_d_scale = rho * Rbcy' * Rwb2' * (twb1 - twb2);

d_X_d_Twb1 = [d_X_d_Rwb1 d_X_d_twb1];
d_X_d_Twb2 = [d_X_d_Rwb2 d_X_d_twb2];


d_err_d_Twb1 = d_bearing_d_X * d_X_d_Twb1;
d_err_d_Twb2 = d_bearing_d_X * d_X_d_Twb2;
d_err_d_scale = d_bearing_d_X * d_X_d_scale;
d_err_d_rho = d_bearing_d_X * d_X_d_rho;





A = Rbcy' * Rwb2' * Rwb1;
d_X_d_Rbcx = -A * SkewSymMat(Rbcx * host_bearing) - rho * A * SkewSymMat(tbcx);


A = Rwb2' * Rwb1 * Rbcx * host_bearing;
B = Rwb2' * (Rwb1 * tbcx + scale * (twb1 - twb2));
d_X_d_Rbcy = Rbcy' * SkewSymMat(A + rho * B);

d_X_d_tbcx = rho * Rbcy' * Rwb2' * Rwb1;
d_X_d_tbcy = -rho * Rbcy';

d_X_d_Tbcx_6dof = zeros(3, 6);
d_X_d_Tbcy_6dof = zeros(3, 6);
d_X_d_Tbcx_6dof(:, 1:3) = d_X_d_Rbcx;
d_X_d_Tbcx_6dof(:, 4:6) = d_X_d_tbcx;
d_X_d_Tbcy_6dof(:,1:3) = d_X_d_Rbcy;
d_X_d_Tbcy_6dof(:,4:6) = d_X_d_tbcy;

if use_5dof
    transform_jac_x = eye(6);
    transform_jac_y = eye(6);
    transform_jac_x(4:6,1:3) = SkewSymMat(tbcx);
    transform_jac_x(4:6,4:6) = -SkewSymMat(tbcx);
    transform_jac_y(4:6,1:3) = SkewSymMat(tbcy);
    transform_jac_y(4:6,4:6) = -SkewSymMat(tbcy);
    
    d_X_d_Tbcx_6dof = d_X_d_Tbcx_6dof * transform_jac_x;
    d_X_d_Tbcy_6dof = d_X_d_Tbcy_6dof * transform_jac_y;
    d_X_d_Tbcx = zeros(3, 5);
    d_X_d_Tbcy = zeros(3, 5);
    d_X_d_Tbcx(:, 1:3) = d_X_d_Tbcx_6dof(:, 1:3);
    d_X_d_Tbcx(:, 4:5) = d_X_d_Tbcx_6dof(:, 4:6) * ProduceOtherOthogonalBasis(tbcx);
    d_X_d_Tbcy(:, 1:3) = d_X_d_Tbcy_6dof(:, 1:3);
    d_X_d_Tbcy(:, 4:5) = d_X_d_Tbcy_6dof(:, 4:6) * ProduceOtherOthogonalBasis(tbcy);
else
    d_X_d_Tbcx = d_X_d_Tbcx_6dof;
    d_X_d_Tbcy = d_X_d_Tbcy_6dof;
end
d_err_d_Tbcx = d_bearing_d_X * d_X_d_Tbcx;
d_err_d_Tbcy = d_bearing_d_X * d_X_d_Tbcy;
if (is_same_cid)
    d_err_d_Tbcx = d_err_d_Tbcx + d_err_d_Tbcy;
    d_err_d_Tbcy = zeros(size(d_err_d_Tbcy));
end


end
function d_bearing_d_pt = compute_d_bearing_d_pt_jac(pt)
d_bearing_d_pt = (norm(pt).*eye(3) - pt * pt'./norm(pt))./(norm(pt)^2);
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