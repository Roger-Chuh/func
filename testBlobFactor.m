function testBlobFactor()

kUseLieAlgebra = true;
add_noise = true;
fix_pose_num = 3;

max_iter = 20;


pose_num = 20;
point_num = 20;

xyz = 10 * rand(point_num, 3);
xyz(:,3) = xyz(:,3) * 20;
xyz(:,3) = xyz(:,3) .* sign(xyz(:,3));

Tgw = [rodrigues(0.1 * (rand(3, 1)-0.5)) 10 * (rand(3, 1)-0.5); 0 0 0 1];
Tbc = [rodrigues(0.1 * (rand(3, 1)-0.5)) 10 * (rand(3, 1)-0.5); 0 0 0 1];
vms = cell(pose_num, 1);
for i = 1 : pose_num
    Twb{i, 1} = [rodrigues(0.5 * (rand(3, 1)-0.5)) 50 * (rand(3, 1)-0.5); 0 0 0 1];
    Tgc =  Tgw * Twb{i, 1} * Tbc;
    Tcg = inv(Tgc);
    vms{i,1} = [];
    for j = 1 : point_num
        reproj = Tcg(1:3,1:3) * xyz(j,:)' + Tcg(1:3,4);
        reproj = reproj'./norm(reproj);
        vms{i,1} = [vms{i,1}; reproj];
    end
end




Tgw_gt = Tgw;
Tbc_gt = Tbc;
Twb = Twb;
scale = 1;

if add_noise
    Tgw = Tgw * [rodrigues(0.1 * (rand(3, 1)-0.5)) 5 * (rand(3, 1)-0.5); 0 0 0 1];
    norm_tbc = norm(Tbc(1:3,4));
    Tbc = Tbc * [rodrigues(0.1 * (rand(3, 1)-0.5)) 5 * (rand(3, 1)-0.5); 0 0 0 1];
    Tbc(1:3,4) = Tbc(1:3,4) / norm(Tbc(1:3,4)) * norm_tbc;
    
    for i = fix_pose_num + 1 : pose_num
        Twb{i,1} = Twb{i,1} * [rodrigues(0.1 * (rand(3, 1)-0.5)) 10 * (rand(3, 1)-0.5); 0 0 0 1];
    end
    scale = 1.2;
end
loss_vec = [];
for iter = 1 : max_iter
    % scale Tgw Tbc / Twb
    H = zeros((pose_num+1) * 6 + 5 + 1 - 6 * fix_pose_num, (pose_num+1) * 6 + 5 + 1 - 6 * fix_pose_num) ;
    b = zeros((pose_num+1) * 6 + 5 + 1 - 6 * fix_pose_num, 1) ;
    
    err_sum = 0;
    err_count = 0;
    for i = 1 : pose_num
        scale_index = 1;
        Tgw_index = 2;
        Tbc_index = 8;
        Twb_index = 13 + (i-1) * 6 - 6 * fix_pose_num;
        for j = 1 : point_num
            [res, d_err_d_Twb, d_err_d_Tbc_5dof, d_err_d_Tgw, d_err_d_scale] = computeResAndJac(scale, Twb{i,1}, Tbc,Tgw, xyz(j,:)', vms{i,1}(j,:)', kUseLieAlgebra);
            err_sum = err_sum + 235 * norm(res);
            err_count = err_count + 1;
            
            H = FillMatrix(H, scale_index, scale_index, 1, 1, d_err_d_scale' * d_err_d_scale);
            H = FillMatrix(H, Tgw_index, Tgw_index, 6, 6, d_err_d_Tgw' * d_err_d_Tgw);
            H = FillMatrix(H, Tbc_index, Tbc_index, 5, 5, d_err_d_Tbc_5dof' * d_err_d_Tbc_5dof);
            if i > fix_pose_num
                H = FillMatrix(H, Twb_index, Twb_index, 6, 6, d_err_d_Twb' * d_err_d_Twb);
            end
            b = FillMatrix(b, scale_index, 1, 1, 1, d_err_d_scale' * res);
            b = FillMatrix(b, Tgw_index, 1, 6, 1, d_err_d_Tgw' * res);
            b = FillMatrix(b, Tbc_index, 1, 5, 1, d_err_d_Tbc_5dof' * res);
            if i > fix_pose_num
                b = FillMatrix(b, Twb_index, 1, 6, 1, d_err_d_Twb' * res);
            end
            % cross entry
            %
            H = FillMatrix(H, scale_index, Tgw_index, 1, 6, d_err_d_scale' * d_err_d_Tgw);
            H = FillMatrix(H, Tgw_index, scale_index, 6, 1, d_err_d_Tgw' * d_err_d_scale);
            
            H = FillMatrix(H, scale_index, Tbc_index, 1, 5, d_err_d_scale' * d_err_d_Tbc_5dof);
            H = FillMatrix(H, Tbc_index, scale_index, 5, 1, d_err_d_Tbc_5dof' * d_err_d_scale);
            if i > fix_pose_num
                H = FillMatrix(H, scale_index, Twb_index, 1, 6, d_err_d_scale' * d_err_d_Twb);
                H = FillMatrix(H, Twb_index, scale_index, 6, 1, d_err_d_Twb'* d_err_d_scale);
            end
            %
            H = FillMatrix(H, Tgw_index, Tbc_index, 6, 5, d_err_d_Tgw' * d_err_d_Tbc_5dof);
            H = FillMatrix(H, Tbc_index, Tgw_index, 5, 6, d_err_d_Tbc_5dof' * d_err_d_Tgw);
            if i > fix_pose_num
                H = FillMatrix(H, Tgw_index, Twb_index, 6, 6, d_err_d_Tgw' * d_err_d_Twb);
                H = FillMatrix(H, Twb_index,Tgw_index, 6, 6, d_err_d_Twb' * d_err_d_Tgw);
            end
            %
            if i > fix_pose_num
                H = FillMatrix(H, Tbc_index, Twb_index, 5, 6, d_err_d_Tbc_5dof' * d_err_d_Twb);
                H = FillMatrix(H, Twb_index, Tbc_index, 6, 5, d_err_d_Twb' * d_err_d_Tbc_5dof);
            end
        end
        
        
        
    end
    
    dx = -inv(H) * b;
    scale = scale + dx(1);
    
    dx_Tgw = dx(2:7);
    if ~kUseLieAlgebra
        dTgw = Exp(dx_Tgw(1:3), dx_Tgw(4:6));
        Tgw = dTgw * Tgw;
    else
        Tgw_log = LogSE3(Tgw);
        Tgw_log_new = Tgw_log + dx_Tgw;
        Tgw = Exp(Tgw_log_new(1:3), Tgw_log_new(4:6));
    end
    
    dx_Tbc = dx(8:12);
    dR = rodrigues(dx_Tbc(1:3));
    Tbc(1:3, 1:3) = dR * Tbc(1:3, 1:3);
    Ag = ProduceOtherOthogonalBasis(Tbc(1:3,4));
    t_trans_rot = Ag * dx_Tbc(4:5);
    Tbc(1:3,4) = rodrigues(t_trans_rot) * Tbc(1:3,4);
    
    for k = fix_pose_num + 1 : pose_num
        dx_Twb = dx(13 + (k-1-fix_pose_num) * 6 : 13 + (k-1-fix_pose_num) * 6 + 5);
        dTwb = Exp(dx_Twb(1:3), dx_Twb(4:6));
        Twb{k,1} = dTwb * Twb{k,1};
    end
    
    loss_vec = [loss_vec; err_sum];
    
    fprintf(sprintf('iter: %d, err_sum: %.15f, err_count: %d, mean_err: %.15f, scale: %.15f, rank: %d, size: %d\n', iter, err_sum, err_count, err_sum / err_count, scale, rank(H), size(H, 1)));
    
    
end





end


function [res, d_err_d_Twb, d_err_d_Tbc_5dof, d_err_d_Tc0h0, d_err_d_scale] = computeResAndJac(scale, Twb, Tbc,Tc0h0, xyz, bearing, kUseLieAlgebra)

d_err_d_Twb = [];
d_err_d_Tbc_5dof = [];
d_err_d_Tc0h0 = [];
d_err_d_scale = [];

Tbc_use = Tbc;
% align visual scale to imu scale
Tbc_use(1:3,4) = Tbc_use(1:3,4) / scale;
Twc = Twb * Tbc_use;
Tgc = Tc0h0 * Twc;
Tcg = inv(Tgc);
reproj = Tcg(1:3,1:3) * (xyz) + Tcg(1:3,4);
res = reproj./norm(reproj) - bearing;

d_bearing_d_X = compute_d_bearing_d_pt_jac(reproj);
d_X_d_Tcg = zeros(3, 6);
d_X_d_Tcg(:,1:3) = -SkewSymMat(reproj);
d_X_d_Tcg(:,4:6) = eye(3);
d_err_d_Tcg = d_bearing_d_X * d_X_d_Tcg;


d_err_d_Twb = d_err_d_Tcg * (-Adj(inv(Twc)));
d_err_d_Tc0h0 = d_err_d_Tcg * (-Adj(Tcg));


d_err_d_Tbc_easy = d_err_d_Tcg * (-Adj(inv(Tbc_use)));

d_err_d_Tbc_true_6dof = d_err_d_Tbc_easy;
d_err_d_Tbc_true_6dof(:,1:3) = d_err_d_Tbc_true_6dof(:,1:3) / scale;


transform_jac = eye(6);
transform_jac(4:6, 1:3) = SkewSymMat(Tbc(1:3,4));
transform_jac(4:6,4:6) = -SkewSymMat(Tbc(1:3,4));
d_err_d_Tbc_true_6dof = d_err_d_Tbc_true_6dof * transform_jac;


d_err_d_Tbc_5dof = zeros(3, 5);
d_err_d_Tbc_5dof(:,1:3) = d_err_d_Tbc_true_6dof(:,1:3);
d_err_d_Tbc_5dof(:,4:5) = d_err_d_Tbc_true_6dof(:,4:6) * ProduceOtherOthogonalBasis(Tbc(1:3,4));


d_err_d_scale = d_err_d_Tbc_easy(:,4:6) * (-Tbc(1:3,4) / scale / scale);


if kUseLieAlgebra
    d_err_d_Tc0h0 = ExpSE3_to_LieAlgebraSE3(d_err_d_Tc0h0, Tc0h0);
    
end


end