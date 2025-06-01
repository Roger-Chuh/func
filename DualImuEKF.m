function DualImuEKF()
global R v p w_head a_head w_carrier w_carrier_cur a_carrier imu_dt use_exact_vel cov aa_head aw_head  aa_carrier aw_carrier real_run Q R_ add_noise disable_jerk
% close all
use_exact_vel = false;
real_run = true;
add_noise = true;
disable_jerk = true;

if ~use_exact_vel
    real_run = false;
    add_noise = false;
end
head_imu = load('G:\matlab\data\direct\gt\D2_011\4\tbc\head_imu_data.txt');
carrier_imu = load('G:\matlab\data\direct\gt\D2_011\4\tbc\carrier_imu_data.txt');

state0_head = eye(6);
state0_carrier = eye(6);
state_cur_relative = eye(6);
imu_dt = 1e-3;
gravity = [-9.7964, 0, 0]';
relative_pose_stack = [];
head_pose_stack = [];
carrier_pose_stack = [];
relative_pose_mat_stack = {};
relative_pose_mat_stack_exact = {};
% generate gt
for id = 1 : size(head_imu, 1)
    gyro_carrier_cur = carrier_imu(id,2:4)';
    if id == 1
        state_cur_head = state0_head;
        state_cur_carrier = state0_carrier;
    else
        gyro_head = head_imu(id-1,2:4)';
        acc_head = head_imu(id-1,5:7)';
        dR_head = rodrigues(gyro_head * imu_dt);
        R_head = state_cur_head(1:3,1:3);
        state_cur_head(1:3,4) = state_cur_head(1:3,4) + 0.5 * (gravity + R_head * acc_head) * imu_dt^2 + state_cur_head(1:3,6) * imu_dt;
        state_cur_head(1:3,6) = state_cur_head(1:3,6) + (gravity + R_head * acc_head) * imu_dt;
        state_cur_head(1:3,1:3) = R_head * dR_head;
        
        
        gyro_carrier = carrier_imu(id-1,2:4)';
        acc_carrier = carrier_imu(id-1,5:7)';
        dR_carrier = rodrigues(gyro_carrier * imu_dt);
        R_carrier = state_cur_carrier(1:3,1:3);
        state_cur_carrier(1:3,4) = state_cur_carrier(1:3,4) + 0.5 * (gravity + R_carrier * acc_carrier) * imu_dt^2 + state_cur_carrier(1:3,6) * imu_dt;
        state_cur_carrier(1:3,6) = state_cur_carrier(1:3,6) + (gravity + R_carrier * acc_carrier) * imu_dt;
        state_cur_carrier(1:3,1:3) = R_carrier * dR_carrier;
    end
    
    Twb = inv(state_cur_carrier(1:4,1:4)) * state_cur_head(1:4,1:4);
    vel = state_cur_head(1:3,6) - state_cur_carrier(1:3,6);
    state_cur_relative(1:4,1:4) = Twb;
    state_cur_relative(1:3,6) = state_cur_carrier(1:3,1:3)' * vel;
    
    vel_exact = state_cur_carrier(1:3,1:3)' * (state_cur_head(1:3,6) - state_cur_carrier(1:3,6)) - SkewSymMat(gyro_carrier_cur) * Twb(1:3, 4);
    
    
    
    relative_pose_mat_stack{id,1} = state_cur_relative;
    state_cur_relative(1:3, 6) = vel_exact;
    relative_pose_mat_stack_exact{id, 1} = state_cur_relative;
    
    relative_pose_stack = [relative_pose_stack;[reshape(state_cur_relative(1:3,1:3), 1, 9) state_cur_relative(1:3,4)' state_cur_relative(1:3,6)']];
    
    head_pose_stack = [head_pose_stack; [reshape(state_cur_head(1:3,1:3), 1, 9) state_cur_head(1:3,4)' state_cur_head(1:3,6)']];
    carrier_pose_stack = [carrier_pose_stack; [reshape(state_cur_carrier(1:3,1:3), 1, 9) state_cur_carrier(1:3,4)' state_cur_carrier(1:3,6)']];
end


R = eye(3);
v = zeros(3, 1);
p = zeros(3, 1);

w_head = head_imu(1, 2:4)';
a_head = head_imu(1, 5:7)';

w_carrier = carrier_imu(1, 2:4)';
a_carrier = carrier_imu(1, 5:7)';

if ~real_run
    cov = 1 * eye(33,33);
    Q = 0.1 * eye(24, 24);
else
    if ~disable_jerk
        cov = 0.001 * eye(27,27);
        Q = 0.002 * eye(18, 18);
        R_ = 0.01 * eye(21, 21);
        R_(1:9, 1:9) = 0.01 * eye(9);
%         R_ = zeros(21, 21);
    else
        cov = 1 * eye(21,21);
        Q = 0.002 * eye(12, 12);
        R_ = 0.01 * eye(21, 21);
        R_(1:9, 1:9) = 0.001 * eye(9);
    end
end

Err = [];
for i = 1 : size(head_imu,1)-1
    
    if ~real_run
        w_head = head_imu(i, 2:4)';
        a_head = head_imu(i, 5:7)';
        w_carrier = carrier_imu(i, 2:4)';
        a_carrier = carrier_imu(i, 5:7)';
        w_carrier_cur = carrier_imu(i+1, 2:4)';
    end
    if ~use_exact_vel
        ProcessOnce(relative_pose_mat_stack{i+1,1}, head_imu(i, 2:4)', head_imu(i, 5:7)', carrier_imu(i, 2:4)', carrier_imu(i, 5:7)', []);
    else
        if i == 1
            aa_head = zeros(3,1);
            aw_head = (head_imu(i+1, 2:4)' - head_imu(i, 2:4)')./imu_dt;
            aa_carrier = zeros(3,1);
            aw_carrier = (carrier_imu(i+1, 2:4)' - carrier_imu(i, 2:4)')./imu_dt;
        end
        err = ProcessOnce(relative_pose_mat_stack_exact{i+1,1}, head_imu(i, 2:4)', head_imu(i, 5:7)', carrier_imu(i, 2:4)', carrier_imu(i, 5:7)',  carrier_imu(i+1, 2:4)');
        Err = [Err; err'];
    end
    
end


figure,plot(Err)
end
function err = ProcessOnce(cur_state, cur_w_head, cur_a_head, cur_w_carrier, cur_a_carrier, cur_w_carrier_next)
global R v p w_head a_head w_carrier a_carrier imu_dt use_exact_vel w_carrier_cur aw_carrier cov aw_head real_run Q R_ add_noise disable_jerk
cur_state_gt = cur_state;
cur_w_head_gt = cur_w_head;
cur_a_head_gt = cur_a_head;
cur_w_carrier_gt = cur_w_carrier;
cur_a_carrier_gt = cur_a_carrier;
% propagate
if ~use_exact_vel
    J3_head = 0.5 * imu_dt^2 * a_head;
    J2_head = imu_dt * a_head;
    J1_head = rodrigues(w_head * imu_dt);
    
    J3_carrier = 0.5 * imu_dt^2 * a_carrier;
    J2_carrier = imu_dt * a_carrier;
    J1_carrier = rodrigues(w_carrier * imu_dt);
    
    R1 = R';
    p1 = -R' * p;
    v1 = -R' * v;
    
    R2_est = J1_head' * R1 * J1_carrier;
    v2_est = J1_head' * (v1 + R1 * J2_carrier - J2_head);
    p2_est = J1_head' * (p1 + v1 * imu_dt + R1 * J3_carrier - J3_head);
    
    R = R2_est';
    p = -R2_est' * p2_est;
    v = -R2_est' * v2_est;
else
    if ~real_run
        p2_est = rodrigues(-w_carrier * imu_dt) * ((eye(3) + SkewSymMat(w_carrier * imu_dt)) * p + imu_dt * v + 0.5 * imu_dt^2 * (R * a_head - a_carrier));
        v2_est = rodrigues(-w_carrier * imu_dt) * (v + SkewSymMat(w_carrier) * p + (R * a_head - a_carrier) * imu_dt) - SkewSymMat(w_carrier_cur) * p2_est;
        R2_est = rodrigues(-w_carrier * imu_dt) * R * rodrigues(w_head * imu_dt);
    else
        if ~disable_jerk
            w_carrier_new = w_carrier + aw_carrier * imu_dt;
            w_carrier_cur_new = w_carrier + 2 * aw_carrier * imu_dt;
            
            w_carrier_dt = w_carrier * imu_dt + 0.5 * aw_carrier * imu_dt^2;
            w_head_dt = w_head * imu_dt + 0.5 * aw_head * imu_dt^2;
            p2_est = rodrigues(-w_carrier_dt) * ((eye(3) + SkewSymMat(w_carrier_dt)) * p + imu_dt * v + 0.5 * imu_dt^2 * (R * a_head - a_carrier));
            v2_est = rodrigues(-w_carrier_dt) * (v + SkewSymMat(w_carrier_new) * p + (R * a_head - a_carrier) * imu_dt) - SkewSymMat(w_carrier_cur_new) * p2_est;
            R2_est = rodrigues(-w_carrier_dt) * R * rodrigues(w_head_dt);
        else
            p2_est = rodrigues(-w_carrier * imu_dt) * ((eye(3) + SkewSymMat(w_carrier * imu_dt)) * p + imu_dt * v + 0.5 * imu_dt^2 * (R * a_head - a_carrier));
            v2_est = rodrigues(-w_carrier * imu_dt) * (v + SkewSymMat(w_carrier) * p + (R * a_head - a_carrier) * imu_dt) - SkewSymMat(w_carrier) * p2_est;
            R2_est = rodrigues(-w_carrier * imu_dt) * R * rodrigues(w_head * imu_dt);
        end
    end
    R = R2_est;
    p = p2_est;
    v = v2_est;
end

temp_pose = eye(6);
temp_pose(1:3,1:3) = R;
temp_pose(1:3,4) = p;
temp_pose(1:3,6) = v;
pose_diff = cur_state - temp_pose;
if max(abs(pose_diff(:))) > 1e-10 && ~real_run
    fprintf(sprintf('pose diff too big: %.20f\n', max(abs(pose_diff(:)))));
end


% compute cov
if ~use_exact_vel
    A_head = computeCov(w_head, a_head);
    A_carrier = computeCov(w_carrier, a_carrier);
    
    cov_A = zeros(9, 9);
    cov_B = zeros(9, 9);
    
    
    cov_A(1:3,1:3) = -R2_est';
    cov_A(4:6,1:3) = J1_head * SkewSymMat(v2_est);
    cov_A(7:9,1:3) = J1_head * SkewSymMat(p2_est);
    cov_A(4:6,4:6) = -J1_head;
    cov_A(7:9,7:9) = -J1_head;
    
    cov_B(1:3,1:3) = eye(3);
    cov_B(4:6,4:6) = J1_head * R2_est;
    cov_B(7:9,7:9) = J1_head * R2_est;
    
    cov_RVP = cov_A * A_head(1:9,1:9) * cov_A' + cov_B * A_carrier(1:9,1:9) * cov_B;
else
    
    [F, G] = computeDualCov(w_carrier, w_head, a_carrier, a_head, p, v, R, aw_carrier);
    cov = F * cov * F' + (F * G * imu_dt) * Q * (F * G * imu_dt)';
    Twb = cur_state(1:4, 1:4);
    vwb = cur_state(1:3, 6);
    if add_noise
        cur_w_head = cur_w_head + 0.01 * (rand(3,1) - 0.5);
        cur_w_carrier = cur_w_carrier + 0.01 * (rand(3,1) - 0.5);
        cur_a_head = cur_a_head + 0.1 * (rand(3,1) - 0.5);
        cur_a_carrier = cur_a_carrier + 0.1 * (rand(3,1) - 0.5);
    end
    % compute jac
    err_stack = [];
    if ~disable_jerk
        for iter = 1 : 1
            [H, err] = computeMeasurementJac(Twb, vwb, cur_w_head, cur_a_head, cur_w_carrier, cur_a_carrier);
            err_stack = [err_stack; norm(err)];
            % update state and cov
            S = H * cov * H' + R_;
            K = cov * H' * inv(S);
            dx = K * err;
            p = p + dx(1:3);
            v = v + dx(4:6);
            R = R * rodrigues(dx(7:9));
            a_carrier = a_carrier + dx(10:12);
            w_carrier = w_carrier + dx(13:15);
            a_head = a_head + dx(16:18);
            w_head = w_head + dx(19:21);
            aw_carrier = aw_carrier + dx(22:24);
            aw_head = aw_head + dx(25:27);
        end
        cov = (eye(27, 27) - K * H) * cov;
    else
        for iter = 1 : 10
            [H, err] = computeMeasurementJacNoJerk(Twb, vwb, cur_w_head, cur_a_head, cur_w_carrier, cur_a_carrier);
            err_stack = [err_stack; norm(err)];
            % update state and cov
            S = H * cov * H' + R_;
            K = cov * H' * inv(S);
            dx = K * err;
            p = p + dx(1:3);
            v = v + dx(4:6);
            R = R * rodrigues(dx(7:9));
            a_carrier = a_carrier + dx(10:12);
            w_carrier = w_carrier + dx(13:15);
            a_head = a_head + dx(16:18);
            w_head = w_head + dx(19:21);
        end
        cov = (eye(21, 21) - K * H) * cov;
    end
    
end
end
function [H, err] = computeMeasurementJac(Twb, vwb, cur_w_head, cur_a_head, cur_w_carrier, cur_a_carrier)
global R v p w_head a_head w_carrier a_carrier imu_dt aw_carrier aw_head
H = zeros(21, 27);
err = zeros(21,1);
%% state order ##### p v R a1 w1 a2 w2 aw1 aw2

%%   err order ##### p v R a1 w1 a2 w2

%% d_errP_d_p
err(1:3) = p - Twb(1:3,4);
H(1:3, 1:3) = eye(3);

%% d_errV_d_v
err(4:6) = v - vwb;
H(4:6, 4:6) = eye(3);

%% d_errR_d_r
err(7:9) = rodrigues(R' * Twb(1:3,1:3));
H(7:9, 7:9) = -JlInv(err(7:9));

%% d_erra1_d_a1
err(10:12) = a_carrier - cur_a_carrier;
H(10:12, 10:12) = eye(3);

%% d_errw1_d_w1
err(13:15) = w_carrier + aw_carrier * imu_dt - cur_w_carrier;
H(13:15, 13:15) = eye(3);
% d_errw1_d_aw1
H(13:15, 22:24) = eye(3) * imu_dt;

%% d_erra2_d_a2
err(16:18) = a_head - cur_a_head;
H(16:18, 16:18) = eye(3);

%% d_errw2_d_w2
err(19:21) = w_head + aw_head * imu_dt - cur_w_head;
H(19:21, 19:21) = eye(3);
% d_errw2_d_aw2
H(19:21, 25:27) = eye(3) * imu_dt;

H = -H;
end
function [H, err] = computeMeasurementJacNoJerk(Twb, vwb, cur_w_head, cur_a_head, cur_w_carrier, cur_a_carrier)
global R v p w_head a_head w_carrier a_carrier imu_dt aw_carrier aw_head
H = zeros(21, 21);
err = zeros(21,1);
%% state order ##### p v R a1 w1 a2 w2 aw1 aw2

%%   err order ##### p v R a1 w1 a2 w2

%% d_errP_d_p
err(1:3) = p - Twb(1:3,4);
H(1:3, 1:3) = eye(3);

%% d_errV_d_v
err(4:6) = v - vwb;
H(4:6, 4:6) = eye(3);

%% d_errR_d_r
err(7:9) = rodrigues(R' * Twb(1:3,1:3));
H(7:9, 7:9) = -JlInv(err(7:9));

%% d_erra1_d_a1
err(10:12) = a_carrier - cur_a_carrier;
H(10:12, 10:12) = eye(3);

%% d_errw1_d_w1
err(13:15) = w_carrier - cur_w_carrier;
H(13:15, 13:15) = eye(3);
% d_errw1_d_aw1
% H(13:15, 22:24) = eye(3) * imu_dt;

%% d_erra2_d_a2
err(16:18) = a_head - cur_a_head;
H(16:18, 16:18) = eye(3);

%% d_errw2_d_w2
err(19:21) = w_head - cur_w_head;
H(19:21, 19:21) = eye(3);
% d_errw2_d_aw2
% H(19:21, 25:27) = eye(3) * imu_dt;

H = -H;
end
function [F, G] = computeDualCov(w1, w2, a1, a2, p12, v12, R12, aw1)
global imu_dt real_run disable_jerk
if 0
    A = zeros(27, 27);
    
    A(1:3, 4:6) = eye(3);
    A(4:6, 1:3) = -SkewSymMat(w1) * SkewSymMat(w1) - SkewSymMat(aw1);
    A(4:6, 4:6) = -2 * SkewSymMat(w1);
    A(4:6, 7:9) = -R12 * SkewSymMat(a2);
    A(4:6, 10:12) = -eye(3);
    A(4:6, 13:15) = SkewSymMat(w1) * SkewSymMat(p12) + SkewSymMat(SkewSymMat(w1) * p12) + 2 * SkewSymMat(v12);
    A(4:6, 16:18) = R12;
    A(4:6, 22:24) = SkewSymMat(p12);
    
    A(7:9, 7:9) = -SkewSymMat(w2);
    A(7:9, 13:15) = -R12';
    A(7:9, 19:21) = eye(3);
    % A(10:12, 22:24) = eye(3);
    % A(13:15, 25:27) = eye(3);
    % A(16:18, 28:30) = eye(3);
    % A(19:21, 31:33) = eye(3);
else
    A = zeros(33, 33);
    A(1:3, 4:6) = eye(3);
    A(4:6, 1:3) = -SkewSymMat(w1) * SkewSymMat(w1) - SkewSymMat(aw1);
    A(4:6, 4:6) = -2 * SkewSymMat(w1);
    A(4:6, 7:9) = -R12 * SkewSymMat(a2);
    A(4:6, 10:12) = -eye(3);
    A(4:6, 13:15) = SkewSymMat(w1) * SkewSymMat(p12) + SkewSymMat(SkewSymMat(w1) * p12) + 2 * SkewSymMat(v12);
    A(4:6, 16:18) = R12;
    A(4:6, 25:27) = SkewSymMat(p12);
    A(7:9, 7:9) = -SkewSymMat(w2);
    A(7:9, 13:15) = -R12';
    A(7:9, 19:21) = eye(3);
    A(10:12, 22:24) = eye(3);
    A(13:15, 25:27) = eye(3);
    A(16:18, 28:30) = eye(3);
    A(19:21, 31:33) = eye(3);
    if ~real_run
        F = expm(A * imu_dt);
    else
        if ~disable_jerk
        AA = [A(1:21, [1:21 25:27 31:33]);...
            A(25:27, [1:21 25:27 31:33]);...
            A(31:33, [1:21 25:27 31:33]);];
        else
            AA = A(1:21, [1:21]);
        end
        F = expm(AA * imu_dt);
    end
    if 1
        G = zeros(33, 24);
        G(4:6, 1:3) = eye(3);
        G(4:6, 4:6) = -2 * SkewSymMat(v12) - SkewSymMat(SkewSymMat(w1) * p12) - SkewSymMat(w1) * SkewSymMat(p12);
        G(4:6, 7:9) = -R12;
        G(4:6, 16:18) = -SkewSymMat(p12);
        G(7:9,4:6) = R12';
        G(7:9,10:12) = -eye(3);
        
        G(10:12,13:15) = -eye(3);
        G(13:15,16:18) = -eye(3);
        G(16:18, 19:21) = -eye(3);
        G(19:21, 22:24) = -eye(3);
    else
        G = zeros(33, 36);
        G(4:6, 1:3) = eye(3);
        G(4:6, 4:6) = -2 * SkewSymMat(v12) - SkewSymMat(SkewSymMat(w1) * p12) - SkewSymMat(w1) * SkewSymMat(p12);
        G(4:6, 7:9) = -R12;
        G(4:6, 16:18) = -SkewSymMat(p12);
        G(7:9,4:6) = R12';
        G(7:9,10:12) = -eye(3);
        
        G(10:12,13:15) = -eye(3);
        G(13:15,16:18) = -eye(3);
        G(16:18, 19:21) = -eye(3);
        G(19:21, 22:24) = -eye(3);
        G(22:24, 25: 27) = eye(3);
        G(25: 27, 28:30) = eye(3);
        G(28:30, 31:33) = eye(3);
        G(31:33, 34:36) = eye(3);
    end
    if real_run
        if ~disable_jerk
            GG = [G(1:21, [1:12 16:18 22:24]);...
                G(25:27, [1:12 16:18 22:24]);...
                G(31:33, [1:12 16:18 22:24]);];
            G = GG;
        else
            G = G(1:21, 1:12);
        end
    end
    
end
end
function A = computeCov(gyro_vehicle, acc_vehicle)
global imu_dt
dt = imu_dt;
dRmat = rodrigues(gyro_vehicle * dt);
rightJmat = Jr(gyro_vehicle * dt);
dRmat_T = dRmat';
rightJmat_T = rightJmat';
leftJmat = Jl(gyro_vehicle * dt);
skew_a = SkewSymMat(acc_vehicle);
skew_w = SkewSymMat(gyro_vehicle);
dt2 = dt * dt;
dt3 = dt2 * dt;
A = zeros(15, 15);
A(10:12, 10:12) = eye(3);
A(13:15, 13:15) = eye(3);

A(1:3,1:3) = dRmat_T;
A(1:3, 10:12) = -(-dRmat_T * leftJmat * dt);
skew_a_skew_w = skew_a * skew_w;
skew_w_skew_a = skew_a_skew_w';


A(4:6, 1:3) = -dRmat_T * SkewSymMat(acc_vehicle * dt);

A(4:6,4:6) = dRmat_T;
A(4:6,10:12) = -(0.5 * dt2 * skew_a - (skew_a_skew_w + skew_w_skew_a) * (dt3) / 6.0);
A(4:6,13:15) = -(-dRmat_T * leftJmat * dt);


A(7:9,1:3) = -dRmat_T * SkewSymMat(0.5 * dt2 * acc_vehicle);

A(7:9, 4:6) = dt * dRmat_T;
A(7:9,7:9) = dRmat_T;
A(7:9,10:12) = -(dRmat_T * SkewSymMat(leftJmat * acc_vehicle * dt3 / 6.0));
A(7:9,13:15) = -(SkewSymMat(gyro_vehicle * dt3 / 3.0) - (dt2 / 2.0) * eye(3));

end