function DualImuEKF()
global R v p w_head a_head w_carrier w_carrier_cur a_carrier imu_dt use_exact_vel

use_exact_vel = true;
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

for i = 1 : size(head_imu,1)-1
    
    if 1
        w_head = head_imu(i, 2:4)';
        a_head = head_imu(i, 5:7)';
        w_carrier = carrier_imu(i, 2:4)';
        a_carrier = carrier_imu(i, 5:7)';
        w_carrier_cur = carrier_imu(i+1, 2:4)';
    end
    if ~use_exact_vel
        ProcessOnce(relative_pose_mat_stack{i+1,1}, head_imu(i, 2:4)', head_imu(i, 5:7)', carrier_imu(i, 2:4)', carrier_imu(i, 5:7)', []);
    else
        ProcessOnce(relative_pose_mat_stack_exact{i+1,1}, head_imu(i, 2:4)', head_imu(i, 5:7)', carrier_imu(i, 2:4)', carrier_imu(i, 5:7)',  carrier_imu(i+1, 2:4)');
    end
    
end



end
function ProcessOnce(cur_state, cur_w_head, cur_a_head, cur_w_carrier, cur_a_carrier, cur_w_carrier_exact)
global R v p w_head a_head w_carrier a_carrier imu_dt use_exact_vel w_carrier_cur
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
    
    p2_est = rodrigues(-w_carrier * imu_dt) * ((eye(3) + SkewSymMat(w_carrier * imu_dt)) * p + imu_dt * v + 0.5 * imu_dt^2 * (R * a_head - a_carrier));
    v2_est = rodrigues(-w_carrier * imu_dt) * (v + SkewSymMat(w_carrier) * p + (R * a_head - a_carrier) * imu_dt) - SkewSymMat(w_carrier_cur) * p2_est;
    R2_est = rodrigues(-w_carrier * imu_dt) * R * rodrigues(w_head * imu_dt);
    R = R2_est;
    p = p2_est;
    v = v2_est;
end

temp_pose = eye(6);
temp_pose(1:3,1:3) = R;
temp_pose(1:3,4) = p;
temp_pose(1:3,6) = v;
pose_diff = cur_state - temp_pose;
if max(abs(pose_diff(:))) > 1e-10
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