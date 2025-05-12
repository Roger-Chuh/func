% =============================
% IESEKF Example with IMU
% iterated error-state extended kalman filter
% =============================
function iesekf()
clear; clc; close all;

% 初始状态：[q, v, p, b_a, b_g]
% 四元数 [x y z w] 表示初始为单位四元数
nominal_state.q = [0; 0; 0; 1];   % 初始姿态
nominal_state.v = [0; 0; 0];      % 初始速度
nominal_state.p = [0; 0; 0];      % 初始位置
nominal_state.ba = [0; 0; 0];     % 加速度偏置
nominal_state.bg = [0; 0; 0];     % 陀螺仪偏置

% 初始误差协方差矩阵
P = eye(15) * 1e-4;   % 误差状态：[dp dv dtheta dba dbg]

% 噪声参数
Qc = diag([1e-4, 1e-4, 1e-4, ...    % 角速度驱动噪声
           1e-4, 1e-4, 1e-4, ...    % 加速度驱动噪声
           1e-6, 1e-6, 1e-6, ...    % ba 漂移噪声
           1e-6, 1e-6, 1e-6]);      % bg 漂移噪声

R_gps = diag([0.1, 0.1, 0.1]);     % GPS 观测噪声
g = [0; 0; -9.81];                 % 重力加速度

% 时间设置
dt = 0.01;                         % IMU 频率 100Hz
N_steps = 1000;
t = (0:N_steps)*dt;

% 存储结果
predicted_positions = zeros(3, N_steps+1);
predicted_orientations = zeros(4, N_steps+1);
predicted_velocities = zeros(3, N_steps+1);

% 初始化
predicted_positions(:,1) = nominal_state.p;
predicted_orientations(:,1) = nominal_state.q;
predicted_velocities(:,1) = nominal_state.v;

% 模拟 IMU 数据
gyro_data = 0.01 * randn(3, N_steps);  % 模拟陀螺仪数据
acc_data = 0.01 * randn(3, N_steps);   % 模拟加速度计数据

% 模拟 GPS 数据（每10步一次）
gps_available = false(1, N_steps);
gps_data = zeros(3, N_steps);
for k = 1:10:N_steps
    gps_available(k) = true;
    gps_data(:,k) = nominal_state.p + 0.1*randn(3,1);  % GPS 带噪声
end

% =============================
% 主循环：IESEKF
% =============================
for k = 1:N_steps
    
    % 当前 IMU 输入
    gyro = gyro_data(:,k);
    acc = acc_data(:,k);
    
    % Step 1: IMU 预测（传播名义状态）
    [nominal_state_next, delta_F, delta_G] = imu_predict(nominal_state, dt, acc, gyro, g);
    
    % Step 2: 协方差传播
    Qd = delta_G * Qc * delta_G' * dt;
    P = delta_F * P * delta_F' + Qd;
    
    % Step 3: 如果有 GPS 观测，进行 Iterated 更新
    if gps_available(k)
        
        % 计算残差
        y = gps_data(:,k) - nominal_state_next.p;
        
        % 迭代更新次数
        max_iter = 5;
        
        for iter = 1:max_iter
            
            % 计算观测雅可比 H
            H = [zeros(3), zeros(3), zeros(3, 3), eye(3), zeros(3)];  % GPS 只观测位置
            
            % 卡尔曼增益
            S = H * P * H' + R_gps;
            K = P * H' / S;
            
            % 更新误差状态
            dx = K * y;
            
            % 更新名义状态
            nominal_state_next = update_nominal_state(nominal_state_next, dx);
            
            % 重新计算残差
            y = gps_data(:,k) - nominal_state_next.p;
            
            % 如果收敛则跳出迭代
            if norm(y) < 1e-4
                break;
            end
        end
        
        % 更新协方差
        P = (eye(15) - K * H) * P;
    end
    
    % 更新名义状态
    nominal_state = nominal_state_next;
    
    % 存储结果
    predicted_positions(:,k+1) = nominal_state.p;
    predicted_orientations(:,k+1) = nominal_state.q;
    predicted_velocities(:,k+1) = nominal_state.v;
    
end

% =============================
% 结果可视化
% =============================
figure;
plot3(predicted_positions(1,:), predicted_positions(2,:), predicted_positions(3,:));
hold on;
scatter3(gps_data(1,gps_available), gps_data(2,gps_available), gps_data(3,gps_available), 'filled');
grid on;
xlabel('X'); ylabel('Y'); zlabel('Z');
legend('Estimated Trajectory', 'GPS Measurements');
title('IESEKF Position Estimation with IMU and GPS');
end
function [state_next, F, G] = imu_predict(state, dt, acc, gyro, g)

q = state.q; q = q / norm(q);
v = state.v;
p = state.p;
ba = state.ba;
bg = state.bg;

% 校正测量值
omega = gyro - bg;
a = acc - ba;

state.a = a;
state.omega = omega;

% 姿态积分
dq = expm(omega_matrix(omega) * dt) * q;
q_next = dq / norm(dq);

% 转换矩阵
R = quat2rot(q_next);
a_world = R * a + g;

% 速度和位置积分
v_next = v + a_world * dt;
p_next = p + v * dt + 0.5 * a_world * dt^2;

% 更新状态
state_next.q = q_next;
state_next.v = v_next;
state_next.p = p_next;
state_next.ba = ba;
state_next.bg = bg;

% 构造 F（状态转移矩阵）
F = build_F_matrix(state, dt, a, omega, R, g);
G = build_G_matrix(state, dt, R);

end
function state = update_nominal_state(state, dx)

delta_p = dx(1:3);
delta_v = dx(4:6);
delta_theta = dx(7:9);
delta_ba = dx(10:12);
delta_bg = dx(13:15);

% 更新姿态
delta_rot = [delta_theta; 0];
state.q = quat_mult(state.q, [0.5*delta_theta; 1]);
state.q = state.q / norm(state.q);

% 更新速度和位置
state.v = state.v + delta_v;
state.p = state.p + delta_p;

% 更新偏置
state.ba = state.ba + delta_ba;
state.bg = state.bg + delta_bg;

end
function S = skew(w)
S = [0, -w(3), w(2);
     w(3), 0, -w(1);
    -w(2), w(1), 0];
end

function Omega = omega_matrix(w)
Omega = zeros(4);
Omega(1:3,1:3) = -skew(w);
Omega(1:3,4) = -w;
Omega(3+1,1:3) = w';
Omega = 0.5 * Omega;
end

function R = quat2rot(q)
q = q / norm(q);
q0 = q(4); qv = q(1:3);
R = (2*q0^2 - 1)*eye(3) + 2*q0*skew(qv) + 2*qv*qv';
end

function q = quat_mult(q1, q2)
v1 = q1(1:3); s1 = q1(4);
v2 = q2(1:3); s2 = q2(4);
v = s1*v2 + s2*v1 + cross(v1, v2);
s = s1*s2 - v1'*v2;
q = [v; s];
end

function F = build_F_matrix(state, dt, a, omega, R, g)

I3 = eye(3);
zero3 = zeros(3);

% 姿态部分
Phi_rr = I3 - 0.5*dt*skew(omega);
Phi_rv = zero3;
Phi_rp = zero3;
Phi_rba = zero3;
Phi_rbg = -0.5*dt*R * skew(a) * dt;

% 速度部分
Phi_vr = -R * skew(a) * dt;
Phi_vv = I3;
Phi_vp = zero3;
Phi_vba = -R * dt;
Phi_vbg = -Phi_vr;

% 位置部分
Phi_pr = Phi_vr * dt;
Phi_pv = I3 * dt;
Phi_pp = I3;
Phi_pba = Phi_vba * dt;
Phi_pbg = Phi_vbg * dt;

% 构建完整 F
F = [
    Phi_rr, Phi_rv, Phi_rp, Phi_rba, Phi_rbg;
    Phi_vr, Phi_vv, Phi_vp, Phi_vba, Phi_vbg;
    Phi_pr, Phi_pv, Phi_pp, Phi_pba, Phi_pbg;
    zero3, zero3, zero3, I3, zero3;
    zero3, zero3, zero3, zero3, I3
];

end

function G = build_G_matrix(state, dt, R)

I3 = eye(3);
zero3 = zeros(3);

% 噪声输入矩阵 G
G_rr = -0.5*R*dt;
G_rv = zero3;
G_rp = zero3;
G_rba = zero3;
G_rbg = -0.5*R * skew(state.a)*dt;

G_vr = -R * skew(state.a) * dt;
G_vv = zero3;
G_vp = zero3;
G_vba = -R * dt;
G_vbg = -G_vr;

G_pr = G_vr * dt;
G_pv = G_vv * dt;
G_pp = zero3;
G_pba = G_vba * dt;
G_pbg = G_vbg * dt;

G = [
    G_rr, G_rv, G_rp, G_rba, G_rbg;
    G_vr, G_vv, G_vp, G_vba, G_vbg;
    G_pr, G_pv, G_pp, G_pba, G_pbg;
    zero3, zero3, zero3, zero3, zero3;
    zero3, zero3, zero3, zero3, zero3
];

end