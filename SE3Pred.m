function SE3Pred()
% 示例：使用 ppval 预测 SE(3) 位姿
% clear; clc;

% 生成模拟数据
[t, p, q, omega, accel] = generate_sample_data(5, 0.2);

% 初始化样条预测器
predictor = SplinePosePredictor(t, p, q, omega, accel);

% 预测未来位姿
t_future = linspace(t(end), t(end) + 1, 6)'; % 预测 1 秒，6 个点
p_pred = zeros(length(t_future), 3);
q_pred = zeros(length(t_future), 4);
for i = 1:length(t_future)
    [p_pred(i, :), q_pred(i, :)] = predict_pose(predictor, t_future(i));
end

% 可视化
visualize_trajectory(t, p, q, t_future, p_pred, q_pred);

% 生成模拟数据
function [t, p, q, omega, accel] = generate_sample_data(n_points, dt)
    t = (0:dt:(n_points-1)*dt)';
    p = [t, sin(t), cos(t)]; % 位置
    angles = deg2rad(t * 20); % 绕 Z 轴旋转
    q = quatFromEuler([zeros(n_points, 2), angles]); % 四元数 [w, x, y, z]
    omega = [zeros(n_points, 2), 20 * ones(n_points, 1) * pi/180]; % 角速度
    accel = [ones(n_points, 1), cos(t), -sin(t)]; % 加速度
end

% 样条预测器初始化
function predictor = SplinePosePredictor(t, p, q, omega, accel)
    predictor.t = t;
    predictor.dt = t(2) - t(1);
    
    % 位置样条
    predictor.splines_p = cell(3, 1);
    v_imu = cumtrapz(t, accel);
    p_imu = cumtrapz(t, v_imu);
    for a = 1:3
        p_fused = 0.8 * p(:, a) + 0.2 * p_imu(:, a); % 融合
        predictor.splines_p{a} = spline(t, p_fused);
    end
    
    % 旋转样条
    rotvecs = rotvecFromQuat(q);
    predictor.splines_rotvec = cell(3, 1);
    rotvec_imu = cumtrapz(t, omega);
    for a = 1:3
        rotvec_fused = 0.8 * rotvecs(:, a) + 0.2 * rotvec_imu(:, a);
        predictor.splines_rotvec{a} = spline(t, rotvec_fused);
    end
    
    predictor.omega = omega;
    predictor.accel = accel;
end

% 预测位姿
function [p_pred, q_pred] = predict_pose(predictor, t_new)
    % 位置预测
    p_pred = zeros(3, 1);
    for a = 1:3
        p_pred(a) = ppval(predictor.splines_p{a}, t_new);
    end
    t_last = predictor.t(end);
    if t_new > t_last
        dt = t_new - t_last;
        v_last = zeros(3, 1);
        for b = 1:3
            v_last(b) = (ppval(predictor.splines_p{b}, t_last + eps) - ...
                         ppval(predictor.splines_p{b}, t_last - eps)) / (2 * eps);
        end
        a_last = predictor.accel(end, :)';
        p_pred = p_pred + v_last * dt + 0.5 * a_last * dt^2;
    end
    
    % 旋转预测
    rotvec_new = zeros(3, 1);
    for c = 1:3
        rotvec_new(c) = ppval(predictor.splines_rotvec{c}, t_new);
    end
    if t_new > t_last
        dt = t_new - t_last;
        omega_last = predictor.omega(end, :)';
        rotvec_new = rotvec_new + omega_last * dt;
    end
    q_pred = quatFromRotvec(rotvec_new);
    q_pred = q_pred / norm(q_pred);
end

% 可视化
function visualize_trajectory(t, p, q, t_future, p_pred, q_pred)
    figure;
    hold on;
    plot3(p(:, 1), p(:, 2), p(:, 3), 'b.-', 'LineWidth', 1.5, 'DisplayName', '历史轨迹');
    plot3(p_pred(:, 1), p_pred(:, 2), p_pred(:, 3), 'r.-', 'LineWidth', 1.5, 'DisplayName', '预测轨迹');
    
    step = max(1, floor(length(t) / 3));
    for a = 1:step:length(t)
        R = quat2rotm([q(a, 1), q(a, 2:4)]);
        dir = R * [0; 0; 0.5];
        quiver3(p(a, 1), p(a, 2), p(a, 3), dir(1), dir(2), dir(3), 'b', 'LineWidth', 1.5);
    end
    for b = 1:length(t_future)
        R = quat2rotm([q_pred(b, 1), q_pred(b, 2:4)]);
        dir = R * [0; 0; 0.5];
        quiver3(p_pred(b, 1), p_pred(b, 2), p_pred(b, 3), dir(1), dir(2), dir(3), 'r', 'LineWidth', 1.5);
    end
    
    grid on; xlabel('X'); ylabel('Y'); zlabel('Z');
    title('SE(3) 位姿预测'); legend('show'); axis equal;
    hold off;
end

% 辅助函数
function q = quatFromEuler(euler)
    r = eul2quat(euler, 'ZYX');
    q = [r(:, 4), r(:, 1:3)]; % [w, x, y, z]
end

function rotvec = rotvecFromQuat(q)
    r = quat2rotm(q);
    rotvec = zeros(size(q, 1), 3);
    for a = 1:size(q, 1)
        rot = rotm2axang(r(:, :, a));
        rotvec(a, :) = rot(1:3) * rot(4);
    end
end

function q = quatFromRotvec(rotvec)
    theta = norm(rotvec);
    if theta < eps
        q = [1, 0, 0, 0]';
    else
        axis = rotvec / theta;
        q = [cos(theta/2); axis * sin(theta/2)];
    end
end
end