function predictPose()

% B-样条外插预测位姿（融合 IMU 数据并可视化）
clear; clc;

% 生成模拟数据（包括 IMU）
[t, p, q, omega, accel] = generate_sample_data(10);

% 初始化 B-样条预测器
predictor = BSplinePosePredictor(t, p, q, omega, accel, 3);

% 外插预测
t_future = linspace(10, 12, 5)'; % 预测未来 5 个点
p_pred = zeros(length(t_future), 3);
q_pred = zeros(length(t_future), 4);
for i = 1:length(t_future)
    [p_pred(i, :), q_pred(i, :)] = predict_pose(predictor, t_future(i));
    fprintf('时间 %.1f:\n', t_future(i));
    fprintf('  位置: [%.4f, %.4f, %.4f]\n', p_pred(i, 1), p_pred(i, 2), p_pred(i, 3));
    fprintf('  四元数: [%.4f, %.4f, %.4f, %.4f]\n', q_pred(i, 1), q_pred(i, 2), q_pred(i, 3), q_pred(i, 4));
end

% 可视化
visualize_trajectory(t, p, q, t_future, p_pred, q_pred);



end
% 生成模拟数据函数（包括 IMU）
function [t, p, q, omega, accel] = generate_sample_data(n_points)
    t = linspace(0, 10, n_points)';  % 时间戳
    p = [sin(t), cos(t), t];         % 位置 [x, y, z]
    angles = deg2rad(t);             % 绕 Z 轴旋转（弧度）
    q = quatFromEuler([zeros(n_points, 2), angles]); % 四元数 [w, x, y, z]
    
    % IMU 数据
    dt = t(2) - t(1);
    omega = [zeros(n_points, 2), ones(n_points, 1) * deg2rad(1)]; % 角速度 [wx, wy, wz]
    accel = [cos(t), -sin(t), ones(n_points, 1)]; % 加速度 [ax, ay, az]
end

% B-样条预测器初始化函数（融合 IMU）
function predictor = BSplinePosePredictor(t, p, q, omega, accel, degree)
    predictor.t = t;
    predictor.degree = degree;
    predictor.dt = t(2) - t(1);
    
    % 扩展节点向量以支持外插
    t_extended = [linspace(t(1) - degree * predictor.dt, t(1), degree + 1)'; 
                  t; 
                  linspace(t(end), t(end) + degree * predictor.dt, degree + 1)'];
    t_extended = unique(t_extended);
    
    % 位置 B-样条拟合（考虑加速度约束）
    predictor.splines_p = cell(3, 1);
    v = cumtrapz(t, accel); % 通过加速度积分估计速度
    p_imu = cumtrapz(t, v); % 通过速度积分估计位置
    for i = 1:3
        % 融合 IMU 和观测位置
        p_fused = 0.7 * p(:, i) + 0.3 * p_imu(:, i); % 加权融合
        predictor.splines_p{i} = spline(t, p_fused);
    end
    
    % 姿态 B-样条拟合（考虑角速度约束）
    rotvecs = rotvecFromQuat(q); % 四元数转旋转向量
    predictor.splines_rotvec = cell(3, 1);
    for i = 1:3
        % IMU 角速度积分得到的旋转向量变化
        rotvec_imu = cumtrapz(t, omega(:, i));
        rotvec_fused = 0.7 * rotvecs(:, i) + 0.3 * rotvec_imu; % 加权融合
        predictor.splines_rotvec{i} = spline(t, rotvec_fused);
    end
end

% 预测位姿函数（融合 IMU 外推）
function [p_pred, q_pred] = predict_pose(predictor, t_new)
    % 位置预测
    p_pred = zeros(3, 1);
    for i = 1:3
        p_pred(i) = ppval(predictor.splines_p{i}, t_new);
    end
    
    % IMU 外推修正（假设匀加速）
    t_last = predictor.t(end);
    if t_new > t_last
        dt = t_new - t_last;
        v_last = (ppval(predictor.splines_p{1}, t_last + eps) - ppval(predictor.splines_p{1}, t_last - eps)) / (2 * eps);
        a_last = (ppval(predictor.splines_p{1}, t_last + eps) - 2 * ppval(predictor.splines_p{1}, t_last) + ...
                  ppval(predictor.splines_p{1}, t_last - eps)) / (eps^2);
        p_pred(1) = p_pred(1) + v_last * dt + 0.5 * a_last * dt^2; % 仅 x 轴示例，实际需扩展
    end
    
    % 姿态预测
    rotvec_new = zeros(3, 1);
    for i = 1:3
        rotvec_new(i) = ppval(predictor.splines_rotvec{i}, t_new);
    end
    q_pred = quatFromRotvec(rotvec_new);
    q_pred = q_pred / norm(q_pred); % 归一化
end

% 可视化函数
function visualize_trajectory(t, p, q, t_future, p_pred, q_pred)
    figure;
    hold on;
    
    % 绘制位置轨迹
    plot3(p(:, 1), p(:, 2), p(:, 3), 'b.-', 'LineWidth', 1.5, 'DisplayName', '原始轨迹');
    plot3(p_pred(:, 1), p_pred(:, 2), p_pred(:, 3), 'r.-', 'LineWidth', 1.5, 'DisplayName', '预测轨迹');
    
    % 绘制姿态（每隔几个点绘制一个箭头）
    step = max(1, floor(length(t) / 5));
    for i = 1:step:length(t)
        rot = quat2rotm([q(i, 4), q(i, 2:3), q(i, 1)]'); % [w, x, y, z] 转旋转矩阵
        dir = rot * [0; 0; 0.5]; % Z 轴方向，长度 0.5
        quiver3(p(i, 1), p(i, 2), p(i, 3), dir(1), dir(2), dir(3), 'b', 'LineWidth', 1.5);
    end
    for i = 1:length(t_future)
        rot = quat2rotm([q_pred(i, 4), q_pred(i, 2:3), q_pred(i, 1)]');
        dir = rot * [0; 0; 0.5];
        quiver3(p_pred(i, 1), p_pred(i, 2), p_pred(i, 3), dir(1), dir(2), dir(3), 'r', 'LineWidth', 1.5);
    end
    
    % 设置图形属性
    grid on;
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('位姿轨迹与预测');
    legend('show');
    axis equal;
    hold off;
end

% 辅助函数：欧拉角转四元数
function q = quatFromEuler(euler)
    r = eul2quat(euler, 'ZYX');
    q = [r(:, 4), r(:, 1:3)]; % [w, x, y, z]
end

% 辅助函数：四元数转旋转向量
function rotvec = rotvecFromQuat(q)
    r = quat2rotm(q');
    rotvec = zeros(size(q, 1), 3);
    for i = 1:size(q, 1)
        rot = rotm2axang(r(:, :, i));
        rotvec(i, :) = rot(1:3) * rot(4);
    end
end

% 辅助函数：旋转向量转四元数
function q = quatFromRotvec(rotvec)
    theta = norm(rotvec);
    if theta < eps
        q = [1, 0, 0, 0]';
    else
        axis = rotvec / theta;
        q = [cos(theta/2); axis * sin(theta/2)]; % [w, x, y, z]
    end
end