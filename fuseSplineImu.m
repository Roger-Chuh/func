function fuseSplineImu()
% 示例：使用 ppval 预测一维位置
clear; clc;
% 生成模拟数据
t = 0:0.5:2; % 时间点 [0, 0.5, 1, 1.5, 2]
p = [0, 0.25, 1, 2.25, 4]; % 位置控制点 (模拟二次函数 p = t^2)
accel = [0, 1, 2, 3, 4]; % 加速度数据
% 拟合三次样条
pp = spline(t, p);
% 预测未来位置
t_future = 2:0.1:3; % 预测时间 [2, 2.1, ..., 3]
p_pred = ppval(pp, t_future);
% IMU 外插（假设匀加速）
dt = t_future - t(end);
v_last = (ppval(pp, t(end)) - ppval(pp, t(end) - 0.01)) / 0.01; % 末端速度
a_last = accel(end); % 末端加速度
p_pred_imu = p(end) + v_last * dt + 0.5 * a_last * dt.^2;
% 融合样条和 IMU
p_pred_fused = 0.7 * p_pred + 0.3 * p_pred_imu;
% 可视化
figure;
plot(t, p, 'bo-', 'LineWidth', 1.5, 'DisplayName', '控制点');
hold on;
plot(t_future, p_pred, 'r--', 'LineWidth', 1.5, 'DisplayName', '样条预测');
plot(t_future, p_pred_imu, 'g--', 'LineWidth', 1.5, 'DisplayName', 'IMU 预测');
plot(t_future, p_pred_fused, 'm-', 'LineWidth', 1.5, 'DisplayName', '融合预测');
legend('show'); grid on;
xlabel('时间 (s)'); ylabel('位置 (m)');
title('一维位置预测示例');
hold off;
end