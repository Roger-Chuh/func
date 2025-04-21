function testPreintSymbol()
close all
% syms t real  % 声明符号变量t（时间步长）
%
% % 定义符号向量
% w = sym('w', [3 1], 'real');
% a = sym('a', [3 1], 'real');
%
% % 构造符号斜对称矩阵
% skew_sym = @(v) [
%     0    -v(3)  v(2)
%     v(3)  0    -v(1)
%    -v(2)  v(1)  0
% ];
%
% % 构建符号矩阵A
% A = [
%     -skew_sym(w)  zeros(3,6)
%     -skew_sym(a) -skew_sym(w) zeros(3,3)
%     zeros(3,3)    eye(3)      -skew_sym(w)
% ]
%
% % 尝试计算矩阵指数（可能需要极长时间）
% expA = expm(A*t)
%
syms wx wy wz ax ay az dt real

w = [wx; wy; wz];
a = [ax; ay; az];


cA = [
    -skew(w), zeros(3, 6), -eye(3), zeros(3, 3);
    -skew(a), -skew(w), zeros(3, 6), -eye(3);
    zeros(3, 3), eye(3), -skew(w), zeros(3, 6);
    zeros(6, 15)
    ];

% 计算 cA^2 和 cA^3
cA2 = cA * cA;
cA3 = cA * cA2;
cA4 = cA2 * cA2;

% 计算 discret_AA  把exp(F * dt)离散化
discret_AA = eye(15) + cA * dt + (cA2 * dt^2)/2 + (cA3 * dt^3)/6% + (cA4 * dt^4)/24

if 0
    exp(sin(1) * 0.01) - (1 + sin(1) * dt + sin(1) * sin(1) * dt^2/2 + sin(1) * sin(1) * sin(1) * dt^3/6 + sin(1) * sin(1) * sin(1) * sin(1)* dt^4/24 + (sin(1)* dt)^5/120)
    exp(cos(1) * dt) - (1 + cos(1) * dt + cos(1)^2 * dt^2/2 + cos(1)^3 * dt^3/6 + cos(1)^4* dt^4/24 + (cos(1)* dt)^5/120)
    
    
    a = 0.01;
    x = 0.2;
    orig = cos(x);
    check = cos(a) - sin(a) * (x-a) - cos(a) * (x-a)^2/(2*1) + sin(a) * (x-a)^3/(2*3) + cos(a) * (x-a)^4/(2*3*4)-sin(a) * (x-a)^5/(2*3*4*5)+cos(a) * (x-a)^6/(2*3*4*5*6)-cos(a) * (x-a)^7/(2*3*4*5*6*7);
    diff = orig - check
    
end


simplified_discret_AA_J5 = simplify(discret_AA(4:6,10:12))
part_J5 = 0.5*dt^2*skew(a) - (skew(a)*skew(w)+skew(w)*skew(a))*dt^3/6

% simplified_discret_AA_J6 = simplify(discret_AA(7:9,10:12))
% skew(leftJmat * acc_minus_ba * dt3 / 6.0)
% part_J6 = rodrigues(w * dt)


simplified_discret_AA_J7 = simplify(discret_AA(7:9,13:15))
part_J7 = skew(w * dt * dt * dt / 3.0) - (dt * dt / 2.0) * eye(3)



if 1
    
    A = -2;
    B = 1;
    dt = 0.1;
    
    % 离散化
    Ad = expm(A * dt); % 矩阵指数
    Bd = (expm(A * dt) - 1) / A * B; % 解析积分
    
    disp(['Ad = ', num2str(Ad)]);
    disp(['Bd = ', num2str(Bd)]);
    
    % 模拟
    x = zeros(1, 50);
    x(1) = 1; % 初始值
    u = ones(1, 50); % 恒定输入
    for k = 1:49
        x(k+1) = Ad * x(k) + Bd * u(k);
    end
    
    % 绘图
    t = 0:dt:dt*49;
    figure, plot(t, x, 'b-', 'LineWidth', 2);
    xlabel('t');
    ylabel('x(t)');
    title('Discrete System Response');
    grid on;
    hold on;plot(t, exp(A * t) / u, 'or')
end

end