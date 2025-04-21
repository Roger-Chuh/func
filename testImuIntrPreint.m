function testImuIntrPreint
% prediction filter
A = zeros(18, 18);
w = rand(3,1);
vb = rand(3,1);
A(1:3,1:3) = -SkewSymMat(w);
A(1:3,4:6) = -eye(3);
A(4:6, 7:9) = eye(3);
A(10:12, 1:3) = -SkewSymMat(vb);
A(10:12, 10:12) = -SkewSymMat(w);
A(10:12, 13:15) = eye(3);
A(13:15, 16:18) = eye(3);

dt = 0.001;
F = expm(A * dt);

F_check = eye(18) + A * dt;


Jl(w * dt);


%% new eskf from ce, use left perturb
A = zeros(18, 18);
R = rodrigues(rand(3,1));
w = rand(3,1);
A(1:3,4:6) = R;
A(4:6, 7:9) = eye(3);
A(10:12,10:12) = SkewSymMat(R * w);
A(10:12,13:15) = eye(3);
A(13:15,13:15) = SkewSymMat(R * w);
A(13:15,16:18) = eye(3);
A(16:18,16:18) = SkewSymMat(R * w);
dt = 0.01;
F = expm(A * dt);

F_check = eye(18) + A * dt;

rodrigues(R * w * dt) * dt - F(10:12,13:15)
rodrigues(R * w * dt) * 0.5 * dt * dt - F(10:12,16:18)


%% new new eskf from ce, use left perturb SE3EKF LeftMulti
A = zeros(18, 18);
R = rodrigues(rand(3,1));
w = rand(3,1);
p = rand(3,1);
v = rand(3,1);
A(1:3,4:6) = R;
A(4:6, 7:9) = eye(3);
A(10:12,4:6) = SkewSymMat(p) * R;
A(10:12,13:15) = eye(3);
A(13:15,4:6) = SkewSymMat(v) * R;
A(13:15,16:18) = eye(3);
A(16:18,16:18) = SkewSymMat(R * w);

dt = 0.001;
F = expm(A * dt);

F_check = eye(18) + A * dt;
F1 = R * dt;
F2 = 0.5 * dt^2 * R;
F3 = dt * eye(3);
F4 = SkewSymMat(p) * R * dt + 0.5 * SkewSymMat(v) * R * dt^2;
F5 = 0.5 * SkewSymMat(p) * R * dt^2 + 1/6 * SkewSymMat(v) * R * dt^3;
F6 = SkewSymMat(v) * R * dt;
F7 = 0.5 * SkewSymMat(v) * R * dt^2;
F8 = eye(3) * dt;
F9 = IntegrateJlTau(R * w, dt);
F10 = Jl(R * w * dt) * dt;
F11 = rodrigues(R * w * dt);

F1 - F(1:3,4:6)
F2 - F(1:3,7:9)
F3 - F(4:6,7:9)
F4 - F(10:12,4:6)
F5 - F(10:12,7:9)
F6 - F(13:15,4:6)
F7 - F(13:15,7:9)
F8 - F(10:12,13:15)
F9 - F(10:12,16:18)
F10 - F(13:15,16:18)
F11 - F(16:18,16:18)

%% new new eskf from ce, use left perturb SE3EKF LeftMulti2
A = zeros(18, 18);
R = rodrigues(rand(3,1));
w = rand(3,1);
p = rand(3,1);
v = rand(3,1);
A(1:3,4:6) = R;
A(4:6, 7:9) = eye(3);
A(10:12,4:6) = SkewSymMat(p) * R;
A(10:12,13:15) = eye(3);

A(13:15,13:15) = SkewSymMat(R * w);
A(13:15,16:18) = R;

% A(16:18,16:18) = SkewSymMat(R * w);

dt = 0.001;
F = expm(A * dt);

F_check = eye(18) + A * dt;

F1 = R * dt;
F2 = 0.5 * dt^2 * R;
F3 = dt * eye(3);
F4 = SkewSymMat(p) * R * dt;
F5 = 0.5 * SkewSymMat(p) * R * dt^2;
F6 = Jl(R * w * dt) * dt;
F7 = IntegrateJlTau(R * w, dt) * R;
F7_check = Jl(R * w * dt) * dt * Jl(-R * w * dt) * dt * R + rodrigues(R * w * dt) * IntegrateJlTau(-R * w, dt) * R;

rodrigues(R * w * dt) * Jl(-R * w * dt) - Jl(R * w * dt)
rodrigues(w * dt) * Jl(-w * dt) * rodrigues(-w * dt) - Jl(w * dt)

F8 = rodrigues(R * w * dt);
F9 = Jl(R * w * dt) * dt * R;

F1 - F(1:3,4:6)
F2 - F(1:3,7:9)
F3 - F(4:6, 7:9)
F4 - F(10:12, 4:6)
F5 - F(10:12, 7:9)
F6 - F(10:12, 13:15)
F7 - F(10:12, 16:18)
F8 - F(13:15, 13:15)
F9 - F(13:15, 16:18)

%% 20250414/1820
A = zeros(18, 18);
R = rodrigues(rand(3,1));
w = rand(3,1);
a = rand(3,1);
A(1:3,1:3) = -SkewSymMat(w);
A(1:3,4:6) = eye(3);
A(4:6,7:9) = eye(3);
A(10:12,10:12) = -SkewSymMat(w);
A(10:12,13:15) = eye(3);
A(13:15,13:15) = -SkewSymMat(w);
A(13:15,16:18) = eye(3);
A(16:18, 1:3) = SkewSymMat(a) * R * SkewSymMat(w);
A(16:18, 4:6) = -SkewSymMat(a) * R;
A(16:18, 16:18) = SkewSymMat(R * w);

dt = 0.001;
F = expm(A * dt);

F_check = eye(18) + A * dt;

%% 20250415/1514
A = zeros(18, 18);
R = rodrigues(rand(3,1));
dt = 0.01;
p = rand(3,1);
A(1:3,4:6) = R;
A(1:3,7:9) = 0.5 * dt * R;
A(4:6, 7:9) = eye(3);
A(10:12,13:15) = eye(3);
A(16:18, 4:6) = SkewSymMat(p) * R;
A(16:18, 7:9) = SkewSymMat(p) * 0.5 * dt * R;
A(16:18, 10:12) = R;
F = expm(A * dt);

F_check = eye(18) + A * dt;


%% 20250415/1550
A = zeros(18, 18);
R = rodrigues(rand(3,1));
dt = 0.01;
v = rand(3,1);
A(1:3,4:6) = R;
A(1:3,7:9) = 0.5 * dt * R;
A(4:6, 7:9) = eye(3);
A(10:12,13:15) = eye(3);
A(16:18, 1:3) = -SkewSymMat(v);
A(16:18, 10:12) = R;
F = expm(A * dt);


%% 20250415/1630
A = zeros(18, 18);
R = rodrigues(rand(3,1));
dt = 0.01;
v = rand(3,1);
w = rand(3,1);
A(1:3,1:3) = -SkewSymMat(w);
A(1:3,4:6) = eye(3);
A(1:3,7:9) = 0.5 * dt * eye(3);
A(4:6,7:9) = eye(3);
A(10:12,13:15) = eye(3);
A(16:18,1:3) = R * SkewSymMat(v);
A(16:18,10:12) = R;
F = expm(A * dt);


%% 20250415/2038
A = zeros(18, 18);
w = rand(3,1);
dt = 0.01;
v = rand(3,1);
p = rand(3,1);
A(1:3,1:3) = -SkewSymMat(w);
A(1:3,4:6) = eye(3);
A(4:6, 7:9) = eye(3);
A(10:12, 4:6) = SkewSymMat(p);
A(10:12, 10:12) = -SkewSymMat(w);
A(10:12, 13:15) = eye(3);
A(13:15, 4:6) = SkewSymMat(v);
A(13:15, 13:15) = -SkewSymMat(w);
A(13:15, 16:18) = eye(3);
F = expm(A * dt);



%% 20250415/2102
A = zeros(18, 18);
R = rodrigues(rand(3, 1));
w = rand(3,1);
dt = 0.01;
v = rand(3,1);
p = rand(3,1);

A(1:3,4:6) = R;
A(4:6, 7:9) = eye(3);
A(10:12, 4:6) = SkewSymMat(p);
A(10:12, 10:12) = -SkewSymMat(w);
A(10:12, 13:15) = eye(3);
A(13:15, 4:6) = SkewSymMat(v);
A(13:15, 13:15) = -SkewSymMat(w);
A(13:15, 16:18) = eye(3);
F = expm(A * dt);


%% 20251110
A = zeros(18, 18);
R = rodrigues(rand(3, 1));
w = rand(3,1);
dt = 0.01;
v = rand(3,1);
p = rand(3,1);

A(1:3,1:3) = -SkewSymMat(w);
A(1:3,4:6) = eye(3);
A(4:6, 7:9) = eye(3);
A(10:12, 4:6) = R * SkewSymMat(p);
A(10:12, 13:15) = eye(3);
A(13:15, 4:6) = R * SkewSymMat(v);
A(13:15, 16:18) = eye(3);
A(16:18, 16:18) = SkewSymMat(R * w);
F = expm(A * dt);

%% 20251116
A = zeros(18, 18);
R = rodrigues(rand(3, 1));
w = rand(3,1);
dt = 0.01;
v = rand(3,1);
p = rand(3,1);

A(1:3,4:6) = R;
A(4:6, 7:9) = eye(3);
A(10:12, 4:6) = R * SkewSymMat(p);
A(10:12, 13:15) = eye(3);
A(13:15, 4:6) = R * SkewSymMat(v);
A(13:15, 16:18) = eye(3);
A(16:18, 16:18) = SkewSymMat(R * w);
F = expm(A * dt);


F12 = R * dt;
F13 = 0.5 * dt^2 * R;
F23 = dt * eye(3);
F42 = R * SkewSymMat(p) * dt  + 0.5 * R * SkewSymMat(v) * dt^2;% (dt + 0.5 * dt^2);
F43 = 0.5 * R * SkewSymMat(p) * dt^2 + 1/6 * R * SkewSymMat(v) * dt^3;
F45 = dt * eye(3);
F46 = IntegrateJlTau(R * w, dt);
F52 = R * SkewSymMat(v) * dt;
F56 = Jl(R * w * dt) * dt;

F66 = rodrigues(R * w * dt);

F12 - F(1:3, 4:6)
F13 - F(1:3, 7:9)
F23 - F(4:6, 7:9)
F42 - F(10:12, 4:6)
F43 - F(10:12, 7:9)
F45 - F(10:12, 13:15)
F46 - F(10:12, 16:18)
F52 - F(13:15, 4:6)
F56 - F(13:15, 16:18)
F66 - F(16:18, 16:18)


%% 20250417/1514
A = zeros(18, 18);
R = rodrigues(rand(3, 1));
w = rand(3,1);
dt = 0.01;
a = rand(3,1);

A(1:3,4:6) = R';
A(4:6,7:9) = eye(3);
A(10:12, 10:12) = SkewSymMat(-R' * w);
A(10:12, 13:15) = eye(3);
A(13:15, 13:15) = SkewSymMat(-R' * w);
A(13:15, 16:18) = eye(3);
A(16:18, 4:6) = -R' * SkewSymMat(a);

F = expm(A * dt);




%% 20250421/1149
A = zeros(18, 18);
R = rodrigues(rand(3, 1));
w = rand(3,1);
dt = 0.01;
A(1:3,4:6) = R';
A(4:6,7:9) = eye(3);
A(10:12,10:12) = SkewSymMat(-R' * w);
A(10:12,13:15) = eye(3);
A(13:15,13:15) = SkewSymMat(-R' * w);
A(13:15,16:18) = eye(3);
A(16:18,16:18) = SkewSymMat(-R' * w);
F = expm(A * dt);

%% single cam relative imu
A = zeros(15,15);
A2 = zeros(15,15);
w = rand(3, 1);
a = rand(3, 1);
dt = 0.01;
A(1:3,1:3) = -SkewSymMat(w);
A(1:3,10:12) = eye(3);
A(4:6, 1:3) = -SkewSymMat(a);
A(4:6, 4:6) = -SkewSymMat(w);
A(4:6, 13:15) = eye(3);
A(7:9, 4:6) = eye(3);
A(7:9, 7:9) = -SkewSymMat(w);
F = expm(A * dt);

A2(1:3,1:3) = -SkewSymMat(w);
A2(1:3,10:12) = -eye(3);
A2(4:6, 1:3) = -SkewSymMat(a);
A2(4:6, 4:6) = -SkewSymMat(w);
A2(4:6, 13:15) = -eye(3);
A2(7:9, 4:6) = eye(3);
A2(7:9, 7:9) = -SkewSymMat(w);
F2 = expm(A2 * dt);



A11 = rodrigues(w * dt)';
A22 = rodrigues(w * dt)';
A33 = rodrigues(w * dt)';
A21 = -rodrigues(w * dt)' * SkewSymMat(a * dt);

%%%%%
X = rodrigues(-w*dt) * (-dt^2 * Jl(w * dt) + IntegrateJlTau(w, dt));
X - F(1:3,7:9)
p1 = -SkewSymMat(w) * Jl(-w * dt) * dt
p1 - F(10:12, 1:3)
p2 = -rodrigues(-w * dt) * SkewSymMat(Jl(w * dt) * dt * vb);
p2 - F(10:12, 4:6)

F12 = -rodrigues(-w * dt) * Jl(w * dt) * dt;
F12 - F(1:3,4:6)
F12 + F(10:12,13:15)
X + F(10:12,16:18)

% preint in orca next
a = rand(3,1);
w = rand(3,1);
F = zeros(30, 30);
m_vec = 0.02 * rand(3, 1);
a_vec = 0.02 * rand(6,1);
d_vec = 0.02 * rand(6,1);
A = zeros(3,3);
D = zeros(3,3);


A(0 + 1, 0 + 1) = a_vec(0 + 1);
A(0 + 1, 1 + 1) = a_vec(1 + 1);
A(0 + 1, 2 + 1) = a_vec(2 + 1);
A(1 + 1, 1 + 1) = a_vec(3 + 1);
A(1 + 1, 2 + 1) = a_vec(4 + 1);
A(2 + 1, 2 + 1) = a_vec(5 + 1);
A = A + eye(3);

D(0 + 1, 0 + 1) = d_vec(0 + 1);
D(0 + 1, 1 + 1) = d_vec(1 + 1);
D(0 + 1, 2 + 1) = d_vec(2 + 1);
D(1 + 1, 1 + 1) = d_vec(3 + 1);
D(1 + 1, 2 + 1) = d_vec(4 + 1);
D(2 + 1, 2 + 1) = d_vec(5 + 1);
D = D + eye(3);

M = rodrigues(m_vec);

F(1:3, 1:3) = -SkewSymMat(M * A * w);
F(1:3, 10:12) = -M * A;
F(1:3, 16:18) = -SkewSymMat(M * A * w);
F(1:3, 19:24) = M * func_K(w);

F(4:6, 1:3) = -SkewSymMat(D * a);
F(4:6, 4:6) = -SkewSymMat(M * A * w);
F(4:6, 13:15) = -D;
F(4:6, 25:30) = func_K(a);

F(7:9, 4:6) = eye(3);
F(7:9, 7:9) = -SkewSymMat(M * A * w);

dt = 0.01;
AA = expm(F * dt);


leftJmat = Jl(M * A * w * dt);
dRmat = rodrigues(M * A * w * dt);
dRmat_T = dRmat';

A11 = dRmat_T;
err_A11 = dRmat_T - AA(1:3,1:3)
A14 = dRmat_T * (-leftJmat * M * A * dt);
err_A14 = A14 - AA(1:3, 10:12)
A16 = dRmat_T * (-leftJmat * dt * SkewSymMat(M * A * w));
err_A16 = A16 - AA(1:3, 16:18)
A17 = dRmat_T * (leftJmat * dt * M * func_K(w));
err_A17 = A17 - AA(1:3, 19:24)

A21 = -dRmat_T * SkewSymMat(leftJmat * dt * D * a);
err_A21 = A21 - AA(4:6,1:3)
A22 = dRmat_T;
err_A22 = A22 - AA(4:6,4:6)
A24_2 = dRmat_T * (-leftJmat * dt * SkewSymMat(D * a) * (-leftJmat * dt * M * A));
A24 = 0.5 * dt * dt * SkewSymMat(D * a) - (SkewSymMat(D * a) * SkewSymMat(M * A * w) + (SkewSymMat(D * a) * SkewSymMat(M * A * w))') * (dt * dt * dt) / 6.0;
err_A24 = A24 - AA(4:6,10:12)
A25 = dRmat_T * (-leftJmat * D * dt);
err_A25 = A25 - AA(4:6,13:15)
A26 = dRmat_T * (-leftJmat * dt * SkewSymMat(D * a));
err_A26 = A26 - AA(4:6,16:18)

A31 = -dRmat_T * SkewSymMat(IntegrateJlt(M * A * w, dt) * D * a);
err_A31 = A31 - AA(7:9,1:3)
A32 = dRmat_T * dt;
err_A32 = A32 - AA(7:9,4:6)
A33 = dRmat_T;
err_A33 = A33 - AA(7:9,7:9)
A34 = dRmat_T * SkewSymMat(leftJmat * D * a * dt * dt * dt / 6.0);
err_A34 = A34 - AA(7:9,10:12)
A35 = SkewSymMat(M * A * w * dt * dt * dt / 3.0) - (dt * dt / 2.0) * eye(3);
err_A35 = A35 - AA(7:9,13:15)
A36 = dRmat_T * (leftJmat * dt);
err_A36 = A36 - AA(7:9,16:18)




%% larvio(qiu xiao chen)'s formulation
FF = zeros(15, 15);
FF(1:3, 10:12) = rodrigues(w);
FF(4:6, 1:3) = -SkewSymMat(rodrigues(w) * a);
FF(4:6, 13:15) = -rodrigues(w);
FF(7:9, 4:6) = eye(3);
AAA = expm(FF * dt);

end