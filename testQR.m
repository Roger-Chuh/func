function testQR()
global use_self_made_givens use_self_made_householder
use_self_made_givens = false;
use_self_made_householder = true;



testA= [0     0    0   0  1 2 3
    0     0    0   0  4 5 6
    7     8    9 10  0 0 0
    10  11 12   0  0 0 0
    13 14 15    16 17 18 19
    20 21 22    23 24 25 26
    27 28 29    30 31 32 33];
testA = [testA [99 : (99 + size(testA,1)-1)]'];
% [Q, R] = QrGivens(testA);
[Q, R] = QrHouseholder(testA);
[Q1, R1] = qr(testA);

x = [3;4];

r = norm(x);

[Q, R] = qr(x);

cc = x(1) / r;
ss = x(2) / r;

G = [cc  ss;
    ss -cc];
x2 = G * [r;0];






Q1 = [cc;
    ss];
Q2 =  [ss;
    -cc];
r2 = G' * x;




Jpose = rand(2, 6);
Jidp = rand(2, 1);
res = rand(2);


J = [Jpose Jidp];

H = J' * J;
b = -J' * res;


r = norm(Jidp);
cc = Jidp(1) / r;
ss = Jidp(2) / r;
Q1 = [cc;
    ss];
Q2 =  [ss;
    -cc];

Q = [cc ss;
    -ss cc];

Jidp_check = Q1 * r

Q1' * Jidp_check




A1= [37 38 39;
    0 40 41;
    0 0 42;
    0 0 0;
    0 0 0;
    0 0 0;
    43 0 0;
    0 44 0;
    0 0 45];
A = zeros(9, 16);
A(:,13:15) = A1;
A(1:3, 1:6) = reshape([1:18], 6, 3)';
A(4:6, 7:12) = reshape([19:36], 6, 3)';
A(1:3, 16) = [46 47 48]';
if ~use_self_made_givens && ~use_self_made_householder
    [Q, R] = qr(A(:,13:15));
    Q = -Q;
    R = -R;
elseif use_self_made_givens
    [Q, ~, R, A_new] = QrGivens(A(:,13:15));
elseif use_self_made_householder
    [Q, R] = QrHouseholder(A(:,13:15));
end
Q11 = Q(:,1:3);
Q22 = Q(:,4:9);
[Q_idp, R_idp] = qr(Jidp);
Q_idp = -Q_idp;
R_idp = -R_idp;

check_landmark_r = Q' * A;
% check
% check_landmark_r(:,13:15) - A_new

Q11 * R(1:3,:)
Q * R

A_check = Q * check_landmark_r;


res_num = 5;
pose_num = 1;


[J_whole, res_whole] = GenJacs(res_num, pose_num);
end


function [J_whole, res_whole] = GenJacs(res_num, pose_num)
global use_self_made_givens use_self_made_householder
use_given_pose_num = true;
if pose_num < 0
    pose_num = res_num+1;
    use_given_pose_num = false;
end

J_whole = zeros(2 * res_num, pose_num * 6 + 1);
res_whole = zeros(2 * res_num,1);

dx_whole = rand(pose_num * 6 + 1,1);



for i = 1 : res_num
    
    J_whole(2 * i - 1 : 2 * i, 1:6) = rand(2, 6);
    if ~use_given_pose_num
        J_whole(2 * i - 1 : 2 * i, 6 * i + 1 : 6 * (i + 1)) = rand(2, 6);
    end
    J_whole(2 * i - 1 : 2 * i, pose_num * 6 + 1) = rand(2, 1);
    res_whole(2 * i - 1 : 2 * i, 1) = rand(2, 1);
end

res_whole = J_whole * dx_whole;

if ~use_self_made_givens && ~use_self_made_householder
    [Q, R] = qr(J_whole(:,pose_num * 6 + 1));
    Q = -Q;
    R = -R;
elseif use_self_made_givens
    [Q, ~, R, ~] = QrGivens(J_whole(:,pose_num * 6 + 1));
elseif use_self_made_householder
    [Q, R] = QrHouseholder(J_whole(:,pose_num * 6 + 1));
end

Q1 = Q(:,1);
Q2 = Q(:,2:2 * res_num);


Q1 * R(1) - J_whole(:,pose_num * 6 + 1)

%% 验证左零空间
left_null= Q2' * J_whole(:,pose_num * 6 + 1);

J_marged = Q' * J_whole;
res_marged = Q' * res_whole;


J_pose_only = Q2' * J_whole(:,1: pose_num * 6);
res_pose_only = Q2' * res_whole;

J_res_once_check = Q' * [J_whole res_whole];
J_res_once =[[ Q1' * J_whole(:,1: pose_num * 6);J_pose_only] R [Q1' * res_whole; res_pose_only]];
J_res_once2 = [Q' * J_whole Q' * res_whole];

H_pose_only = J_pose_only' * J_pose_only;
b_pose_only = J_pose_only' * res_pose_only;
dx_pose_only_solve = inv(H_pose_only) * b_pose_only;
dx_idp_only_solve = (1/R(1)) * Q1' * (res_whole - J_whole(:,1:pose_num * 6) * dx_pose_only_solve);

res_diff = res_pose_only -J_pose_only * dx_whole(1:6)

J_res_mat = [J_whole res_whole];
if ~use_self_made_givens && ~use_self_made_householder
    [Q_big, R_big] = qr(J_res_mat);
    Q_big = -Q_big;
    R_big = -R_big;
elseif use_self_made_givens
    [Q_big, ~, R_big, ~] = QrGivens(J_res_mat);
elseif use_self_made_householder
    [Q_big, R_big] = QrHouseholder(J_res_mat);
end

%使用这个trick的隐含条件是：marg后剩余的状态量一定要是排在后面的，在这个情形下，应该是把pose信息marg到idp上才对。
J_idp_QR_form = R_big(pose_num * 6+1:pose_num * 6+1, pose_num * 6+1:pose_num * 6+1);
res_idp_QR_form = R_big(pose_num * 6+1:pose_num * 6+1, pose_num * 6 + 2);
H_idp_QR_form = J_idp_QR_form' * J_idp_QR_form;
b_idp_QR_form = J_idp_QR_form' * res_idp_QR_form;
dx_idp_QR_form = inv(H_idp_QR_form) * b_idp_QR_form;
dx_pose_QR_form = (inv(R_big(1:pose_num * 6, 1:pose_num * 6))) * (R_big(1:pose_num*6, pose_num*6+2) -R_big(1:pose_num*6,pose_num*6+1) * dx_idp_QR_form);

dx_pose_idp_QR_form2 = (inv(R_big(1:pose_num * 6+1,1:pose_num * 6+1)' * R_big(1:pose_num * 6+1, 1:pose_num * 6+1))) * R_big(1:pose_num * 6+1,1:pose_num * 6+1)'* R_big(1:pose_num * 6+1,pose_num * 6+2);

err_dx = [dx_pose_QR_form; dx_idp_QR_form] - dx_pose_idp_QR_form2;

res_marged_check = J_marged * dx_whole;

H_whole = J_whole' * J_whole;
rank(H_whole)
b_whole = J_whole' * res_whole;
dx_whole_solve = inv(H_whole) * b_whole;

err_dx2 = dx_whole_solve - dx_pose_idp_QR_form2;

H11 = H_whole(1:pose_num * 6, 1:pose_num * 6) - H_whole(1:pose_num * 6, pose_num * 6 + 1) * (1/H_whole(pose_num * 6+1,pose_num * 6+1)) * H_whole(1:pose_num * 6, pose_num * 6 + 1)';

b1 = b_whole(1:pose_num * 6) - H_whole(1:pose_num * 6, pose_num * 6 + 1) * (1/H_whole(pose_num * 6+1,pose_num * 6+1)) * b_whole(pose_num * 6 + 1);
H11_lazy = H_whole(1:pose_num * 6, 1:pose_num * 6) - (Q1' * J_whole(:,1: pose_num * 6))' * (Q1' * J_whole(:,1: pose_num * 6));
b1_lazy = b_whole(1:pose_num * 6) - (Q1' * J_whole(:,1: pose_num * 6))' * (Q1' * res_whole);
%% 重要验证项
H11_diff = H11_lazy - H11;
b1_diff = b1 - b1_lazy;



dx_sc_pose = inv(H11) * b1;
dx_sc_idp = (1/H_whole(pose_num * 6+1,pose_num * 6+1)) * ( b_whole(pose_num * 6 + 1) - H_whole(1:pose_num * 6, pose_num * 6 + 1)' * dx_sc_pose);

dx_pose_diff = dx_sc_pose - dx_pose_only_solve

dx_diff = dx_whole_solve - dx_whole

method = 1;
J_res_accum = DoQrIncrementally(J_whole, res_whole, use_given_pose_num, res_num, pose_num, method)
qr_diff = J_res_accum - J_res_once
qr_diff2 = J_res_once2 - J_res_once
qr_diff3 = J_res_once_check - J_res_once

H11_diff = J_res_once(2:end,1:pose_num * 6)' * J_res_once(2:end,1:pose_num * 6) - J_res_accum(2:end,1:pose_num * 6)' * J_res_accum(2:end,1:pose_num * 6)
b1_diff = J_res_once(2:end,1:pose_num * 6)' * J_res_once(2:end,pose_num * 6 + 2) - J_res_accum(2:end,1:pose_num * 6)' * J_res_accum(2:end,pose_num * 6 + 2)
end
function J_res_accum = DoQrIncrementally(J_whole, res_whole, use_given_pose_num, res_num, pose_num, method)
global use_self_made_givens use_self_made_householder
if use_given_pose_num
    J_res_accum = [];zeros(1,pose_num * 6 + 1);
    
end
for i = 1 : res_num
    J_each = J_whole(2 * i-1:2*i,:);
    res_each = res_whole(2 * i-1:2*i);
    if ~use_self_made_givens && ~use_self_made_householder
        [Q_each, R_each] = qr(J_each(:,pose_num * 6 + 1));
        Q_each = -Q_each;
        R_each = -R_each;
    elseif use_self_made_givens
        [Q_each,~, R_each, ~] = QrGivens(J_each(:,pose_num * 6 + 1));
    elseif use_self_made_householder
        [Q_each, R_each] = QrHouseholder(J_each(:,pose_num * 6 + 1));
    end
    Q1_each = Q_each(:,1);
    Q2_each = Q_each(:,2:2);
    
    Q1t_Jpose_each = Q1_each' * J_each(:,1:pose_num * 6);
    Q1t_res_each = Q1_each' * res_each;
    Q2t_Jpose_each = Q2_each' * J_each(:,1:pose_num * 6);
    Q2t_res_each = Q2_each' * res_each;
    
    if method == 0
        if 1%size(J_res_accum,1) <= 2
            J_res_accum = [J_res_accum; [[Q1t_Jpose_each;Q2t_Jpose_each] R_each] [Q1t_res_each;Q2t_res_each]];
        end
        if size(J_res_accum,1) > 2
            if ~use_self_made_givens && ~use_self_made_householder
                [Q_accum, R_accum] = qr(J_res_accum(:,pose_num * 6 + 1));
                Q_accum = -Q_accum;
                R_accum = -R_accum;
            elseif use_self_made_givens
                [Q_accum,~, R_accum] = QrGivens(J_res_accum(:,pose_num * 6 + 1));
            elseif use_self_made_householder
                [Q_accum, R_accum] = QrHouseholder(J_res_accum(:,pose_num * 6 + 1));
            end
            Q1_accum = Q_accum(:,1);
            Q2_accum = Q_accum(:,2:size(Q_accum, 2));
            Q1_accum * R_accum(1)
            % 验证Q2左零空间
            Q2_accum' * J_res_accum(:,pose_num * 6 + 1);
            Q1t_Jpose_accum = Q1_accum' * J_res_accum(:,1:pose_num * 6);
            Q1t_res_accum = Q1_accum' * J_res_accum(:,pose_num * 6 + 2);
            Q2t_Jpose_accum = Q2_accum' * J_res_accum(:,1:pose_num * 6);
            Q2t_res_accum = Q2_accum' * J_res_accum(:,pose_num * 6 + 2);
            J_res_accum = [[[Q1t_Jpose_accum;Q2t_Jpose_accum] R_accum] [Q1t_res_accum;Q2t_res_accum]];
        end
    elseif 1
        if size(J_res_accum,1) < 2
            J_res_accum = [J_res_accum; [[Q1t_Jpose_each;Q2t_Jpose_each] R_each] [Q1t_res_each;Q2t_res_each]];
        else
            J_res_comb = [J_res_accum(1,:);[Q1t_Jpose_each R_each(1) Q1t_res_each]];
            if ~use_self_made_givens && ~use_self_made_householder
                [Q_comb, R_comb] = qr(J_res_comb(:,pose_num * 6 + 1));
                Q_comb = -Q_comb;
                R_comb = -R_comb;
            elseif use_self_made_givens
                [Q_comb, ~,R_comb] = QrGivens(J_res_comb(:,pose_num * 6 + 1));
            elseif use_self_made_householder
                [Q_comb, R_comb] = QrHouseholder(J_res_comb(:,pose_num * 6 + 1));
            end
            Q1_comb = Q_comb(:,1);
            Q2_comb = Q_comb(:,2:size(Q_comb, 2));
            Q1_comb * R_comb(1)
            % 验证Q2左零空间
            Q2_comb' * J_res_comb(:,pose_num * 6 + 1);
            Q1t_Jpose_comb = Q1_comb' * J_res_comb(:,1:pose_num * 6);
            Q1t_res_comb = Q1_comb' * J_res_comb(:,pose_num * 6 + 2);
            Q2t_Jpose_comb = Q2_comb' * J_res_comb(:,1:pose_num * 6);
            Q2t_res_comb = Q2_comb' * J_res_comb(:,pose_num * 6 + 2);
            
            J_res_accum(1,:) = [Q1t_Jpose_comb R_comb(1) Q1t_res_comb];
            
            J_res_accum = [J_res_accum; [Q2t_Jpose_comb 0 Q2t_res_comb]];
            J_res_accum = [J_res_accum; [Q2t_Jpose_each 0 Q2t_res_each]];
            
            %             J_res_accum = [[[Q1t_Jpose_comb;Q2t_Jpose_comb] R_comb] [Q1t_res_comb;Q2t_res_comb]];
        end
        
    end
end

end