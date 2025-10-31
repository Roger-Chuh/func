function testErrorStateTransitionMatrix()
aw1 = rand(3,1);

a1 = zeros(3, 1);rand(3,1);
a2 = zeros(3, 1);rand(3,1);
w1 = rand(3,1);
w2 = rand(3,1);

R12 = rodrigues(rand(3,1));
p12 = rand(3,1);
v12 = rand(3,1);


dt = 0.01;

A_w_a = makeA(true, true, aw1,a1,a2,w1,w2,R12,p12,v12);

F_w_a = expm(A_w_a * dt);
F_w_a_mod = F_w_a;
F_w_a_mod([10:12 16:18 28:33],1:33) = zeros(12,33);
F_w_a_mod(1:33,[10:12 16:18 28:33]) = zeros(12,33)';
F_w_a_mod(10:12, 10:12) = eye(3);
F_w_a_mod(16:18, 16:18) = eye(3);
F_w_a_mod(28:33, 28:33) = eye(6);

A_wo_a = makeA(true, false, aw1,a1,a2,w1,w2,R12,p12,v12);
% F_wo_a = zeros(33,33);
% F_wo_a(1:21,1:21) = expm(A_wo_a(1:21,1:21) * dt);
F_wo_a_33 = expm(A_wo_a * dt);

% diff = F_wo_a - F_w_a;
diff_mod = F_wo_a_33 - F_w_a_mod;

max(abs(diff_mod(:)))



end
function A = makeA(use_jerk, use_a, aw1,a1,a2,w1,w2,R12,p12,v12)
state_pose_start_index = 0;
state_vel_start_index = 3;
state_rot_start_index = 6;

state_ba1_start_index = 9;
state_bg1_start_index = 12;
state_ba2_start_index = 15;
state_bg2_start_index = 18;

state_aw1_start_index = 21;
state_aw2_start_index = 24;
state_aa1_start_index = 27;
state_aa2_start_index = 30;




A = zeros(33,33);
% // state err p
A(state_pose_start_index+1:state_pose_start_index+3, state_vel_start_index+1:state_vel_start_index+3) = eye(3);

%     // state err v
A(state_vel_start_index+1:state_vel_start_index+3, state_pose_start_index+1:state_pose_start_index+3) = -SkewSymMat(w1) * SkewSymMat(w1) - SkewSymMat(aw1);
A(state_vel_start_index+1:state_vel_start_index+3, state_vel_start_index+1:state_vel_start_index+3) = -2 * SkewSymMat(w1);
if 1 %use_a
    A(state_vel_start_index+1:state_vel_start_index+3, state_rot_start_index+1:state_rot_start_index+3) = -R12 * SkewSymMat(a2);
end
if use_a
    A(state_vel_start_index+1:state_vel_start_index+3, state_ba1_start_index+1:state_ba1_start_index+3) = -eye(3);
end
A(state_vel_start_index+1:state_vel_start_index+3, state_bg1_start_index+1:state_bg1_start_index+3) = SkewSymMat(w1) * SkewSymMat(p12) + SkewSymMat(SkewSymMat(w1) * p12) + 2 * SkewSymMat(v12);
if use_a
    A(state_vel_start_index+1:state_vel_start_index+3, state_ba2_start_index+1:state_ba2_start_index+3) = R12;
end
if use_jerk
    A(state_vel_start_index+1:state_vel_start_index+3, state_aw1_start_index+1:state_aw1_start_index+3) = SkewSymMat(p12);
end
%     // state err R
A(state_rot_start_index+1:state_rot_start_index+3, state_rot_start_index+1:state_rot_start_index+3) = -SkewSymMat(w2);
A(state_rot_start_index+1:state_rot_start_index+3, state_bg1_start_index+1:state_bg1_start_index+3) = -R12';

A(state_rot_start_index+1:state_rot_start_index+3, state_bg2_start_index+1:state_bg2_start_index+3) = eye(3);


if use_jerk
    A(state_bg1_start_index+1:state_bg1_start_index+3, state_aw1_start_index+1 : state_aw1_start_index+3) = eye(3);
    
    A(state_bg2_start_index+1:state_bg2_start_index+3, state_aw2_start_index+1:state_aw2_start_index+3) =eye(3);
    if use_a
        A(state_ba1_start_index+1:state_ba1_start_index+3, state_aa1_start_index+1:state_aa1_start_index+3) = eye(3);
    end
    if use_a
        A(state_ba2_start_index+1:state_ba2_start_index+3, state_aa2_start_index+1:state_aa2_start_index+3) =eye(3);
    end
end
end