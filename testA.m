function testA()

aw1 = rand(3,1);

a1 = rand(3,1);
a2 = rand(3,1);
w1 = rand(3,1);
w2 = rand(3,1);

R12 = rodrigues(rand(3,1));
p12 = rand(3,1);
v12 = rand(3,1);



dt = 0.01;

A_w_jerk = makeA(true, aw1,a1,a2,w1,w2,R12,p12,v12);

F_w_jerk = expm(A_w_jerk * dt);
F_w_jerk_mod = F_w_jerk;
F_w_jerk_mod(22:33,1:33) = zeros(12,33);
F_w_jerk_mod(1:33,22:33) = zeros(12,33)';

A_wo_jerk = makeA(false, aw1,a1,a2,w1,w2,R12,p12,v12);
F_wo_jerk = zeros(33,33);
F_wo_jerk(1:21,1:21) = expm(A_wo_jerk(1:21,1:21) * dt);
F_wo_jerk_33 = expm(A_wo_jerk * dt);

diff = F_wo_jerk - F_w_jerk;
diff_mod = F_wo_jerk - F_w_jerk_mod;

max(abs(diff_mod(:)))

end
function A = makeA(use_jerk, aw1,a1,a2,w1,w2,R12,p12,v12)
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
A(state_vel_start_index+1:state_vel_start_index+3, state_rot_start_index+1:state_rot_start_index+3) = -R12 * SkewSymMat(a2);
A(state_vel_start_index+1:state_vel_start_index+3, state_ba1_start_index+1:state_ba1_start_index+3) = -eye(3);

A(state_vel_start_index+1:state_vel_start_index+3, state_bg1_start_index+1:state_bg1_start_index+3) = SkewSymMat(w1) * SkewSymMat(p12) + SkewSymMat(SkewSymMat(w1) * p12) + 2 * SkewSymMat(v12);

A(state_vel_start_index+1:state_vel_start_index+3, state_ba2_start_index+1:state_ba2_start_index+3) = R12;

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
    
    A(state_ba1_start_index+1:state_ba1_start_index+3, state_aa1_start_index+1:state_aa1_start_index+3) = eye(3);
    
    A(state_ba2_start_index+1:state_ba2_start_index+3, state_aa2_start_index+1:state_aa2_start_index+3) =eye(3);
end
end