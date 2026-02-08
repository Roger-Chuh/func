function testMergePreint()
global dt

dt = 0.003;
acc1 = rand(3, 100);
acc2 = rand(3, 100);
gyro1 = rand(3, 50);
gyro2 = rand(3, 50);
acc3 = [acc1 acc2];
gyro3 = [gyro1 gyro2];

[J1_1, J2_1, J3_1] = Preint(gyro1, acc1);
[J1_2, J2_2, J3_2] = Preint(gyro2, acc2);
[J1_3, J2_3, J3_3] = Preint(gyro3, acc3);


J2_3_check = J2_1 + J1_1 * J2_2;
J2_diff = J2_3_check - J2_3;

J3_3_check = J3_1 + J2_1 * ((size(gyro2,2) + 0) * dt) + J1_1 * (J3_2 + J2_2 * ((size(gyro2,2) + 0) * dt));
J3_diff = J3_3_check - J3_3;

J1_3_check = J1_1 * J1_2;
J1_diff = J1_3_check - J1_3;
end

function [J1, J2, J3] = Preint(gyro, acc)
global dt

J1 = eye(3);
J2 = zeros(3, 1);
J3 = zeros(3, 1);
for i = 1 : size(gyro, 2)
    J3 = J3 + J2 * dt + 0.5 * dt^2 * J1 * acc(:,i);
    J2 = J2 + J1 * acc(:,i) * dt;
    J1 = J1 * rodrigues(gyro(:,i) * dt);
    
end


end