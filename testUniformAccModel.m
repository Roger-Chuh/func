function testUniformAccModel()
vals = [1 4];
max_x = 0.8 * 0.5 * 5;
max_y = 0.8 * 0.5 * 5;
max_z = 0.8 * 0.5 * 5;
k = 2 * pi / 20;
time_list = 0.033:0.033:20;
cnt = 1;
for t = time_list
    p(cnt,:) = [max_x * cos(k * t), max_y * sin(k * t), max_z * sin(10 * k * t)];
    cnt = cnt + 1;
end

pos = sin(time_list);
velo = cos(time_list);
acc = -sin(time_list);

p0 = pos(1);
velo0 = velo(1);
acc0 = acc(1);

for i = 1 : length(acc)
    
end


dx1 = pos(10) - pos(2);
dx2 = pos(20) - pos(10);

dt1 = time_list(10) - time_list(2);
dt2 = time_list(20) - time_list(10);

[acc, v0, v1, v2] = CalcAccVelo(dx1, dx2, dt1, dt2);


dx1_check = v0 * dt1 + 0.5 * acc * dt1^2;
dx1_check - dx1

dx2_check = (v0 + acc * dt1)* dt2 + 0.5 * acc * dt2^2;
dx2_check - dx2

dx_check = v0 * (dt1 + dt2) + 0.5 * acc * (dt1 + dt2)^2;
dx_check - (dx1 + dx2)
end
function [acc, v0, v1, v2] = CalcAccVelo(dx1, dx2, dt1, dt2)
  acc = 2 * (dx2 - dt2 / dt1 * dx1) / (dt2 * (dt1 + dt2));
  v0 = (dx1 - 0.5 * acc * dt1 * dt1) / dt1;
  v1 = v0 + acc * (dt1);
  v2 = v0 + acc * (dt1 + dt2);
end