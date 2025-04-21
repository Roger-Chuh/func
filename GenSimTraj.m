function GenSimTraj()

max_x = 0.8 * 0.5 * 5;
max_y = 0.8 * 0.5 * 5;
max_z = 0.8 * 0.5 * 3;
k = 2 * pi / 20;
p = [max_x * cos(k * t), max_y * sin(k * t), max_z * sin(10 * k * t)];
dp = [-max_x * sin(k * t) * k, max_y * cos(kt) * k, max_z * cos(10 * k * t) * 10 * k];
ddp = [-max_x * cos(k * t) * k * k, -max_y * sin(k * t) * k * k, -max_z * sin(10 * k * t) * 10 * 10 * k * k];

end