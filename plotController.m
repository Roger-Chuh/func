function plotController()

close all;

% a = load('G:\matlab\data\direct\gt\D2_011\4\tbc\ekf\lba_info.txt');
prefix = 'G:\matlab\data\direct\gt\D2_011\4\tbc\ekf';
prefix = '//192.168.9.225/original_hand_dataset/zrj';

for jid = 0 : 1
a = load(fullfile(prefix, sprintf('lba_info_%d.txt', jid)));

idx1 = find(a(:,2) == 0);
idx2 = find(a(:,2) == 4);

figure,subplot(1,2,1),plot(a(idx1, [4:6]));legend('prior','reproj','trifocal');
       subplot(1,2,2),plot(a(idx2, [3:6]));legend('imu','prior','reproj','trifocal');
end
end