function plotWorldEKFJoystick()

close all;

inputDir = '\\192.168.9.225\original_hand_dataset\zrj\';
% inputDir = 'G:\matlab\data\direct\gt\D2_011\4\tbc\ekf\';
% inputDir = 'G:\matlab\data\direct\gt\D2_011\4\tbc\ekf\Download\';

inputDir = 'G:\matlab\data\direct\gt\D2_011\4\tbc\ekf\Download\shared\';

for jid = 1 %: 1
    
    cur_a = load(fullfile(inputDir, sprintf('world_cur_a_%d.txt',jid)));
    updated_a = load(fullfile(inputDir, sprintf('world_updated_a_%d.txt',jid)));
    cur_w = load(fullfile(inputDir, sprintf('world_cur_w_%d.txt',jid)));
    updated_w = load(fullfile(inputDir, sprintf('world_updated_w_%d.txt',jid)));
    cur_v = load(fullfile(inputDir, sprintf('world_cur_v_%d.txt',jid)));
    updated_v = load(fullfile(inputDir, sprintf('world_updated_v_%d.txt',jid)));
    
    updated_aw = load(fullfile(inputDir, sprintf('world_updated_aw_%d.txt',jid)));
    
    pose_inte = load(fullfile(inputDir, sprintf('world_hf_inte_results_%d.txt',jid)));
    pose_ekf = load(fullfile(inputDir, sprintf('world_hf_ekf_results_%d.txt',jid)));
    
    len2 = min([size(pose_inte,1) size(pose_ekf,1)]);
    
    
    len = min([size(cur_a,1) size(cur_w,1) size(cur_v,1) size(updated_a,1) size(updated_w,1) size(updated_v,1)]);
    
    figure,
    subplot(3,3,1);plot(cur_a(1:len,2:4),'-r','LineWidth', 2);title('cur a');hold on;plot(updated_a(1:len,2:4),'-b');
    subplot(3,3,4);plot(updated_a(1:len,2:4));title('updated a');
    subplot(3,3,7);plot(updated_a(1:len,2:4) - cur_a(1:len,2:4));title('diff a');
    
    subplot(3,3,2);plot(cur_w(1:len,2:4),'-r','LineWidth', 2);title('cur w');hold on;plot(updated_w(1:len,2:4),'-b');
    subplot(3,3,5);plot(updated_w(1:len,2:4));title('updated w');
    subplot(3,3,8);plot(updated_w(1:len,2:4) - cur_w(1:len,2:4));title('diff w');
    
    subplot(3,3,3);plot(cur_v(1:len,2:4),'-r','LineWidth', 2);title('cur v');hold on;plot(updated_v(1:len,2:4),'-b');
    subplot(3,3,6);plot(updated_v(1:len,2:4));title('updated v');
    subplot(3,3,9);plot(updated_v(1:len,2:4) - cur_v(1:len,2:4));title('diff v');
    
    figure,plot(updated_aw(:,2:4));title('aw');
    
    figure,subplot(3,1,1);plot(pose_inte(1:len2, 2:4),'-r','LineWidth', 2);title('inte trans');hold on;plot(pose_ekf(1:len2, 2:4),'-b');
    subplot(3,1,2);plot(pose_ekf(1:len2, 2:4));title('ekf');
    subplot(3,1,3);plot(pose_ekf(1:len2, 2:4) - pose_inte(1:len2, 2:4));title('diff');
    
    figure,subplot(3,1,1);plot(pose_inte(1:len2, 5:8),'-r','LineWidth', 2);title('inte rot');hold on;plot(pose_ekf(1:len2, 5:8),'-b');
    subplot(3,1,2);plot(pose_ekf(1:len2, 5:8));title('ekf');
    subplot(3,1,3);plot(pose_ekf(1:len2, 5:8) - pose_inte(1:len2, 5:8));title('diff');
    
end


end