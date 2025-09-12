function plotEKFJoystick()

close all;

inputDir = '\\192.168.9.225\original_hand_dataset\zrj\';
% inputDir = 'G:\matlab\data\direct\gt\D2_011\4\tbc\ekf\';
% inputDir = 'G:\matlab\data\direct\gt\D2_011\4\tbc\ekf\Download\';

for jid = 0 : 1
    
    aa = load(strcat(inputDir,sprintf('delayed_time_%d.txt', jid)));
    figure,plot([aa(:,2) [-aa(:,3)] [-aa(:,2) - aa(:,3)] [-aa(:,3) + aa(:,2) + aa(:,3)]]);legend('head - carrier','head - img','carrier - img', 'head - carrier');
    
    %     bb = load(strcat(inputDir,'diff_in_dt.txt'));
    %     figure,plot(bb * 1000)
    
    ba_head = load(strcat(inputDir,sprintf('ba_%d.txt', jid)));
    bg_head = load(strcat(inputDir,sprintf('bg_%d.txt', jid)));
    ba_carrier = load(strcat(inputDir,sprintf('ba_carrier_%d.txt', jid)));
    bg_carrier = load(strcat(inputDir,sprintf('bg_carrier_%d.txt', jid)));
    figure,subplot(2,2,1);plot(bg_head(:,2:4));title('bg head');subplot(2,2,2);plot(ba_head(:,2:4));title('ba head');
    subplot(2,2,3);plot(bg_carrier(:,2:4));title('bg carrier');subplot(2,2,4);plot(ba_carrier(:,2:4));title('ba carrier');
    
    
    try
        a = load(strcat(inputDir,sprintf('hf_predict_output_%d.txt',jid)));
        b = load(strcat(inputDir,sprintf('hf_inte_output_%d.txt',jid)));
        c = load(strcat(inputDir,sprintf('hf_ekf_output_%d.txt',jid)));
        % figure;hold on;plot(a(:,1),[a(:,2:4)],'-r');legend('predict','hf','ekf');plot(b(:,1),[b(:,2:4)],'-b');plot(c(:,1),[c(:,2:4)],'-g');
        % figure;hold on;plot(a(:,1),[a(:,2:4)],'-r');plot(b(:,1),[b(:,2:4)],'-b');plot(c(:,1),[c(:,2:4)],'-g');legend('predict','hf','ekf')
        figure;hold on;plot(b(:,1),[b(:,2)],'-b');plot(b(:,1),[b(:,3)],'-b');plot(b(:,1),[b(:,4)],'-b');
        plot(a(:,1),[a(:,2)],'-r');plot(a(:,1),[a(:,3)],'-r');plot(a(:,1),[a(:,4)],'-r');
        plot(c(:,1),[c(:,2)],'-g');plot(c(:,1),[c(:,3)],'-g');plot(c(:,1),[c(:,4)],'-g');
        legend('hf','hf','hf', 'predict','predict','predict', 'ekf','ekf','ekf');
    catch
        fprintf('sth wrong1\n');
    end
    
    trans_before = load(strcat(inputDir,sprintf('cur_trans_%d.txt',jid)));
    trans_after = load(strcat(inputDir,sprintf('updated_trans_%d.txt',jid)));
    len = min([size(trans_before, 1) size(trans_after, 1)]);
    figure,subplot(2,1,1);plot([trans_before(1:len,2:4)],'-r','LineWidth', 2);hold on;plot(trans_after(1:len,2:4),'-b');subplot(2,1,2);plot([trans_before(1:len,2:4) - trans_after(1:len,2:4)]);title('trans');
    
    
    rot_before = load(strcat(inputDir,sprintf('cur_rot_%d.txt',jid)));
    rot_after = load(strcat(inputDir,sprintf('updated_rot_%d.txt',jid)));
    len = min([size(rot_before, 1) size(rot_after, 1)]);
    figure,subplot(2,1,1);plot([rot_before(1:len,2:4)],'-r','LineWidth', 2);hold on;plot(rot_after(1:len,2:4),'-b');subplot(2,1,2);plot([rot_before(1:len,2:4) - rot_after(1:len,2:4)]);title('rot');
    
    
    vel_before = load(strcat(inputDir,sprintf('cur_v_%d.txt',jid)));
    vel_after = load(strcat(inputDir,sprintf('updated_v_%d.txt',jid)));
    len = min([size(vel_before, 1) size(vel_after, 1)]);
    figure,subplot(2,1,1);plot([vel_before(1:len,2:4)],'-r','LineWidth', 2);hold on;plot(vel_after(1:len,2:4),'-b');subplot(2,1,2);plot([vel_before(1:len,2:4) - vel_after(1:len,2:4)]);title('vel');
    
    
    acc_head_before = load(strcat(inputDir,sprintf('cur_acc_head_%d.txt',jid)));
    acc_head_after = load(strcat(inputDir,sprintf('updated_acc_head_%d.txt',jid)));
    ba_head = load(strcat(inputDir,sprintf('updated_ba_head_%d.txt',jid)));
    len = min([size(acc_head_before, 1) size(acc_head_after, 1) size(ba_head,1)]);
    figure,subplot(2,1,1);plot([acc_head_before(1:len,2:4)],'-r','LineWidth', 2);hold on;
%     plot(acc_head_after(1:len,2:4) + ba_head(1:len,2:4), '-g');
    plot(acc_head_after(1:len,2:4),'-b');subplot(2,1,2);plot([acc_head_before(1:len,2:4) - acc_head_after(1:len,2:4)]);title('acc head');
    
    acc_carrier_before = load(strcat(inputDir,sprintf('cur_acc_carrier_%d.txt',jid)));
    acc_carrier_after = load(strcat(inputDir,sprintf('updated_acc_carrier_%d.txt',jid)));
    len = min([size(acc_carrier_before, 1) size(acc_carrier_after, 1)]);
    figure,subplot(2,1,1);plot([acc_carrier_before(1:len,2:4)],'-r','LineWidth', 2);hold on;plot(acc_carrier_after(1:len,2:4),'-b');subplot(2,1,2);plot([acc_carrier_before(1:len,2:4) - acc_carrier_after(1:len,2:4)]);title('acc carrier');
    
    
    gyro_head_before = load(strcat(inputDir,sprintf('cur_gyro_head_%d.txt',jid)));
    gyro_head_after = load(strcat(inputDir,sprintf('updated_gyro_head_%d.txt',jid)));
    len = min([size(gyro_head_before, 1) size(gyro_head_after, 1)]);
    figure,subplot(2,1,1);plot([gyro_head_before(1:len,2:4)],'-r','LineWidth', 2);hold on;plot(gyro_head_after(1:len,2:4),'-b');subplot(2,1,2);plot([gyro_head_before(1:len,2:4) - gyro_head_after(1:len,2:4)]);title('gyro head');
    
    
    gyro_carrier_before = load(strcat(inputDir,sprintf('cur_gyro_carrier_%d.txt',jid)));
    gyro_carrier_after = load(strcat(inputDir,sprintf('updated_gyro_carrier_%d.txt',jid)));
    len = min([size(gyro_carrier_before, 1) size(gyro_carrier_after, 1)]);
    figure,subplot(2,1,1);plot([gyro_carrier_before(1:len,2:4)],'-r','LineWidth', 2);hold on;plot(gyro_carrier_after(1:len,2:4),'-b');subplot(2,1,2);plot([gyro_carrier_before(1:len,2:4) - gyro_carrier_after(1:len,2:4)]);title('gyro carrier');
    
    
    aa_head = load(strcat(inputDir,sprintf('updated_aa_head_%d.txt',jid)));
    aw_head = load(strcat(inputDir,sprintf('updated_aw_head_%d.txt',jid)));
    aa_carrier = load(strcat(inputDir,sprintf('updated_aa_carrier_%d.txt',jid)));
    aw_carrier = load(strcat(inputDir,sprintf('updated_aw_carrier_%d.txt',jid)));
    figure,subplot(2,2,1);plot([aa_head(:,2:4)]);title('aa head');subplot(2,2,2);plot([aw_head(:,2:4)]);title('aw head');
    subplot(2,2,3);plot([aa_carrier(:,2:4)]);title('aa carrier');subplot(2,2,4);plot([aw_carrier(:,2:4)]);title('aw carrier');
    
    
    
    a = load(fullfile(inputDir, sprintf('lba_info_%d.txt', jid)));
    idx1 = find(a(:,2) == 0);
    idx2 = find(a(:,2) == 1);
    figure,subplot(1,2,1),plot(a(idx1, [4:6]));legend('prior','reproj','trifocal');
    subplot(1,2,2),plot(a(idx2, [3:6]));legend('imu','prior','reproj','trifocal');
    
    
    err_pvr = load(strcat(inputDir,sprintf('err_pvr_%d.txt',jid)));
    figure,subplot(3,1,1);plot(err_pvr(:,2:4));title('err p');
    subplot(3,1,2);plot(err_pvr(:,5:7));title('err v');
    subplot(3,1,3);plot(err_pvr(:,8:10));title('err r');
    
    figure,subplot(2,2,1);plot(err_pvr(:,11:13));title('err gyro carrier');
    subplot(2,2,2);plot(err_pvr(:,14:16));title('err gyro head');
    subplot(2,2,3);plot(err_pvr(:,17:19));title('err acc carrier');
    subplot(2,2,4);plot(err_pvr(:,20:22));title('err acc head');
    
    
    pose_diff = load(strcat(inputDir,sprintf('pose_diff_%d.txt',jid)));
    figure, subplot(4,1,1);plot(pose_diff(:,2:4));title(sprintf('rot diff'));
            subplot(4,1,2);plot(pose_diff(:,5:7));title(sprintf('trans diff'));
            subplot(4,1,3);plot(pose_diff(:,8:9));legend('old reproj','new reproj');
            subplot(4,1,4);plot(pose_diff(:,10));title('vm num');
    
end

a = load(strcat(inputDir,sprintf('processImgOnce.txt')));
b = load(strcat(inputDir,sprintf('processImuOnce_2.txt')));
c = load(strcat(inputDir,sprintf('processImuOnce_0.txt')));
d = load(strcat(inputDir,sprintf('processImuOnce_1.txt')));

figure,subplot(2,4,1);plot(a(:,2));title('img (ms)');
       subplot(2,4,2);plot(b(:,2));title('head imu (ms)');
       subplot(2,4,3);plot(c(:,2));title('controller imu0 (ms)');
       subplot(2,4,4);plot(d(:,2));title('controller imu1 (ms)');
       subplot(2,4,5);hist(a(:,2), 100);title('img (ms)');
       subplot(2,4,6);hist(b(:,2), 100);title('head imu (ms)');
       subplot(2,4,7);hist(c(:,2), 100);title('controller imu0 (ms)');
       subplot(2,4,8);hist(d(:,2), 100);title('controller imu1 (ms)');
       
figure,subplot(4,1,1),plot(a(:,2:14));legend('all','prep','ProvideAllImus', 'GetPrediction2','front end','feed vm to estimator','orca','lba','marg','OptBatch','imu graph','append','trim');
       subplot(4,1,2),hist(a(:,2:14), 20);legend('all','prep','ProvideAllImus', 'GetPrediction2','front end','feed vm to estimator','orca','lba','marg','OptBatch','imu graph','append','trim');
       subplot(4,1,3),hist(a(:,[8 9 12 13 14]), 100);legend('orca(ms)','lba','imu graph','append','trim');
       subplot(4,1,4),plot(a(:,[8 9 12 13 14]));legend('orca(ms)','lba','imu graph','append','trim');
end