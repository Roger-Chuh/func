function plotEKFJoystick()

close all;

inputDir = '\\192.168.9.225\original_hand_dataset\zrj\';
% inputDir = 'G:\matlab\data\direct\gt\D2_011\4\tbc\ekf\';

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
        b = load(strcat(inputDir,sprintf('hf_filter_output_%d.txt',jid)));
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
    
    
    vel_before = load(strcat(inputDir,sprintf('cur_v_%d.txt',jid)));
    vel_after = load(strcat(inputDir,sprintf('updated_v_%d.txt',jid)));
    len = min([size(vel_before, 1) size(vel_after, 1)]);
    figure,subplot(2,1,1);plot([vel_before(1:len,2:4)],'-r','LineWidth', 2);hold on;plot(vel_after(1:len,2:4),'-b');subplot(2,1,2);plot([vel_before(1:len,2:4) - vel_after(1:len,2:4)]);title('vel');
    
    
    acc_head_before = load(strcat(inputDir,sprintf('cur_acc_head_%d.txt',jid)));
    acc_head_after = load(strcat(inputDir,sprintf('updated_acc_head_%d.txt',jid)));
    ba_head = load(strcat(inputDir,sprintf('updated_ba_head_%d.txt',jid)));
    len = min([size(acc_head_before, 1) size(acc_head_after, 1) size(ba_head,1)]);
    figure,subplot(2,1,1);plot([acc_head_before(1:len,2:4)],'-r','LineWidth', 2);hold on;
    plot(acc_head_after(1:len,2:4) + ba_head(1:len,2:4), '-g');
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
    idx2 = find(a(:,2) == 4);
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
    
    
end

end