function plotEKF()

close all;

% inputDir = '\\192.168.9.225\original_hand_dataset\zrj\';
inputDir = 'G:\matlab\data\direct\gt\D2_011\4\tbc\ekf\';

aa = load(strcat(inputDir,'delayed_time.txt'));
figure,plot([aa(:,2:3) aa(:,2) + aa(:,3)])
bb = load(strcat(inputDir,'diff_in_dt.txt'));
figure,plot(bb * 1000)

ba_head = load(strcat(inputDir,'ba.txt')); 
bg_head = load(strcat(inputDir,'bg.txt')); 
ba_carrier = load(strcat(inputDir,'ba_carrier.txt')); 
bg_carrier = load(strcat(inputDir,'bg_carrier.txt')); 
figure,subplot(2,2,1);plot(bg_head(:,2:4));title('bg head');subplot(2,2,2);plot(ba_head(:,2:4));title('ba head');
subplot(2,2,3);plot(bg_carrier(:,2:4));title('bg carrier');subplot(2,2,4);plot(ba_carrier(:,2:4));title('ba carrier');


a = load(strcat(inputDir,'hf_predict_output.txt'));
b = load(strcat(inputDir,'hf_filter_output.txt'));
c = load(strcat(inputDir,'hf_ekf_output.txt'));
% figure;hold on;plot(a(:,1),[a(:,2:4)],'-r');legend('predict','hf','ekf');plot(b(:,1),[b(:,2:4)],'-b');plot(c(:,1),[c(:,2:4)],'-g');
% figure;hold on;plot(a(:,1),[a(:,2:4)],'-r');plot(b(:,1),[b(:,2:4)],'-b');plot(c(:,1),[c(:,2:4)],'-g');legend('predict','hf','ekf')
figure;hold on;plot(b(:,1),[b(:,2)],'-b');plot(b(:,1),[b(:,3)],'-b');plot(b(:,1),[b(:,4)],'-b');
               plot(a(:,1),[a(:,2)],'-r');plot(a(:,1),[a(:,3)],'-r');plot(a(:,1),[a(:,4)],'-r');
               plot(c(:,1),[c(:,2)],'-g');plot(c(:,1),[c(:,3)],'-g');plot(c(:,1),[c(:,4)],'-g');
               legend('hf','hf','hf', 'predict','predict','predict', 'ekf','ekf','ekf');



vel_before = load(strcat(inputDir,'cur_v.txt'));
vel_after = load(strcat(inputDir,'updated_v.txt'));
len = min([size(vel_before, 1) size(vel_after, 1)]);
figure,subplot(2,1,1);plot([vel_before(1:len,2:4)],'-r','LineWidth', 2);hold on;plot(vel_after(1:len,2:4),'-b');subplot(2,1,2);plot([vel_before(1:len,2:4) - vel_after(1:len,2:4)]);title('vel');


acc_head_before = load(strcat(inputDir,'cur_acc_head.txt'));
acc_head_after = load(strcat(inputDir,'updated_acc_head.txt'));
ba_head = load(strcat(inputDir,'updated_ba_head.txt'));
len = min([size(acc_head_before, 1) size(acc_head_after, 1) size(ba_head,1)]);
figure,subplot(2,1,1);plot([acc_head_before(1:len,2:4)],'-r','LineWidth', 2);hold on;
plot(acc_head_after(1:len,2:4) + ba_head(1:len,2:4), '-g');
plot(acc_head_after(1:len,2:4),'-b');subplot(2,1,2);plot([acc_head_before(1:len,2:4) - acc_head_after(1:len,2:4)]);title('acc head');

acc_carrier_before = load(strcat(inputDir,'cur_acc_carrier.txt'));
acc_carrier_after = load(strcat(inputDir,'updated_acc_carrier.txt'));
len = min([size(acc_carrier_before, 1) size(acc_carrier_after, 1)]);
figure,subplot(2,1,1);plot([acc_carrier_before(1:len,2:4)],'-r','LineWidth', 2);hold on;plot(acc_carrier_after(1:len,2:4),'-b');subplot(2,1,2);plot([acc_carrier_before(1:len,2:4) - acc_carrier_after(1:len,2:4)]);title('acc carrier');


gyro_head_before = load(strcat(inputDir,'cur_gyro_head.txt'));
gyro_head_after = load(strcat(inputDir,'updated_gyro_head.txt'));
len = min([size(gyro_head_before, 1) size(gyro_head_after, 1)]);
figure,subplot(2,1,1);plot([gyro_head_before(1:len,2:4)],'-r','LineWidth', 2);hold on;plot(gyro_head_after(1:len,2:4),'-b');subplot(2,1,2);plot([gyro_head_before(1:len,2:4) - gyro_head_after(1:len,2:4)]);title('gyro head');


gyro_carrier_before = load(strcat(inputDir,'cur_gyro_carrier.txt'));
gyro_carrier_after = load(strcat(inputDir,'updated_gyro_carrier.txt'));
len = min([size(gyro_carrier_before, 1) size(gyro_carrier_after, 1)]);
figure,subplot(2,1,1);plot([gyro_carrier_before(1:len,2:4)],'-r','LineWidth', 2);hold on;plot(gyro_carrier_after(1:len,2:4),'-b');subplot(2,1,2);plot([gyro_carrier_before(1:len,2:4) - gyro_carrier_after(1:len,2:4)]);title('gyro carrier');


aa_head = load(strcat(inputDir,'updated_aa_head.txt'));
aw_head = load(strcat(inputDir,'updated_aw_head.txt'));
aa_carrier = load(strcat(inputDir,'updated_aa_carrier.txt'));
aw_carrier = load(strcat(inputDir,'updated_aw_carrier.txt'));
figure,subplot(2,2,1);plot([aa_head(:,2:4)]);title('aa head');subplot(2,2,2);plot([aw_head(:,2:4)]);title('aw head');
subplot(2,2,3);plot([aa_carrier(:,2:4)]);title('aa carrier');subplot(2,2,4);plot([aw_carrier(:,2:4)]);title('aw carrier');


err_pvr = load(strcat(inputDir,'err_pvr.txt'));
figure,subplot(3,1,1);plot(err_pvr(:,2:4));title('err p');
subplot(3,1,2);plot(err_pvr(:,5:7));title('err v');
subplot(3,1,3);plot(err_pvr(:,8:10));title('err r');

end