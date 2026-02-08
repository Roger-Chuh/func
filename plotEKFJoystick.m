function plotEKFJoystick()

close all;

inputDir = '\\192.168.9.225\original_hand_dataset\zrj\';
% inputDir = 'G:\matlab\data\direct\gt\D2_011\4\tbc\ekf\';
% inputDir = 'G:\matlab\data\direct\gt\D2_011\4\tbc\ekf\Download\';

inputDir = 'G:\matlab\data\direct\gt\D2_011\4\tbc\ekf\Download\shared\';



w_mean = [0.0001 : 0.0001 : 0.8];a = 0.0;k = 0.2;scale = 1-(1-a) * exp(-k .* w_mean);figure,plot(w_mean, scale);

for jid = 0 : 1
    if 0
    a = load(fullfile(inputDir, sprintf('gyro_data_%d.txt', jid)));
    
    imu_window_size = 10;
    gyro_data_filtered = FilterData((a(:,2:4)), imu_window_size);
    gyro_data_variance = FilterDataVariance(a(:,2:4), imu_window_size);
    acc_data_filtered = FilterData((a(:,5:7)), imu_window_size);
    acc_data_variance = FilterDataVariance(a(:,5:7), imu_window_size);
    vel_data_filtered = FilterData((a(:,8:10)), imu_window_size);
    vel_data_variance = FilterDataVariance(a(:,8:10), imu_window_size);
    [~,gyro_raw_norm] = NormalizeVector((a(:,2:4)));
    [~,acc_raw_norm] = NormalizeVector((a(:,5:7)));
    [~,vel_raw_norm] = NormalizeVector((a(:,8:10)));
    [~,gyro_filtered_norm] = NormalizeVector(gyro_data_filtered);
    [~,gyro_variance_norm] = NormalizeVector(gyro_data_variance);
    [~,acc_filtered_norm] = NormalizeVector(acc_data_filtered);
    [~,acc_variance_norm] = NormalizeVector(acc_data_variance);
    [~,vel_filtered_norm] = NormalizeVector(vel_data_filtered);
    [~,vel_variance_norm] = NormalizeVector(vel_data_variance);
    
%     [~,gyro_pc_filtered_norm] = NormalizeVector(a(:,5:7));
    
    figure,subplot(3,1,1),plot(gyro_raw_norm,'-r');hold on;plot(gyro_filtered_norm,'-b');plot(gyro_variance_norm,'-g');legend('raw','filtered','variance');title('gyro data');
           subplot(3,1,2),plot(acc_raw_norm,'-r');hold on;plot(acc_filtered_norm,'-b');plot(acc_variance_norm,'-g');legend('raw','filtered','variance');title('acc data');
           subplot(3,1,3),plot(vel_raw_norm,'-r');hold on;plot(vel_filtered_norm,'-b');plot(vel_variance_norm,'-g');legend('raw','filtered','variance');title('vel data');
    end

    aa = load(strcat(inputDir,sprintf('delayed_time_%d.txt', jid)));
    a = abs([-aa(:,2) - aa(:,3)]);
    index = a < 1000;
    figure,subplot(3,2,[1 2]),plot([aa(index,2) [-aa(index,3)] [-aa(index,2) - aa(index,3)] [-aa(index,3) + aa(index,2) + aa(index,3)] ones(sum(index), 1)]);legend('head - carrier','head - img','carrier - img', 'head - carrier','zero');
    
    %     bb = load(strcat(inputDir,'diff_in_dt.txt'));
    %     figure,plot(bb * 1000)
    
    ba_head = load(strcat(inputDir,sprintf('ba_%d.txt', jid)));
    bg_head = load(strcat(inputDir,sprintf('bg_%d.txt', jid)));
    ba_carrier = load(strcat(inputDir,sprintf('ba_carrier_%d.txt', jid)));
    bg_carrier = load(strcat(inputDir,sprintf('bg_carrier_%d.txt', jid)));
    subplot(3,2,3);plot(bg_head(:,1), bg_head(:,2:4),'x-');title('bg head');subplot(3,2,4);plot(ba_head(:,1),ba_head(:,2:4),'x-');title('ba head');
    subplot(3,2,5);plot(bg_carrier(:,1),bg_carrier(:,2:4),'x-');title('bg carrier');subplot(3,2,6);plot(ba_carrier(:,1),ba_carrier(:,2:4),'x-');title('ba carrier');
    
    zupt = load(strcat(inputDir,sprintf('use_zupt_%d.txt',jid)));
    try
        a = load(strcat(inputDir,sprintf('hf_predict_output_%d.txt',jid)));
        b = load(strcat(inputDir,sprintf('hf_inte_output_%d.txt',jid)));
        c = load(strcat(inputDir,sprintf('hf_ekf_output_%d.txt',jid)));
        d = load(strcat(inputDir,sprintf('output_%d.txt',jid)));
        update_intervals = diff(d(:,1));
        
        aa = GetPoseMat(a);
        bb = GetPoseMat(b);
        cc = GetPoseMat(c);
        
        
        % figure;hold on;plot(a(:,1),[a(:,2:4)],'-r');legend('predict','hf','ekf');plot(b(:,1),[b(:,2:4)],'-b');plot(c(:,1),[c(:,2:4)],'-g');
        % figure;hold on;plot(a(:,1),[a(:,2:4)],'-r');plot(b(:,1),[b(:,2:4)],'-b');plot(c(:,1),[c(:,2:4)],'-g');legend('predict','hf','ekf')
        
        figure;subplot(4,1,1),hold on;plot(b(:,1),[b(:,2)],'-b');plot(b(:,1),[b(:,3)],'-b');plot(b(:,1),[b(:,4)],'-b');
        plot(a(:,1),[a(:,2)],'-r');plot(a(:,1),[a(:,3)],'-r');plot(a(:,1),[a(:,4)],'-r');
        plot(c(:,1),[c(:,2)],'-g');plot(c(:,1),[c(:,3)],'-g');plot(c(:,1),[c(:,4)],'-g');
        legend('hf','hf','hf', 'predict','predict','predict', 'ekf','ekf','ekf');
        subplot(4,1,2);hold on;plot(d(:,1), d(:, 2:4),'-r');plot(b(:,1), b(:, 2:4),'-b');legend('lf','lf','lf','hf inte','hf inte','hf inte');
        subplot(4,1,3);hold on;plot(d(2:end,1), 10 .* update_intervals,'-g');legend('lf update');
        subplot(4,1,4),hold on;plot(bb(:,1),[bb(:,2)],'-b');plot(bb(:,1),[bb(:,3)],'-b');plot(bb(:,1),[bb(:,4)],'-b');
        plot(aa(:,1),[aa(:,2)],'-r');plot(aa(:,1),[aa(:,3)],'-r');plot(aa(:,1),[aa(:,4)],'-r');
        plot(cc(:,1),[cc(:,2)],'-g');plot(cc(:,1),[cc(:,3)],'-g');plot(cc(:,1),[c(:,4)],'-g');
        legend('hf','hf','hf', 'predict','predict','predict', 'ekf','ekf','ekf');
    catch
        fprintf('sth wrong1\n');
    end
    
    trans_before = load(strcat(inputDir,sprintf('cur_trans_%d.txt',jid)));
    trans_after = load(strcat(inputDir,sprintf('updated_trans_%d.txt',jid)));
    len = min([size(trans_before, 1) size(trans_after, 1)]);
    figure,subplot(3,3,1);plot([trans_before(1:len,2:4)],'-r','LineWidth', 2);hold on;plot(trans_after(1:len,2:4),'-b');subplot(3,3,4);plot([trans_before(1:len,2:4) - trans_after(1:len,2:4)]);title('trans');
    
    
    rot_before = load(strcat(inputDir,sprintf('cur_rot_%d.txt',jid)));
    rot_after = load(strcat(inputDir,sprintf('updated_rot_%d.txt',jid)));
    len = min([size(rot_before, 1) size(rot_after, 1)]);
    err_rot = [];
    err_rot_vec_diff = [];
    for id = 1 : len
        err_rot = [err_rot;[(norm(rodrigues(rodrigues(rot_before(id,2:4)) * rodrigues(rot_after(id,2:4))')))]];
        err_rot_vec_diff = [err_rot_vec_diff;[rodrigues(rodrigues(rot_before(id,2:4)) * rodrigues(rot_after(id,2:4))')]'];
    end
    %     figure,subplot(2,1,1);plot([rot_before(1:len,2:4)],'-r','LineWidth', 2);hold on;plot(rot_after(1:len,2:4),'-b');subplot(2,1,2);plot([rot_before(1:len,2:4) - rot_after(1:len,2:4)]);title('rot');
    subplot(3,3,2);plot([rot_before(1:len,2:4)],'-r','LineWidth', 2);hold on;plot(rot_after(1:len,2:4),'-b');subplot(3,3,5);plot(err_rot);title('rot (rad)');subplot(3,3,8);plot(err_rot_vec_diff);title('rot vec diff (rad)');
    
    
    vel_before = load(strcat(inputDir,sprintf('cur_v_%d.txt',jid)));
    vel_after = load(strcat(inputDir,sprintf('updated_v_%d.txt',jid)));
    len = min([size(vel_before, 1) size(vel_after, 1)]);
    subplot(3,3,3);plot([vel_before(1:len,2:4)],'-r','LineWidth', 2);hold on;plot(vel_after(1:len,2:4),'-b');subplot(3,3,6);plot([vel_before(1:len,2:4) - vel_after(1:len,2:4)]);title('vel');
    
    
    acc_head_before = load(strcat(inputDir,sprintf('cur_acc_head_%d.txt',jid)));
    acc_head_after = load(strcat(inputDir,sprintf('updated_acc_head_%d.txt',jid)));
    ba_head = load(strcat(inputDir,sprintf('updated_ba_head_%d.txt',jid)));
    len = min([size(acc_head_before, 1) size(acc_head_after, 1) size(ba_head,1)]);
    figure,subplot(5,4,1);plot([acc_head_before(1:len,2:4)],'-r','LineWidth', 2);hold on;
    %     plot(acc_head_after(1:len,2:4) + ba_head(1:len,2:4), '-g');
    plot(acc_head_after(1:len,2:4),'-b');
    subplot(5,4,5);plot([acc_head_before(1:len,2:4) - acc_head_after(1:len,2:4)]);title('acc head');
    
    acc_carrier_before = load(strcat(inputDir,sprintf('cur_acc_carrier_%d.txt',jid)));
    acc_carrier_after = load(strcat(inputDir,sprintf('updated_acc_carrier_%d.txt',jid)));
    len = min([size(acc_carrier_before, 1) size(acc_carrier_after, 1)]);
    subplot(5,4,2);plot([acc_carrier_before(1:len,2:4)],'-r','LineWidth', 2);hold on;plot(acc_carrier_after(1:len,2:4),'-b');subplot(5,4,6);plot([acc_carrier_before(1:len,2:4) - acc_carrier_after(1:len,2:4)]);title('acc carrier');
    
    
    gyro_head_before = load(strcat(inputDir,sprintf('cur_gyro_head_%d.txt',jid)));
    gyro_head_after = load(strcat(inputDir,sprintf('updated_gyro_head_%d.txt',jid)));
    len = min([size(gyro_head_before, 1) size(gyro_head_after, 1)]);
    subplot(5,4,3);plot([gyro_head_before(1:len,2:4)],'-r','LineWidth', 2);hold on;plot(gyro_head_after(1:len,2:4),'-b');subplot(5,4,7);plot([gyro_head_before(1:len,2:4) - gyro_head_after(1:len,2:4)]);title('gyro head');
    
    
    gyro_carrier_before = load(strcat(inputDir,sprintf('cur_gyro_carrier_%d.txt',jid)));
    gyro_carrier_after = load(strcat(inputDir,sprintf('updated_gyro_carrier_%d.txt',jid)));
    len = min([size(gyro_carrier_before, 1) size(gyro_carrier_after, 1)]);
    subplot(5,4,4);plot([gyro_carrier_before(1:len,2:4)],'-r','LineWidth', 2);hold on;plot(gyro_carrier_after(1:len,2:4),'-b');subplot(5,4,8);plot([gyro_carrier_before(1:len,2:4) - gyro_carrier_after(1:len,2:4)]);title('gyro carrier');
    
    
    aa_head = load(strcat(inputDir,sprintf('updated_aa_head_%d.txt',jid)));
    aw_head = load(strcat(inputDir,sprintf('updated_aw_head_%d.txt',jid)));
    aa_carrier = load(strcat(inputDir,sprintf('updated_aa_carrier_%d.txt',jid)));
    aw_carrier = load(strcat(inputDir,sprintf('updated_aw_carrier_%d.txt',jid)));
    subplot(5,4,9);plot([aa_head(:,2:4)]);title('aa head');subplot(5,4,11);plot([aw_head(:,2:4)]);title('aw head');
    subplot(5,4,10);plot([aa_carrier(:,2:4)]);title('aa carrier');subplot(5,4,12);plot([aw_carrier(:,2:4)]);title('aw carrier');
    
    
    err_pvr = load(strcat(inputDir,sprintf('err_pvr_%d.txt',jid)));
    subplot(5,4,16);plot(err_pvr(:,11:13));title('err gyro carrier');
    subplot(5,4,15);plot(err_pvr(:,14:16));title('err gyro head');
    subplot(5,4,14);plot(err_pvr(:,17:19));title('err acc carrier');
    subplot(5,4,13);plot(err_pvr(:,20:22));title('err acc head');
    
    subplot(5,4,17);plot(300 * diff(acc_head_before(:,2:4)));title('acc head diff');
    subplot(5,4,18);plot(300 * diff(acc_carrier_before(:,2:4)));title('acc carrier diff');
    subplot(5,4,19);plot(300 * diff(gyro_head_before(:,2:4)));title('gyro head diff');
    subplot(5,4,20);plot(300 * diff(gyro_carrier_before(:,2:4)));title('gyro carrier diff');
    
    
    
    figure,subplot(4,1,1);plot(err_pvr(:,2:4));title('err p');
    subplot(4,1,2);plot(err_pvr(:,5:7));title('err v');
    subplot(4,1,3);plot(err_pvr(:,8:10));title('err r (rad)');
    subplot(4,1,4);plot(zupt(:,2));title('use zupt');
    
    
    dx_pvr = load(strcat(inputDir,sprintf('dx_pvr_babg1_babg2_abg12_aba12_%d.txt',jid)));
    figure,subplot(3,4,1);plot(dx_pvr(:,2:4));title('dp');
    subplot(3,4,2);plot(dx_pvr(:,5:7));title('dv');
    subplot(3,4,3);plot(dx_pvr(:,8:10));title('dr (rad)');
    subplot(3,4,5);plot(dx_pvr(:,11:13));title('dba1');
    subplot(3,4,6);plot(dx_pvr(:,14:16));title('dbg1');
    subplot(3,4,7);plot(dx_pvr(:,17:19));title('dba2');
    subplot(3,4,8);plot(dx_pvr(:,20:22));title('dbg2');
    subplot(3,4,9);plot(dx_pvr(:,29:31));title('daba1');
    subplot(3,4,10);plot(dx_pvr(:,23:25));title('dabg1');
    subplot(3,4,11);plot(dx_pvr(:,32:34));title('daba2');
    subplot(3,4,12);plot(dx_pvr(:,26:28));title('dabg2');
    
    
    
    
    
    pose_diff = load(strcat(inputDir,sprintf('pose_diff_%d.txt',jid)));
    figure, subplot(4,1,1);plot(pose_diff(:,1), pose_diff(:,2:4));title(sprintf('rot diff'));
    subplot(4,1,2);plot(pose_diff(:,1), pose_diff(:,5:7));title(sprintf('trans diff'));
    subplot(4,1,3);plot(pose_diff(:,1), pose_diff(:,8:9));legend('old reproj','new reproj');
    subplot(4,1,4);plot(pose_diff(:,1), pose_diff(:,10));title('vm num');
    
    a = load(fullfile(inputDir, sprintf('lba_info_%d.txt', jid)));
    idx1 = find(a(:,2) == 0);
    idx2 = find(a(:,2) == 1);
    [timestamps, x, y] = unique(a(:,1));
    iter = x - 1;
    figure,subplot(5,2,1),plot(a(idx1,1), a(idx1, [4:6]));legend('prior','reproj','trifocal');title('iter 0');
    subplot(5,2,2),plot(a(idx2,1), a(idx2, [3 4 6 7]));legend('imu','prior','reproj','trifocal');title('iter n');
    subplot(5,2,3),plot(a(idx1,1), a(idx1, [5]));title('iter 0 reproj, all mean');
    subplot(5,2,4),hist(a(idx1, [5]), 30);title('iter 0 reproj, all mean');
    subplot(5,2,5),plot(a(idx1,1), a(idx1, [8]));title('iter 0 reproj, cur mean');
    subplot(5,2,6),hist( a(idx1, [8]), 30);title('iter 0 reproj, cur mean');
    subplot(5,2,7),plot(a(idx2,1), a(idx2, [8]));title('iter n reproj, cur mean');
    subplot(5,2,8),hist(a(idx2, [8]), 30);title('iter n reproj, cur mean');
    subplot(5,2,9),plot(a(iter(3:end), 1), a(iter(3:end),2) + 1);title('iter num');
    
    window_size = 50;
    
    w_filtered = FilterData(gyro_carrier_after(:,2:4), window_size);
    vel_filtered = FilterData(vel_after(:,2:4), window_size);
    
    w_filtered_variance = FilterDataVariance(gyro_carrier_after(:,2:4), window_size);
    vel_filtered_variance = FilterDataVariance(vel_after(:,2:4), window_size);
    try
        [~,w_norm] = NormalizeVector(w_filtered);
        [~,vel_norm] = NormalizeVector(vel_filtered);
        [~,w_norm_variance] = NormalizeVector(w_filtered_variance);
        [~,vel_norm_variance] = NormalizeVector(vel_filtered_variance);
        figure,subplot(2,2,1);plot(w_norm);title('w filtered');
        subplot(2,2,2);plot(vel_norm);title('vel filtered');
        subplot(2,2,3);plot(w_norm_variance);title('w filtered variance');
        subplot(2,2,4);plot(vel_norm_variance);title('vel filtered variance');
    catch
        fprintf('no enough data to filter\n');
    end
    
    a = load(fullfile(inputDir, sprintf('Tjxyz_info_%d.txt', jid)));
    figure,subplot(2,1,1);plot(rad2deg(a(:,2:4)),'-b');hold on;plot(rad2deg(a(:,8:10)),'-r');title('rot(deg)');
    subplot(2,1,2);plot(1000 .* a(:,5:7),'-b');hold on;plot(1000 .* a(:,11:13),'-r');title('trans(mm)');
    
    a = load(strcat(inputDir,sprintf('hf_inte_output_%d.txt',jid)));
    b = load(strcat(inputDir,sprintf('output_%d.txt',jid)));
    
    
end
if 1
    
    
    a = load(strcat(inputDir,sprintf('head_imu_minus_image.txt')));
    
    figure,plot(1000 * a(:,2:3));title('head imu minus image (ms)');legend('first', 'second');
    
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
end
function w_filtered = FilterData(gyro_carrier_after, window)

w_filtered = [];
for i = window + 1 : size(gyro_carrier_after,1)
    w_sum = zeros(1,3);
    for j = i - window : i
        w_sum = w_sum + abs(gyro_carrier_after(j, 1:3));
    end
    w_filtered =[w_filtered; w_sum / (window + 1)];
end
end
function w_filtered = FilterDataVariance(gyro_carrier_after, window)

w_filtered = [];
for i = window + 1 : size(gyro_carrier_after,1)
    w_sum = zeros(1,3);
    for j = i - window : i
        w_sum = w_sum + (gyro_carrier_after(j, 1:3));
    end
    w_mean = w_sum ./ (window + 1);
    w_sum_square = zeros(1,3);
    for j = i - window : i
        w_sum_square = w_sum_square + (gyro_carrier_after(j, 1:3) - w_mean).^2;
    end
    w_filtered =[w_filtered; sqrt(w_sum_square / (window + 1))];
end
end
function poseMat = GetPoseMat(data)

poseMat = [];
Twc_stack = {};
for i = 1 : size(data,1)
    data1 = data(i,:);
    xyzw = data1(5:8);
    trans = data1(2:4);
    R = quat2rotm(xyzw([4 1 2 3]));
    rot = rodrigues(R);
    %     poseMat = [poseMat; [data1(1) reshape(R,1,9), trans]];
    poseMat = [poseMat; [data1(1) rot(1) rot(2) rot(3)]];
    Twc_stack{i,1} = [R trans';0 0 0 1];
end
end