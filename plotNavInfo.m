function plotNavInfo()

close all;

inputDir = '\\192.168.9.225\original_hand_dataset\zrj\';
% inputDir = 'G:\matlab\data\direct\gt\D2_011\4\tbc\ekf\';
% inputDir = 'G:\matlab\data\direct\gt\D2_011\4\tbc\ekf\Download\';

for jid = 0 : 1
    
    
    ba_head = load(strcat(inputDir,sprintf('ba_%d.txt', jid)));
    bg_head = load(strcat(inputDir,sprintf('bg_%d.txt', jid)));
    ba_carrier = load(strcat(inputDir,sprintf('ba_carrier_%d.txt', jid)));
    bg_carrier = load(strcat(inputDir,sprintf('bg_carrier_%d.txt', jid)));
    figure(100),subplot(2,2,1);plot(bg_head(:,2:4));title('bg head');subplot(2,2,2);plot(ba_head(:,2:4));title('ba head');
    subplot(2,2,3);plot(bg_carrier(:,2:4));title('bg carrier');subplot(2,2,4);plot(ba_carrier(:,2:4));title('ba carrier');
    
    
    data = load(strcat(inputDir,sprintf('nav_info_%d.txt', jid)));
    [a, b, c] = unique(data(:,3));
    
    figure(1);subplot(2,2,1);hold on;title('rot');
    subplot(2,2,2);hold on;title('trans');
    subplot(2,2,3);hold on;title('vel');
    for fid = 1 : length(a)
        idx = find(abs(data(:, 3) - a(fid)) < 0.00001);
        temp = data(idx,:);
        [aa, bb, cc] = unique(temp(:,1));
        timestamps = temp(bb, 1);
        figure(1), subplot(2,2,1);hold on;plot(timestamps, temp(bb,4:6),'-r');
        subplot(2,2,2);hold on;plot(timestamps, temp(bb,7:9),'-r');
        subplot(2,2,3);hold on;plot(timestamps, temp(bb,10:12),'-r');
        subplot(2,2,4);%plot(timestamps, temp(:,10:12));title(sprintf('timestamp: %0.2f, duration: %0.2f(s)',a(fid), a(fid) - data(1,1)));
    end
    
end



end