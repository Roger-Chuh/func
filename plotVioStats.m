function a = plotVioStats()

close all;
a = load('G:\matlab\data\direct\gt\D2_011\4\tbc\align_frame_info.txt');figure, plot(a(a(:,2) == 0 & a(:,3) == 3,4));b = load('G:\matlab\data\direct\gt\D2_011\4\tbc\lba_info.txt');figure,subplot(2,1,1); plot(b(b(:,1) == 1 & b(:,2) == 2,3));subplot(2,1,2); plot(b(b(:,1) == 1 & b(:,2) == 2,4:8));legend('imu','prior','blob','trifocal blob','reproj')
try
    ba = load('G:\matlab\data\direct\gt\D2_011\4\tbc\ba.txt');
    if (size(ba,2) <= 8)
        figure,subplot(2,2,1);plot(ba(:,2:4));subplot(2,2,2);plot(ba(:,5:7));subplot(2, 2, 3);plot(ba(:,8),'-m');subplot(2, 2, 4);hist(ba(:,8),50);
    else
        figure,subplot(3,2,1);plot(ba(:,2:4));subplot(3,2,3);plot(ba(:,5:7));subplot(3, 2, 2);plot(ba(:,11),'-m');subplot(3, 2, 4);hist(ba(:,11),50);subplot(3, 2, 5);plot(ba(:,8:10));
    end
    figure,plot(diff(ba(:,1)))
    bg = load('G:\matlab\data\direct\gt\D2_011\4\tbc\bg.txt');
    if (size(bg,2) <= 7)
        figure,subplot(2,1,1);plot(bg(:,2:4));subplot(2,1,2);plot(bg(:,5:7))
    else
        figure,subplot(3,1,1);plot(bg(:,2:4));subplot(3,1,2);plot(bg(:,5:7));subplot(3,1,3);plot(bg(:,8:10))
    end
catch
    fprintf('sth wrong\n');
    ba = load('G:\matlab\data\direct\gt\D2_011\4\tbc\ba.txt');
    bg = load('G:\matlab\data\direct\gt\D2_011\4\tbc\bg.txt');
    figure,plot(ba(:,1), ba(:,2:4));
    figure,plot(bg(:,1), bg(:,2:4));
end
a = load('G:\matlab\data\direct\gt\D2_011\4\tbc\scale.txt');figure,plot([a(:,2:3)]);figure,subplot(2,2,1);hold on;plot(a(:,[7 9]));plot(a(:,[8 10]));title('cam0');subplot(2,2,2);hold on;plot(a(:,[12 14]));plot(a(:,[13 15]));title('cam1');subplot(2,2,3);hold on;plot(a(:,[17 19]));plot(a(:,[18 20]));title('cam2');subplot(2,2,4);hold on;plot(a(:,[22 24]));plot(a(:,[23 25]));title('cam3');figure,plot(a(:,[5 26]))
[poseMat0, time0] = plotTraj('G:\matlab\data\direct\gt\D2_011\4\tbc\raw_output.txt');close;
[poseMat, time] = plotTraj('G:\matlab\data\direct\gt\D2_011\4\tbc\all_output.txt');close;
[a1, a2, a3] = intersect(time0, time, 'rows');
figure,subplot(2,2,1),plot(time0, poseMat0(:,10:12),'-r');hold on;plot(time, poseMat(:,10:12),'-b');title('trans');subplot(2,2,2);plot(time0, poseMat0(:,13:15),'-r');hold on;plot(time, poseMat(:,13:15),'-b');title('vel');
subplot(2,2,3),plot(a1, poseMat0(a2,10:12) - poseMat(a3,10:12));title('trans');subplot(2,2,4);plot(a1, poseMat0(a2,13:15) - poseMat(a3,13:15));title('vel');
figure(99),subplot(4,1,1),plot(time, poseMat(:,10:12));title('trans');subplot(4,1,3);plot(time, poseMat(:,13:15));title('vel');
[poseMat1, time1] = plotTraj('G:\matlab\data\direct\gt\D2_011\4\tbc\full_output.txt');close;
try
[poseMat2, time2] = plotTraj('G:\matlab\data\direct\gt\D2_011\4\tbc\short_output.txt');close;
figure(99),subplot(4,1,2),hold on;plot(time1,poseMat1(:,10:12),'-r');plot(time2,poseMat2(:,10:12), '-b');title('trans');
subplot(4,1,4);hold on;plot(time1, poseMat1(:,13:15), '-r');plot(time2, poseMat2(:,13:15), '-b');title('vel');
catch
   sprintf('no short frame pose\n'); 
end
% poseMat = plotTraj('G:\matlab\data\direct\gt\D2_011\4\tbc\global_output.txt');
% figure,plot(poseMat(:,10:12));
a = load('G:\matlab\data\direct\gt\D2_011\4\tbc\elevator_acc.txt');
figure,subplot(2,1,1),plot(a(:,2));subplot(2,1,2),hist(a(:,2), 40);

% poseMat = plotTraj('G:\matlab\data\direct\gt\D2_011\4\tbc\hf_filter_output.txt');
% figure,plot(poseMat(:,10:12));

a = load('G:\matlab\data\direct\gt\D2_011\4\tbc\extr.txt');
cols = (size(a, 2)-1)/2;
est = a(:,2:1+cols);
prior = a(:,2+cols:end);
if 0
    figure;subplot(2,1,1);plot(a(:,2:end));subplot(2,1,2);plot(est - prior);
else
    figure,
    for i = 1 : 6
       subplot(2, 3, i);plot([est(:,i) prior(:,i)]);
    end
end

a = load('G:\matlab\data\direct\gt\D2_011\4\tbc\vehicle_imu.txt');
if size(a,1) > 10
figure,subplot(2,2,1);plot(a(:,1), a(:,2:4)); title('gyro');
       subplot(2,2,2);plot(a(:,1), a(:,5:7)); title('acc');
       subplot(2,2,3);plot(a(:,1), a(:,8:10)); title('d_gyro');
       subplot(2,2,4);plot(a(:,1), a(:,11:13)); title('d_acc');
end
% clear;

end