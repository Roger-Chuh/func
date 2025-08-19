function plotDMVIO
close all;
a = load('G:\matlab\data\direct\gt\D2_011\4\tbc\ekf\output_rmse.txt');figure,subplot(1,2,1),plot(a(:,3));subplot(1,2,2),plot(a(:,2));
a = load('G:\matlab\data\direct\gt\D2_011\4\tbc\ekf\output_frameEnergyTh.txt');figure,plot(a(:,2:3))


end