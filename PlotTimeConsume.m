function PlotTimeConsume()

a = load('G:\matlab\data\direct\gt\D2_011\4\align1d.txt');

figure(1),subplot(2,2,1);hold on;plot(a(a(:,2) == 0,3));title('align1d');
subplot(2,2,2);hold on;plot(a(a(:,2) == 1,3));title('align frame');
subplot(2,2,3);hold on;plot(a(a(:,2) == 2,3));title('lba');
subplot(2,2,4);hold on;plot(a(a(:,2) == 3,3));title('marg');


b = load('G:\matlab\data\direct\gt\D2_011\4\align_frame.txt');

end