function testDirectH22()


close all;

a = load('G:\matlab\data\direct\gt\D2_011\4\align1d_H22_statistics.txt');

id = find(a(:,1) == 1);
aa  = a(id,:);
length(unique(aa(:,2)))


pids = unique(a(:,2));

figure, subplot(1,2,1),hist(a(:,3), 100);subplot(1,2,2),plot(sort(a(:,3)))
figure,plot(unique(a(:,3)))
H22 = [];
for i = 1 : length(pids)
   pid = pids(i); 
   index = find(a(:,2) == pid); 
   temp = a(index,:);
   fid_num = length(unique(temp(:,1)));
%    assert(fid_num == size(temp,1));
   H22 = [H22; unique(temp(:,3))];
end

figure,subplot(1,2,1),plot(H22);subplot(1,2,2);hist(H22, 50);

end