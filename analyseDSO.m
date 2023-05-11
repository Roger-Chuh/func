function analyseDSO()

inputDir = 'G:\matlab\data\dso';
data = loadData(fullfile(inputDir, 'idpLog.txt'));



pc_size = 1;


pc_init = cell2mat(data(:,2));
pc_opt = cell2mat(data(:,3));



if 0
    figure, plot3(data{10000,1}{3,1}(:,1), data{10000,1}{3,1}(:,2), data{10000,1}{3,1}(:,3), '.r');axis equal;
    hold on; plot3(data{10000,1}{3,1}(1,1), data{10000,1}{3,1}(1,2), data{10000,1}{3,1}(1,3), 'og','MarkerSize', 5, 'LineWidth',5);
    plot3(data{10000,1}{3,1}(end,1), data{10000,1}{3,1}(end,2), data{10000,1}{3,1}(end,3), 'ob','MarkerSize', 5, 'LineWidth',5);
    legend('trace','init', 'opt');
end


[~,dist1] =NormalizeVector(pc_init);
[~,dist2] =NormalizeVector(pc_opt);
plotQuolity = cell2mat(data(:,5));

trace_num = cell2mat(data(:,6));


figure,hist(trace_num((plotQuolity == 3)),100)

index3 = (dist1 < 5) & (dist2 < 5) & (plotQuolity == 3);

figure,subplot(2,2,1);pcshow(pc_init(index3,:),'MarkerSize', pc_size);
       subplot(2,2,2);pcshow(pc_opt(index3,:),'MarkerSize', pc_size);
subplot(2,2,3);pcshow(pc_init(index3,:), [1 0 0],'MarkerSize', 1);hold on;pcshow(pc_opt(index3,:),[0 1 0],'MarkerSize', 1);legend('init','opt');
idepth_diff = cell2mat(data(:,4));
subplot(2,2,4);hist(idepth_diff(abs(idepth_diff) < 0.3), 200);


index1 = (dist1 < 5) & (dist2 < 5) & (plotQuolity == 1);

figure,subplot(2,2,1);pcshow(pc_init(index1,:),'MarkerSize', pc_size);
       subplot(2,2,2);pcshow(pc_opt(index1,:),'MarkerSize', pc_size);
subplot(2,2,3);pcshow(pc_init(index1,:), [1 0 0],'MarkerSize', 1);hold on;pcshow(pc_opt(index1,:),[0 1 0],'MarkerSize', 1);legend('init','opt');
idepth_diff = cell2mat(data(:,4));
subplot(2,2,4);hist(idepth_diff(abs(idepth_diff) < 0.3), 200);



       

figure,subplot(1,3,1);pcshow(pc_init(index3,:),'MarkerSize', pc_size);subplot(1,3,2);pcshow(pc_opt(index3,:),'MarkerSize', pc_size);
      subplot(1,3,3);pcshow(pc_opt(index1,:),'MarkerSize', pc_size);
      


figure,subplot(3,2,1);pcshow(pc_init(index3,:),'MarkerSize', pc_size);subplot(3,2,2);pcshow(pc_opt(index3,:),'MarkerSize', pc_size);
      subplot(3,2,3);pcshow(pc_opt(index3,:),'MarkerSize', pc_size);subplot(3,2,4);pcshow(pc_init(index1,:),'MarkerSize', pc_size);
      subplot(3,2,5);pcshow(pc_init(index1,:),'MarkerSize', pc_size);subplot(3,2,6);pcshow(pc_opt(index1,:),'MarkerSize', pc_size);

end
function [Data] = loadData(fileName)


fid=fopen(fileName);       %首先打开文本文件coordinate.txt
temp = [];

names = {};
FrameCam = [];
Data = {};
cnt = 1;
while ~feof(fid)    % while循环表示文件指针没到达末尾，则继续
    % 每次读取一行, str是字符串格式
    str = fgetl(fid);     
    a = strsplit(str,' ');
%     idx = find(str == '/');
% %     names = [names; {str(end-21:end)}];
%     names = [names; {str(idx(end-1)+1:end)}];
%    FrameCam = [FrameCam;[str2num(str(1:5)) str2num(str(8))] ];
data = zeros(length(a)-1,1);
for i = 1 : length(a)-1
    data(i,1) = str2double(a{i});
end
tmp = cell(3,1);
tmp{1,1} = data(1);
tmp{2,1} = data(2:3)';
tmp{3,1} = reshape(data(5:end),3,[])';
Data{cnt,1} = tmp;%data;
Data{cnt,2} = tmp{3,1}(1,:);
Data{cnt,3} = tmp{3,1}(end,:);
Data{cnt,5} = data(4);
Data{cnt,6} = size(tmp{3,1},1);
[~,dist] = NormalizeVector(tmp{3,1});
Data{cnt,4} = diff(dist);
cnt = cnt + 1;
end
fclose(fid);





end