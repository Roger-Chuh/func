function drawGradDir()
inputDir = 'G:\matlab\data\direct\gt\1';
filename = fullfile(inputDir, 'grad_dir.txt');
idx_sfst = 1;

camInfo = dir(fullfile(inputDir, 'Camera*'));

dirCam0 = dir(fullfile(inputDir, camInfo(1).name,'images','*.bmp'));
dirCam1 = dir(fullfile(inputDir, camInfo(2).name,'images','*.bmp'));
dirCam2 = dir(fullfile(inputDir, camInfo(3).name,'images','*.bmp'));
dirCam3 = dir(fullfile(inputDir, camInfo(4).name,'images','*.bmp'));
dirCams = {dirCam0, dirCam1, dirCam2, dirCam3};

timestamp1 = zeros(length(dirCam0),1);
for i = 1 : length(dirCam0)
    timestamp1(i,1) = str2double(dirCam0(i).name(1:end-4));
end
timestamp2 = zeros(length(dirCam1),1);
for i = 1 : length(dirCam1)
    timestamp2(i,1) = str2double(dirCam1(i).name(1:end-4));
end
timestamp3 = zeros(length(dirCam2),1);
for i = 1 : length(dirCam2)
    timestamp3(i,1) = str2double(dirCam2(i).name(1:end-4));
end
timestamp4 = zeros(length(dirCam3),1);
for i = 1 : length(dirCam3)
    timestamp4(i,1) = str2double(dirCam3(i).name(1:end-4));
end


Data = LoadData(filename);

trace = cell2mat(Data(:,1));

seeds = unique(trace(:,3));


for i = 1 : length(Data)
   if Data{i,2}.target_fid == 2
      cur_time =  Data{i,2}.cur_time;
      [camIds_target] = findClosestFrame(cur_time, timestamp1, timestamp2, timestamp3, timestamp4);
      camIds_first_host = camIds_target-1;
      break;
   end
end
check_seed_id = 32;
ids = find(trace(:,7) == check_seed_id);
for i = 1 : length(seeds)
   seed = seeds(i);
   if seed ~= check_seed_id
       continue;
   end
   idx = find(trace(:,3) == seed);
   trace_list = Data(idx,2);
   host_fid = trace_list{1}.host_fid;
   host_cid = trace_list{1}.host_cid;
   for j = 1 : size(trace_list,1)
       temp = trace_list{j};
       host_fid = trace_list{j}.host_fid;
       target_fid = trace_list{j}.target_fid;
       host_cid = trace_list{j}.host_cid;
       target_cid = trace_list{j}.target_cid;
       dirCam_host = dirCams{host_cid};
       dirCam_target = dirCams{target_cid};
       cur_time = trace_list{j}.cur_time;
       host_fid_delta = host_fid - camIds_first_host(host_cid);
       [camIds_target] = findClosestFrame(cur_time, timestamp1, timestamp2, timestamp3, timestamp4);
       %         [camIds_host] = findClosestFrame(timestamp_host, timestamp1, timestamp2, timestamp3, timestamp4);
       host_image = imread(fullfile(inputDir, camInfo(host_cid).name, 'images',dirCam_host(camIds_first_host(host_cid) + host_fid_delta).name));
       target_image = imread(fullfile(inputDir, camInfo(target_cid).name, 'images',dirCam_target(camIds_target(target_cid)).name));
      if ~(temp.host_pix(1) == 331 && temp.host_pix(2) == 50)
         % continue;
      end
      if abs(temp.angle - temp.angle0) < 20
%           continue;
      end
      if temp.seed_id ~= check_seed_id
          continue;
      end
       
      if target_fid ~= 192 + 1
%          continue; 
      end
      
       figure(100),subplot(1,2,1);imshow(host_image); hold on;plot(temp.host_pix(1), temp.host_pix(2), '.r');title(sprintf('grad angle: [%f %f] deg, res: %d\nhost fid: %d, host cid: %d', temp.angle, temp.angle0, temp.res, host_fid-1, host_cid-1));
       if(temp.res)
           subplot(1,2,2);imshow(target_image); hold on;plot(temp.trace(:,1), temp.trace(:,2), '.r'); hold on;plot(temp.target_pix_best(1), temp.target_pix_best(:,2), '.y');plot(temp.target_pix(1), temp.target_pix(:,2), '.g');title(sprintf('zmssd: %d, seed id: %d, is fully converged: %d\ntarget fid: %d, target cid: %d\ntrace: [%d / %d]',temp.zmssd, temp.seed_id, temp.is_fully_converged, target_fid-1, target_cid-1, j, length(ids)));
       else
           subplot(1,2,2);imshow(target_image); hold on;plot(temp.trace(:,1), temp.trace(:,2), '.r'); hold on;plot(temp.target_pix_best(1), temp.target_pix_best(:,2), '.y');plot(temp.target_pix(1), temp.target_pix(:,2), 'xg');title(sprintf('zmssd: %d, seed id: %d, is fully converged: %d\ntarget fid: %d, target cid: %d\ntrace: [%d / %d]',temp.zmssd, temp.seed_id, temp.is_fully_converged, target_fid-1, target_cid-1, j, length(ids)));
       end
   end
end

end
function Data = LoadData(fileName)
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

offset = 3;
back_offset = 4;

if(~(str2double(a{end-back_offset}) == 0 || str2double(a{end-back_offset}) == 1))
    continue;
end



data.cur_time = str2double(a{1})+1;
data.angle0 = str2double(a{2})+1;
data.target_fid = str2double(a{0+offset})+1;
data.host_fid = str2double(a{1+offset})+1;
data.seed_id = str2double(a{2+offset});
data.host_cid = str2double(a{3+offset})+1;
data.host_pix = [str2double(a{4+offset}) str2double(a{5+offset})]+[1 1];
data.angle = str2double(a{6+offset});
data.target_cid = str2double(a{7+offset})+1;
data.res = str2double(a{end-back_offset});
data.target_pix = [str2double(a{end-2-back_offset}) str2double(a{end-1-back_offset})]+[1 1];
data.trace = str2double(reshape(a(8+offset:end-3-back_offset),2,[])') + 1;


data.is_fully_converged = [str2double(a{end-3})];
data.target_pix_best = [str2double(a{end-2}) str2double(a{end-1})] + [1 1];
data.zmssd = str2double(a{end});

Data{cnt,1} = [data.host_fid data.target_fid data.seed_id data.res data.angle data.angle0 data.seed_id];
Data{cnt,2} = data;
cnt = cnt + 1;
end
fclose(fid);

end
function [camIds] = findClosestFrame(timestamp, timestamp1, timestamp2, timestamp3, timestamp4)

camId1 = findClosest(timestamp, timestamp1);
camId2 = findClosest(timestamp, timestamp2);
camId3 = findClosest(timestamp, timestamp3);
camId4 = findClosest(timestamp, timestamp4);

camIds = [camId1 camId2 camId3 camId4];

end
function ind = findClosest(timestamp, timestamp1)
delta = abs(timestamp - timestamp1);
[~, ind] = min(delta);
end