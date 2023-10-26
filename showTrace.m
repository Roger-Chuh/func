function showTrace(varargin)



% showTrace('point_trace_statistics.txt')
save_vid = false; true;

close all;

inputDir = 'G:\matlab\data\direct\gt\D2_001';
save_name = 'G:\matlab\data\direct\gt\D2_001\imgs';

inputDir = 'G:\matlab\data\direct\gt\D2_011';
save_name = 'G:\matlab\data\direct\gt\D2_011\imgs';


imgDir = 'G:\matlab\data\direct\gt\D2_011';
inputDir = 'G:\matlab\data\direct\gt\D2_011\4';
save_name = 'G:\matlab\data\direct\gt\D2_011\4\imgs';


% imgDir = 'G:\matlab\data\direct\gt\D2_011';
% inputDir = 'G:\matlab\data\direct\gt\D2_011\3';
% save_name = 'G:\matlab\data\direct\gt\D2_011\3\imgs';


% inputDir = 'G:\matlab\data\direct\gt\D2_004';
% save_name = 'G:\matlab\data\direct\gt\D2_004\imgs';

if nargin < 1
    filename = fullfile(inputDir, 'trace_info.txt');
    mat_name = 'trace_info.mat';
    % trace_info = load(filename);
    % save(mat_name,'trace_info');
elseif nargin == 1
    filename =  fullfile(inputDir,varargin{1});
    mat_name = 'point_trace.mat';
    % point_trace = load(filename);
    % save(mat_name,'point_trace');
else
end





% mat_name = 'trace_info.mat';


load(mat_name);

if(~exist('trace_info','var'))
   trace_info = point_trace; 
end

trace_info0 = trace_info;

trace_info = trace_info(trace_info(:,29) > 0,:);



vkf_id = find(trace_info(:,29) == 2);
vkf_timestamp = unique(trace_info(:,3));


trace_info = trace_info0;


vig{1,1} = double(imread('G:\matlab\data\direct\gt\D2_002\vignette_0.png'))./65535;
vig{2,1} = double(imread('G:\matlab\data\direct\gt\D2_002\vignette_1.png'))./65535;
vig{3,1} = double(imread('G:\matlab\data\direct\gt\D2_002\vignette_2.png'))./65535;
vig{4,1} = double(imread('G:\matlab\data\direct\gt\D2_002\vignette_3.png'))./65535;

for i = 1 : 4
    vig{i,1}(vig{i,1} < 0.15) = nan;
end


try
    gt = load(fullfile(inputDir, 'output.txt'));
    gt_first_time = gt(1,1);
catch
    saglhj = 1;
end

MakeDirIfMissing(save_name);
camInfo = dir(fullfile(imgDir, 'Camera*'));
dirCam0 = dir(fullfile(imgDir, camInfo(1).name,'images','*.bmp'));
dirCam1 = dir(fullfile(imgDir, camInfo(2).name,'images','*.bmp'));
dirCam2 = dir(fullfile(imgDir, camInfo(3).name,'images','*.bmp'));
dirCam3 = dir(fullfile(imgDir, camInfo(4).name,'images','*.bmp'));
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


%% cur_time  |  cur_fid   | host_time   | host_fid   |    Twb    | seed_id   | rho   | rho1d  |  host_cid |  target_cid  | is_outlier |    host_dir   |    host_uv |  target_dir  |   target_uv | first_frame_loss  | cur_frame_loss |  frame_kind  |
%%    1      |     2      |      3      |       4    |   5:10    |   11      |   12  |  13    |    14     |      15      |       16   |    17:19      |   20:21    |    22:24     |     25:26   |        27         |       28       |       29     |


trace_info(:,[14 15 20 21 25 26]) = trace_info(:,[14 15 20 21 25 26]) + 1;
cur_frames = unique(trace_info(:,1));


ArrangeTrace(trace_info);

max_len = 3;
max_time_diff = 10;20;1;100; 0.2; 1; 20; 0.2;


first_time = trace_info(1,3);
last_time = trace_info(end,1);
total_duration = (last_time - first_time)*1e-9;



colorMat = [1 0 0; 0 1 0; 0 0 1; 0 1 1; 1 0 1; 1 1 0;];
colorMat = [colorMat;colorMat;colorMat;colorMat;colorMat;colorMat];



check_time_window = [70 90];
check_time_window = [10 90];
check_time_window = [13 90];
check_time_window = [0 90];
% check_time_window = [10 20];
check_time_window = [-10 11110];
close_window_interval = 20;

figure(98);

has_shown = false;


if save_vid
v = VideoWriter(fullfile(inputDir, 'trace.avi'),'Uncompressed AVI');
v.FrameRate = 30;60; 30; % 10;
open(v);
end


frame_cnt = 1;
for(i = 1 : length(cur_frames))
    
    
    cur_trace = trace_info(find(trace_info(:,1) == cur_frames(i)), :);
    
    
    duration = (cur_trace(1,1) - first_time)*1e-9;
    
    if(duration < check_time_window(1) || duration > check_time_window(2))
        continue;
        
    end
    if( mod(frame_cnt,close_window_interval) == 0)
        figure(98);
        close ;
        figure(98);
        %         set(gcf,'outerposition',get(0,'screensize'));
        has_shown = false;
    end
    
    frame_cnt = frame_cnt + 1;
    
    [camIds_target] = findClosestFrame(cur_frames(i), timestamp1, timestamp2, timestamp3, timestamp4);
    
    host_cids = unique(cur_trace(:,14));
    cur_cids = unique(cur_trace(:,15));
    
    for(k = 1 : 4)
        % 按照相机id去组织trace, 这个相机下共有多少个vm，这里要组织trace矩阵了
        id = find(cur_trace(:, 15) == k);
        trace_temp =  cur_trace(id, :);
        cur_uv{k,1} = {};
        for j = 1 : size(trace_temp,1)
            pid = trace_temp(j, 11);
            cur_time = trace_temp(j, 1);
%             all_vms_by_pid = trace_info(find(trace_info(:,11) == pid & trace_info(:,1) <= cur_time & trace_info(:,15) == k & trace_info(:,14) == k & trace_info(:,16) == 0),:);
            all_vms_by_pid = trace_info(find(trace_info(:,11) == pid & trace_info(:,1) <= cur_time & trace_info(:,15) == k & trace_info(:,16) == 0),:);
            within_range = cur_trace(1,1) - all_vms_by_pid(:,1) < max_time_diff*1e9;
            %             all_vms_by_pid = all_vms_by_pid(max([1 size(all_vms_by_pid,1)-max_len]):end,:);
            all_vms_by_pid = all_vms_by_pid(within_range,:);
            host_time = all_vms_by_pid(1,3);
            
            color_id = find(vkf_timestamp == host_time);
            
            assert(~isempty(color_id))
            
            
            uv = all_vms_by_pid(:,[25 26 11 12 13]);
            cur_uv{k,1}{j,1} = [uv color_id.*ones(size(uv,1),1)];
        end
        
        dirCam_target = dirCams{k};
        target_imgs{k,1} = (imread(fullfile(imgDir, camInfo(k).name, 'images',dirCam_target(camIds_target(k)).name)));
    end
    img_big = uint8(zeros(480*2, 640*2));
    img_big(1:480,1:640) = target_imgs{2};
    img_big(1:480,641:end) = target_imgs{3};
    img_big(481:end,1:640) = target_imgs{1};
    img_big(481:end,641:end) = target_imgs{4};
    imshow(img_big);hold on;
    if ~has_shown
        set(gcf,'outerposition',get(0,'screensize'));
        has_shown = true;
    end
    for k = 1 : 4
        if(k == 1)
            for jj = 1 : length(cur_uv{k,1})
                if(~isempty(cur_uv{k,1}{jj,1}))
                    plot(cur_uv{k,1}{jj,1}(:,1), 480 + cur_uv{k,1}{jj,1}(:,2),'-g');
                    plot(cur_uv{k,1}{jj,1}(end,1), 480 + cur_uv{k,1}{jj,1}(end,2),'*','Color', colorMat(cur_uv{k,1}{jj,1}(end,6),:));
                    text(cur_uv{k,1}{jj,1}(end,1)+2,480 + cur_uv{k,1}{jj,1}(end,2)+2, num2str(cur_uv{k,1}{jj,1}(end,3)),'Color', [1 1 0]);
                end
            end
        elseif k == 2
            for jj = 1 : length(cur_uv{k,1})
                if(~isempty(cur_uv{k,1}{jj,1}))
                    plot(cur_uv{k,1}{jj,1}(:,1), cur_uv{k,1}{jj,1}(:,2),'-g');
                    plot(cur_uv{k,1}{jj,1}(end,1), cur_uv{k,1}{jj,1}(end,2),'*','Color', colorMat(cur_uv{k,1}{jj,1}(end,6),:));
                    text(cur_uv{k,1}{jj,1}(end,1)+2,cur_uv{k,1}{jj,1}(end,2)+2, num2str(cur_uv{k,1}{jj,1}(end,3)),'Color', [1 1 0]);
                end
            end
        elseif k == 3
            for jj = 1 : length(cur_uv{k,1})
                if(~isempty(cur_uv{k,1}{jj,1}))
                    plot(640+cur_uv{k,1}{jj,1}(:,1), cur_uv{k,1}{jj,1}(:,2),'-g');
                    plot(640+cur_uv{k,1}{jj,1}(end,1), cur_uv{k,1}{jj,1}(end,2),'*','Color', colorMat(cur_uv{k,1}{jj,1}(end,6),:));
                    text(640+cur_uv{k,1}{jj,1}(end,1)+2, cur_uv{k,1}{jj,1}(end,2)+2, num2str(cur_uv{k,1}{jj,1}(end,3)),'Color', [1 1 0]);
                end
            end
        else
            for jj = 1 : length(cur_uv{k,1})
                if(~isempty(cur_uv{k,1}{jj,1}))
                    plot(640+cur_uv{k,1}{jj,1}(:,1), 480+cur_uv{k,1}{jj,1}(:,2),'-g');
                    plot(640+cur_uv{k,1}{jj,1}(end,1), 480+cur_uv{k,1}{jj,1}(end,2),'*','Color', colorMat(cur_uv{k,1}{jj,1}(end,6),:));
                    text(640+cur_uv{k,1}{jj,1}(end,1)+2,480 + cur_uv{k,1}{jj,1}(end,2)+2, num2str(cur_uv{k,1}{jj,1}(end,3)),'Color', [1 1 0]);
                end
            end
        end
    end
    title(sprintf('duration: %fs', duration));
%     saveas(gcf, fullfile(save_name, sprintf('%010d.png', round(duration*10000))));
    drawnow;
    
    frame = getframe(gcf);
    image =  imresize(frame.cdata(22:904, 361:1530,:),1);
    imwrite(image,fullfile(save_name,sprintf('point_trace_%010d.png', round(duration*10000))));
    %        image =  imresize(frame.cdata(28:922, 493:1406,:),1);
    if save_vid
        writeVideo(v,image);
    end
    
end

if save_vid
    close(v);
end

end

function ArrangeTrace(trace_info)

% seed_ids = unique(trace_info(:,11));
%
% frame_data = cell(max(seed_ids),1);
% for i = 1 : length(seed_ids)
%     seed_id = seed_ids(i);
%     seed_trace = trace_info(find(trace_info(:,11) == seed_id),:);
%     target_time = unique(seed_trace(:,1));
%     for j = 1 : size(seed_trace,1)
%        target_cid =  seed_trace(j, 15);
%        target_uv = seed_trace(j, 25:26);
%
%     end

% end

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