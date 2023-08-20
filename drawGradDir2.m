function drawGradDir2(varargin)

global extra_patch_coord patch_offset
close all;

patch_offset = [[0, 0]; [-1, -1]; [1, -1]; [-1, 1]; [-2, 0]; [0, -2]; [2, 0]; [0, 2]];
patch_offset = [[0, 0]; [-1, -1]; [-2, -2]; [1, -1]; [2, -2]; [-1, 1]; [-2, 2]; [2, 2]];

inputDir = 'G:\matlab\data\direct\gt\1';
save_name = 'G:\matlab\data\direct\gt\1\imgs';

inputDir = 'G:\matlab\data\direct\gt\D2_001';
save_name = 'G:\matlab\data\direct\gt\D2_001\imgs';


% inputDir = 'G:\matlab\data\direct\gt\5';
% save_name = 'G:\matlab\data\direct\gt\5\imgs';
%
%
% inputDir = 'G:\matlab\data\direct\gt\D2_011';
% save_name = 'G:\matlab\data\direct\gt\D2_011';
%
inputDir = 'G:\matlab\data\direct\gt\D2_002';
save_name = 'G:\matlab\data\direct\gt\D2_002\imgs';


if nargin < 1
    filename = fullfile(inputDir, 'grad_dir.txt');
    mat_name = 'grad.mat';
    extra_patch_coord = false;
elseif nargin == 1
    filename =  fullfile(inputDir,varargin{1});
    mat_name = 'depth.mat';
    extra_patch_coord = true;
else
end


% inputDir = 'G:\matlab\data\direct\gt\55\5';
% save_name = 'G:\matlab\data\direct\gt\55\5\imgs';


MakeDirIfMissing(save_name);
% filename = fullfile(inputDir, 'grad_dir.txt');
idx_sfst = 1;

draw_pt_only = 0;

trace_single_pid = true; false; true; false; true;
save_trace = false; true;



check_seed_id = 15; 4756;168; 566;485;840;1390;1296; 956; 1110; 236;172; 651; 595;610; 602;16;
check_seed_id = 148; 333;6;4372;333;20;18; 395; 26; 531; 100;1050;

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


% Data = LoadData(filename);
% save(mat_name,'Data');

load(mat_name);

trace = cell2mat(Data(:,1));


cur_fids = unique(trace(:,4));
% cur_fids = cur_fids(end:-1:1);

seeds = unique(trace(:,3));


for i = 1 : size(Data,1)
    if Data{i,2}.target_fid == 2
        cur_time =  Data{i,2}.cur_time;
        [camIds_target] = findClosestFrame(cur_time, timestamp1, timestamp2, timestamp3, timestamp4);
        camIds_first_host = camIds_target-1;
        break;
    end
end


warped_intensities_best_zm = zeros(8,1);
warped_intensities_second_best_zm = zeros(8,1);
if ~trace_single_pid
    plot_vec = [3 1 2 4];
    for i = 1 : length(cur_fids)
        index_in_cur_fid = find(trace(:,4) == cur_fids(i));
        pids_in_cur_fid = trace(index_in_cur_fid, 3);
        [camIds_target] = findClosestFrame(Data{index_in_cur_fid(1),2}.cur_time, timestamp1, timestamp2, timestamp3, timestamp4);
        [camIds_host] = findClosestFrame(Data{index_in_cur_fid(1),2}.host_time, timestamp1, timestamp2, timestamp3, timestamp4);
        target_imgs = cell(4,1);
        figure(1),
        for(k = 1 : 4)
            dirCam_target = dirCams{k};
            target_imgs{k,1} = imread(fullfile(inputDir, camInfo(k).name, 'images',dirCam_target(camIds_target(k)).name));
            subplot(2,2,plot_vec(k));imshow(target_imgs{k,1});hold on;
        end
        
        
        
        for(j = 1 : length(pids_in_cur_fid))
            trace_info = Data{index_in_cur_fid(j),2};
            figure(1),
            dist = zeros(4,1);
            for jj = 1 : 4
                trace_list = squeeze(trace_info.trace(:,jj,:));
                
                if(trace_list(2,1) ~= 0 && trace_list(2,1) ~= 1)
                    dist(jj) = norm(trace_info.A(1:2, jj) - trace_info.B(1:2, jj));
                    
                    if extra_patch_coord
                        patches = trace_info.patch{jj,1};
                        patch_coord_x = squeeze(patches(:,1,:));
                        patch_coord_y = squeeze(patches(:,2,:));
                        %                         subplot(2,2,plot_vec(jj));hold on;plot(patch_coord_x(:), patch_coord_y(:) ,'.w');
                    end
                    
                    if ~draw_pt_only
                        subplot(2,2,plot_vec(jj));hold on;plot(trace_list(2,:),trace_list(3,:),'.r');
                        if(size(trace_list,2)>=3)
                            %                             if(trace_info.is_epi_search_success && size(trace_list,2)>=3)
                            plot(trace_list(2,end-2),trace_list(3,end-2),'.y');
                            if extra_patch_coord
                                if trace_info.is_epi_search_success
                                    offset = 0;
                                else
                                    offset = 1;
                                end
                                subplot(2,2,plot_vec(jj));hold on;plot(patch_coord_x(:, end-2+offset), patch_coord_y(:,end-2+offset) ,'.y');
                                warped_intensities_second_best = interp2(double(target_imgs{jj,1}), patch_coord_x(:, end-2+offset), patch_coord_y(:, end-2+offset), 'linear', 0);
                                warped_intensities_second_best_zm = warped_intensities_second_best - mean(warped_intensities_second_best);
                            end
                        end
                        if extra_patch_coord
                            subplot(2,2,plot_vec(jj));hold on;plot(patch_coord_x(:, end-1+offset), patch_coord_y(:,end-1+offset) ,'.c');
                            %                             subplot(2,2,plot_vec(jj));hold on;plot(patch_coord_x(:, end), patch_coord_y(:,end) ,'.g');
                            warped_intensities_best = interp2(double(target_imgs{jj,1}), patch_coord_x(:, end-1+offset), patch_coord_y(:, end-1+offset), 'linear', 0);
                            warped_intensities_best_zm = warped_intensities_best - mean(warped_intensities_best);
                        end
                        plot(trace_list(2,end-1),trace_list(3,end-1),'.c');
                        plot(trace_list(2,end),trace_list(3,end),'.g');
                        text(trace_list(2,end)+2,trace_list(3,end)-2, num2str(trace_info.seed_id),'Color', [1 1 0]);
                        plot(trace_info.A(1, jj), trace_info.A(2, jj),'om');
                        plot(trace_info.B(1, jj), trace_info.B(2, jj),'ow');
                        if(trace_info.seed_id == check_seed_id)
                            asdfgkj = 1;
                        end
                        
                        
                        
                    else
                        subplot(2,2,plot_vec(jj));hold on;plot(trace_list(2,end),trace_list(3,end),'.g');
                    end
                end
            end
            host_cid = trace_info.host_cid;
            dirCam_host = dirCams{host_cid};
            if ~save_trace
                host_imgs = imread(fullfile(inputDir, camInfo(host_cid).name, 'images',dirCam_host(camIds_host(host_cid)).name));
                host_patch_value = interp2(double(host_imgs), trace_info.host_pix(:,1), trace_info.host_pix(:,2), 'linear', 0);
                host_patch_value_zm = host_patch_value - mean(host_patch_value);
                zmssd_best = sum((host_patch_value_zm - warped_intensities_best_zm).^2);
                zmssd_second_best = sum((host_patch_value_zm - warped_intensities_second_best_zm).^2);
                figure(2),clf;subplot(2,1,1);imshow(host_imgs);hold on;plot(trace_info.host_pix(:,1), trace_info.host_pix(:,2),'.g');title(sprintf('seed id:%d, epi search suc: %d\nfinal suc: %d, [max min] epi len: [%f %f], fully conv: %d\nrho: %f, rho align1d: %f\nquality: %f, step: %f deg\nsigma: %f, b: %f, fail count: %f, [angle-max angle angle0]: [%f %f %f]\n[best_i steps0 steps suc-steps]: [%d %d %d %d]\n[valid / full = ratio] angle: [%f / %f = %f]\nepi len = [%.2f, %.2f, %.2f, %.2f]\n[second best / best] = [%f / %f] = %f',trace_info.seed_id, trace_info.is_epi_search_success, trace_info.is_success, trace_info.max_epi_len, trace_info.min_epi_len, trace_info.is_fully_conv, trace_info.rho, trace_info.rho_align1d, trace_info.quality, trace_info.step, trace_info.sigma, trace_info.b,trace_info.fail_count, trace_info.angle_max, trace_info.angle, trace_info.angle0, trace_info.best_i, trace_info.n_step0, trace_info.n_step, trace_info.success_step, trace_info.epi_angle_valid, trace_info.epi_angle_full, trace_info.epi_angle_valid/trace_info.epi_angle_full,dist(1),dist(2),dist(3),dist(4),zmssd_second_best,zmssd_best, (zmssd_second_best/zmssd_best)));
                subplot(2,1,2);plot(trace_info.rho_list);hold on;plot(length(trace_info.rho_list),trace_info.rho_list(end),'or');
            end
            
            drawnow;
        end
        figure(1),subplot(2,2,1);hold on;title(sprintf('pids; %d',  length(pids_in_cur_fid)));
        saveas(gcf, fullfile(save_name, sprintf('%018d.png', trace_info.cur_time)));
    end
else
    plot_vec = [4 1 2 5];
    ids = find(trace(:,3) == check_seed_id);
    for i = 1 : length(seeds)
        seed = seeds(i);
        if check_seed_id > 0
            if seed ~= check_seed_id
                continue;
            end
        end
        idx = find(trace(:,3) == seed);
        trace_list = Data(idx,2);
        host_time = trace_list{1}.host_time;
        [camIds_host] = findClosestFrame(host_time, timestamp1, timestamp2, timestamp3, timestamp4);
        host_cid = trace_list{1}.host_cid;
        dirCam_host = dirCams{host_cid};
        host_imgs = imread(fullfile(inputDir, camInfo(host_cid).name, 'images',dirCam_host(camIds_host(host_cid)).name));
        figure(1),clf;subplot(2,3,6);imshow(host_imgs);hold on;plot(trace_list{1}.host_pix(:,1), trace_list{1}.host_pix(:,2),'.g');
        host_patch_value = interp2(double(host_imgs), trace_list{1}.host_pix(:,1), trace_list{1}.host_pix(:,2), 'linear', 0);
        host_patch_value_zm = host_patch_value - mean(host_patch_value);
        %             for j = length(trace_list) : -1 : 1
        for j = 1 : length(trace_list)
            
            trace_info = trace_list{j};
            target_time = trace_info.cur_time;
            [camIds_target] = findClosestFrame(target_time, timestamp1, timestamp2, timestamp3, timestamp4);
            target_imgs = cell(4,1);
            figure(1),
            for(k = 1 : 4)
                dirCam_target = dirCams{k};
                target_imgs{k,1} = imread(fullfile(inputDir, camInfo(k).name, 'images',dirCam_target(camIds_target(k)).name));
                subplot(2,3,plot_vec(k));imshow(target_imgs{k,1});hold on;
            end
            dist = zeros(4,1);
            for jj = 1 : 4
                trace_mat = squeeze(trace_info.trace(:,jj,:));
                if(trace_mat(2,1) ~= 0 && trace_mat(2,1) ~= 1)
                    dist(jj) = norm(trace_info.A(1:2, jj) - trace_info.B(1:2, jj));
                    
                    if extra_patch_coord
                        patches = trace_info.patch{jj,1};
                        patch_coord_x = squeeze(patches(:,1,:));
                        patch_coord_y = squeeze(patches(:,2,:));
                        %                     subplot(2,2,plot_vec(jj));hold on;plot(patch_coord_x(:), patch_coord_y(:) ,'.w');
                    end
                    
                    
                    
                    subplot(2,3,plot_vec(jj));hold on;plot(trace_mat(2,:),trace_mat(3,:),'.r');
                    if(size(trace_mat,2)>=3)
                        %                         if(trace_info.is_epi_search_success && size(trace_mat,2)>=3)
                        plot(trace_mat(2,end-2),trace_mat(3,end-2),'.y');
                        if extra_patch_coord
                            if trace_info.is_epi_search_success
                                offset = 0;
                            else
                                offset = 1;
                            end
                            subplot(2,3,plot_vec(jj));hold on;plot(patch_coord_x(:, end-2+offset), patch_coord_y(:,end-2+offset) ,'.y');
                            warped_intensities_second_best = interp2(double(target_imgs{jj,1}), patch_coord_x(:, end-2+offset), patch_coord_y(:, end-2+offset), 'linear', 0);
                            warped_intensities_second_best_zm = warped_intensities_second_best - mean(warped_intensities_second_best);
                        end
                    end
                    if extra_patch_coord
                        subplot(2,3,plot_vec(jj));hold on;plot(patch_coord_x(:, end-1+offset), patch_coord_y(:,end-1+offset) ,'.c');
                        %                         subplot(2,2,plot_vec(jj));hold on;plot(patch_coord_x(:, end), patch_coord_y(:,end) ,'.g');
                        warped_intensities_best = interp2(double(target_imgs{jj,1}), patch_coord_x(:, end-1+offset), patch_coord_y(:, end-1+offset), 'linear', 0);
                        warped_intensities_best_zm = warped_intensities_best - mean(warped_intensities_best);
                    end
                    plot(trace_mat(2,end-1),trace_mat(3,end-1),'.c');
                    plot(trace_mat(2,end),trace_mat(3,end),'.g');
                    plot(trace_info.A(1, jj), trace_info.A(2, jj),'om');
                    plot(trace_info.B(1, jj), trace_info.B(2, jj),'ow');
                end
            end
            zmssd_best = sum((host_patch_value_zm - warped_intensities_best_zm).^2);
            zmssd_second_best = sum((host_patch_value_zm - warped_intensities_second_best_zm).^2);
            
            subplot(2,3,6);title(sprintf('seed id： %d, trace len: [%d / %d], epi search suc: %d, final suc: %d\n[max min] epi len: [%f %f], fully conv: %d\nrho: %f, rho align1d: %f\nquality: %f, step: %f deg\nsigma: %f, b: %f\nfail count: %f, [angle-max angle angle0]: [%f %f %f]\n[best_i steps0 steps suc-steps]: [%d %d %d %d]\n[valid / full = ratio] angle: [%f / %f = %f]\nepi len = [%.2f, %.2f, %.2f, %.2f]\n[second best / best] = [%f / %f] = %f',seed, j, length(trace_list), trace_info.is_epi_search_success, trace_info.is_success, trace_info.max_epi_len, trace_info.min_epi_len, trace_info.is_fully_conv, trace_info.rho, trace_info.rho_align1d, trace_info.quality, trace_info.step, trace_info.sigma,trace_info.b,trace_info.fail_count, trace_info.angle_max, trace_info.angle, trace_info.angle0, trace_info.best_i, trace_info.n_step0, trace_info.n_step, trace_info.success_step, trace_info.epi_angle_valid, trace_info.epi_angle_full, trace_info.epi_angle_valid/trace_info.epi_angle_full,dist(1),dist(2),dist(3),dist(4),zmssd_second_best,zmssd_best, zmssd_second_best/zmssd_best));
            subplot(2,3,3);plot(trace_info.rho_list);hold on;plot(length(trace_info.rho_list),trace_info.rho_list(end),'or');
            saveas(gcf, fullfile(save_name, sprintf('trace_%018d.png', trace_info.cur_time)));
        end
    end
    
    
end




% check_seed_id = 32;
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
        
        figure(100),subplot(1,2,1);imshow(host_image); hold on;plot(temp.host_pix(:,1), temp.host_pix(:,2), '.r');title(sprintf('grad angle: [%f %f] deg, res: %d\nhost fid: %d, host cid: %d', temp.angle, temp.angle0, temp.res, host_fid-1, host_cid-1));
        if(temp.res)
            subplot(1,2,2);imshow(target_image); hold on;plot(temp.trace(:,1), temp.trace(:,2), '.r'); hold on;plot(temp.target_pix_best(1), temp.target_pix_best(:,2), '.y');plot(temp.target_pix(1), temp.target_pix(:,2), '.g');title(sprintf('zmssd: %d, seed id: %d, is fully converged: %d\ntarget fid: %d, target cid: %d\ntrace: [%d / %d]',temp.zmssd, temp.seed_id, temp.is_fully_converged, target_fid-1, target_cid-1, j, length(ids)));
        else
            subplot(1,2,2);imshow(target_image); hold on;plot(temp.trace(:,1), temp.trace(:,2), '.r'); hold on;plot(temp.target_pix_best(1), temp.target_pix_best(:,2), '.y');plot(temp.target_pix(1), temp.target_pix(:,2), 'xg');title(sprintf('zmssd: %d, seed id: %d, is fully converged: %d\ntarget fid: %d, target cid: %d\ntrace: [%d / %d]',temp.zmssd, temp.seed_id, temp.is_fully_converged, target_fid-1, target_cid-1, j, length(ids)));
        end
    end
end

end
function Data = LoadData(fileName)
global extra_patch_coord patch_offset

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
    
    offset = 19;
    back_offset = 4;
    
    % if(~(str2double(a{end-back_offset}) == 0 || str2double(a{end-back_offset}) == 1))
    %     continue;
    % end
    
    
    
    data.cur_time = str2double(a{1});
    data.target_fid = str2double(a{2})+1;
    data.host_fid = str2double(a{3})+1;
    data.seed_id = str2double(a{4});
    data.host_cid = str2double(a{5})+1;
    data.host_pix = [str2double(a{6}) str2double(a{7})]+[1 1] + patch_offset;
    data.rho = 1/str2double(a{8});
    data.rho_align1d = 1/str2double(a{9});
    data.do_search = str2double(a{10});
    
    data.is_success = str2double(a{end});
    data.max_epi_len = str2double(a{end-1});
    data.is_fully_conv = str2double(a{end-2});
    data.is_epi_search_success = str2double(a{end-3});
    data.host_time = str2double(a{end-4});
    data.quality = str2double(a{end-5});
    data.step = rad2deg(str2double(a{end-6}));
    data.min_epi_len = str2double(a{end-7});
    data.b = str2double(a{end-8});
    data.fail_count = str2double(a{end-9});
    data.angle0 = str2double(a{end-10});
    data.angle = str2double(a{end-11});
    data.angle_max = str2double(a{end-12});
    data.success_step = str2double(a{end-13});
    data.n_step = str2double(a{end-14});
    data.n_step0 = str2double(a{end-15});
    data.best_i = str2double(a{end-16});
    data.epi_angle_full = rad2deg(str2double(a{end-17}));
    data.epi_angle_valid = rad2deg(str2double(a{end-18}));
    data.sigma = sqrt(str2double(a{end-19}));
    
    epi_search_data = str2double(a(11:end-(offset+1))) + 1;
    if ~extra_patch_coord
        epi_search_data_mat = reshape(epi_search_data, 5+4, 4, []);
        epi_search_data_mat(4:5,:,:) = epi_search_data_mat(4:5,:,:)-1;
    else
        epi_search_data_mat = reshape(epi_search_data, 5+4 + 16, 4, []);
        epi_search_data_mat(4:5,:,:) = epi_search_data_mat(4:5,:,:)-1;
        patchs = epi_search_data_mat(10:25,:,:);
        for k = 1 : 4
            data.patch{k,1} = reshape(squeeze(patchs(:,k,:)),8,2,[]);
        end
    end
    data.trace = epi_search_data_mat;
    
    data.rho_list = squeeze(epi_search_data_mat(5,1,:));
    
    data.A = epi_search_data_mat(6:7,:,1);
    data.B = epi_search_data_mat(8:9,:,1);
    
    
    % if(data.do_search)
    %     data.rho_list = 1./squeeze(epi_search_data_mat(5,1,:));
    % else
    %     afs = 1;
    %
    % end
    
    
    Data{cnt,1} = [data.host_fid data.target_fid data.seed_id data.cur_time];
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