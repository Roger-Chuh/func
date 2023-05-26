function AnalyseDepthFilter2()

% close all;

inputDir = 'C:\Users\Roger\Desktop\vins-fusion-fisheye\yvr\april\controller_2022-04-12-05-46-16';
inputDir = 'G:\matlab\data\tag\controller_2022-05-31-06-33-32';
inputDir = 'G:\matlab\data\direct\gt\1';
% inputDir = 'G:\matlab\data\direct\gt\ke\2';

% inputDir = 'G:\matlab\data\direct\gt\D2_001';
% inputDir = 'G:\matlab\data\direct\gt\self';
% inputDir = 'G:\matlab\data\direct\gt\ke\3';

scale = 1;
camInfo = dir(fullfile(inputDir, 'Camera*'));

% depth_filter_data = load(fullfile(inputDir, 'depth_filter_statistics.txt'));

%  [Data] = readNames(fullfile(inputDir, 'depth_filter_statistics.txt'));
 
% save('Data.mat','Data');

check_pid = 3077;  2454;106;454;
check_fid = 14 - 2;
check_fid = max([0 check_fid]);


draw = 0;1; 0;
draw2 = 0; 1; 0; 1; 0; 1; 0; 1;
draw3 = 0; 1; 0;
big_fig = 0;1;0;

 load('Data.mat');
 
 pc = load(fullfile(inputDir, 'landmark_statistics.txt'));
 idepth = pc(:,2);
 landmark = pc(:,3:5);
 sigma = pc(:,6);
 goodness = pc(:,7);
 [~,depths_] = NormalizeVector(landmark);
%  figure(4000),pcshow(landmark(goodness == 1 & depths_<2 & sigma < 1,:),'MarkerSize', 100);
 figure(4000),pcshow(landmark(goodness == 1 & depths_<5.5 & sigma < 0.1 & idepth > -0.8,:),'MarkerSize', 100);
 
 depth_filter_data = extractData(Data);

% depth_filter_data(:,4:12) = depth_filter_data(:,4:12) + 1;


check_cam_id = 1000; -1; 1; -1; 

dirCam0 = dir(fullfile(inputDir, camInfo(1).name,'images','*.bmp'));
dirCam1 = dir(fullfile(inputDir, camInfo(2).name,'images','*.bmp'));
dirCam2 = dir(fullfile(inputDir, camInfo(3).name,'images','*.bmp'));
dirCam3 = dir(fullfile(inputDir, camInfo(4).name,'images','*.bmp'));
dirCams = {dirCam0, dirCam1, dirCam2, dirCam3};
% assert(length(dirCam0) == length(dirCam3));
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

% % frameNum = unique(depth_filter_data(:,1));
% 
% % featIndAll = find(abs(depth_filter_data(:,1) - timestamp(1)) < 1000);
% % featIdAll = depth_filter_data(featIndAll, 4);
% % % hostCoord = depth_filter_data(featIndAll, 5:6);
% % % hostCoord = zeros
% % featMatX = zeros(max(depth_filter_data(featIndAll,4)), length(frameNum));
% % featMatY = zeros(max(depth_filter_data(featIndAll,4)), length(frameNum));
% % featMatAX = zeros(max(depth_filter_data(featIndAll,4)), length(frameNum));
% % featMatBX = zeros(max(depth_filter_data(featIndAll,4)), length(frameNum));
% % featMatX(featIdAll,1) = hostCoord(:,1);
% % featMatY(featIdAll,1) = hostCoord(:,2);
% % featMatAX(featIdAll,1) = hostCoord(:,1);
% % featMatAY(featIdAll,1) = hostCoord(:,2);
% % featMatBX(featIdAll,1) = hostCoord(:,1);
% % featMatBY(featIdAll,1) = hostCoord(:,2);
% 
% for imgId = 1 : length(frameNum) %length(dirCam0)
% %     img0 = imread(fullfile(inputDir, camInfo(1).name, 'images',dirCam0(imgId).name));
% %     img3 = imread(fullfile(inputDir, camInfo(4).name, 'images',dirCam3(imgId).name));
%     feats = find(abs(depth_filter_data(:,1) - timestamp(imgId)) < 1000);
%     featid = depth_filter_data(feats, 4);
%     hostCoord(featid,:) = depth_filter_data(feats, 5:6);
%     px_cur = depth_filter_data(feats, 11:12);
%     px_A = depth_filter_data(feats, 7:8);
%     px_B = depth_filter_data(feats, 9:10);
%     zList = depth_filter_data(feats, 13);
%     featMatX(featid,imgId) = px_cur(:,1);
%     featMatY(featid,imgId) = px_cur(:,2);
%     featMatAX(featid,imgId) = px_A(:,1);
%     featMatAY(featid,imgId) = px_A(:,2);
%     featMatBX(featid,imgId) = px_B(:,1);
%     featMatBY(featid,imgId) = px_B(:,2);
%     zMat(featid,imgId) = zList;
%     %     for feat_id = 1 : feats
%     %     end
%     
%     
%     
% end
% 
% 
% validId = find(sum(featMatX > 0,2) > 0);



index = [];
for k = 1 : length(depth_filter_data)
    
    index = [index; [k, depth_filter_data{k}.frame_id;]];
    
    
end






frameNum = unique(index(:,2))';
depth_filter_data_all = depth_filter_data;




pids = [];
for kkk = 1 : size(depth_filter_data_all,1)
pids = [pids; depth_filter_data_all{kkk}.pid];
end


check_pid_count = 0;
inds_ = find(pids == check_pid)

zMat = [];
xyzCell =cell(9999999,1);
pair_num = 0;
base_fig = 200;
depth_bounds = cell(1000000,1);
for frame_id = frameNum % 1 : length(depth_filter_data)
    if(frame_id < check_fid)
        continue;
    end
    ind = find(index(:,2) == frame_id);
    depth_filter_data = depth_filter_data_all(ind);
    host = [];
    target = [];
    imgPair = {};
%     host_points = cell(4, 1000);target_points = cell(4,1000);
    track_stack = {};
    for feat_id =  1 : length(depth_filter_data)
        
        temp = depth_filter_data{feat_id};
        timestamp = temp.ts;
        timestamp_host = temp.host_ts;
        [camIds] = findClosestFrame(timestamp, timestamp1, timestamp2, timestamp3, timestamp4);
        [camIds_host] = findClosestFrame(timestamp_host, timestamp1, timestamp2, timestamp3, timestamp4);
        imgId = temp.frame_id;
        host_cid = temp.host_cid;
        target_cid = temp.target_cid;
        host_level = temp.host_level;
        target_search_level = temp.target_search_level;
        
        dirCam_host = dirCams{host_cid};
        dirCam_target = dirCams{target_cid};
        pt3d = temp.xyz_w';
        
        if (sum(abs(temp.px_opt - temp.host)) == 0)
            
            slhf = 1;
          
        end
        
        if temp.warp == 0
            draw2 = 0;
        end
        
        
        xyzCell{temp.pid} = [xyzCell{temp.pid}; [pt3d temp.frame_id temp.host_cid temp.target_cid temp.host temp.px_opt, temp.pxA temp.pxB temp.host_cid temp.target_cid]];
        
        img0 = imread(fullfile(inputDir, camInfo(host_cid).name, 'images',dirCam_host(camIds_host(host_cid)).name));
        img3 = imread(fullfile(inputDir, camInfo(target_cid).name, 'images',dirCam_target(camIds(target_cid)).name));
        host = [host;temp.host ];
        target = [target;temp.px_opt ];
        
        
        depth_bounds{temp.pid,1} = [depth_bounds{temp.pid}; [temp.z 1/temp.z_fused pt3d temp.px_opt temp.pxA temp.pxB temp.ab  ]];
        
        %         imgPair{}
%         host_points{host_cid, camIds_host(host_cid)} = [host_points{host_cid, camIds_host(host_cid)}; [temp.host]];
%         target_points{target_cid, camIds(target_cid)} = [target_points{target_cid, camIds(target_cid)}; [temp.px_opt]];
        track_stack = [track_stack; {[host_cid, camIds_host(host_cid) temp.host ],[target_cid, camIds(target_cid) temp.px_opt ], [temp]}];
%         track_stack = [track_stack; [{host_cid, camIds_host(host_cid) temp.host }],[{target_cid, camIds(target_cid) temp.px_opt }]];
        warp_mat_row = size(temp.warp,1);
        zMat(temp.pid, temp.batch_counter) = temp.z;
        if draw
            figure(1),subplot(2,2,1);cla;imshow(img0);hold on;plot(temp.host(1), temp.host(2),'or', 'MarkerSize',3,'LineWidth',3);
          try
            line([temp.warp(1,1) temp.warp(1,warp_mat_row)],[temp.warp(1,warp_mat_row+1) temp.warp(1,2*warp_mat_row)]);
            line([temp.warp(1,1) temp.warp(warp_mat_row,1)],[temp.warp(1,warp_mat_row+1) temp.warp(warp_mat_row,warp_mat_row+1)]);
            line([temp.warp(warp_mat_row,warp_mat_row) temp.warp(1,warp_mat_row)],[temp.warp(warp_mat_row,2*warp_mat_row) temp.warp(1,2*warp_mat_row)]);
            line([temp.warp(warp_mat_row,warp_mat_row) temp.warp(warp_mat_row,1)],[temp.warp(warp_mat_row,2*warp_mat_row) temp.warp(warp_mat_row,warp_mat_row+1)]);
          catch
              sdglkh=1;
          end
            title(sprintf('depth reproj err: %f, det: %f', temp.reproj_error(1), temp.reproj_error(2)));
            figure(1),subplot(2,2,2);cla;imshow(img3);hold on;plot(temp.pxA(1), temp.pxA(2),'og', 'MarkerSize',3,'LineWidth',3);
            plot(temp.pxB(1), temp.pxB(2),'ob', 'MarkerSize',3,'LineWidth',3);
            figure(1),subplot(2,2,[3 4]);cla;
        end
        z_buf = [];
        
        if 1
            
%             figure(1),subplot(2,2,1);imshow(img0);hold on;subplot(2,2,2);
            for step = 1 : size(temp.px_step,1)
                if draw
                    figure(1),subplot(2,2,2);
                    %plot(temp.warp(:,1:10),temp.warp(:,11:20),'.c');
                    plot(temp.px_step(step,1), temp.px_step(step,2),'oy','MarkerSize',1,'LineWidth',1');
                    
                    if(step == size(temp.px_step,1)  - 0 )
                        opt_px = temp.px_step(step,1:2);
                        try
                        line([opt_px(1)-warp_mat_row/2 opt_px(1)+warp_mat_row/2],[opt_px(2)-warp_mat_row/2 opt_px(2)-warp_mat_row/2]);
                        line([opt_px(1)-warp_mat_row/2 opt_px(1)+warp_mat_row/2],[opt_px(2)+warp_mat_row/2 opt_px(2)+warp_mat_row/2]);
                        line([opt_px(1)-warp_mat_row/2 opt_px(1)-warp_mat_row/2],[opt_px(2)-warp_mat_row/2 opt_px(2)+warp_mat_row/2]);
                        line([opt_px(1)+warp_mat_row/2 opt_px(1)+warp_mat_row/2],[opt_px(2)-warp_mat_row/2 opt_px(2)+warp_mat_row/2]);
                        catch
                            sgkhj = 1;
                        end
                    end
                    
                    if(step == size(temp.px_step,1))
                        assert(temp.px_step(step,1) == temp.px_opt(1));
                        subplot(2,2,2);plot(temp.px_step(step,1), temp.px_step(step,2),'om','MarkerSize',2,'LineWidth',2');
                        
                    end
                    title(sprintf('step id: [%d / %d]', step,  size(temp.px_step,1)));
                    drawnow;
                end
                %             plot(featMatBX(feat_id,imgid), featMatBY(feat_id,imgid),'ob','MarkerSize',3,'LineWidth',3');
                %             plot(featMatX(feat_id,imgid), featMatY(feat_id,imgid),'or','MarkerSize',3,'LineWidth',3');title(sprintf('trial id: %d\n', imgid));
                %             figure(1),subplot(2,2,[3 4]);plot(z_buf);drawnow;
                ashk = 1;
                
            end
        end
        
    end
    if 0
        figure(base_fig + pair_num),subplot(2,2,1);imshow(img0);hold on;plot(host(:,1), host(:,2),'or', 'MarkerSize',3,'LineWidth',3);
        subplot(2,2,2);imshow(img3);hold on;plot(target(:,1), target(:,2),'or', 'MarkerSize',3,'LineWidth',3);drawnow;
        pair_num = pair_num + 1;
    end
    
    host_data_ = cell2mat(track_stack(:,1));
    temp_data_ = track_stack(:,3);
    target_data_ = cell2mat(track_stack(:,2));
    host_data = host_data_(target_data_(1,2) - host_data_(:,2)< 11000000,:);
    target_data = target_data_(target_data_(1,2) - host_data_(:,2)< 1100000,:);
    
    host_comb = intersect(host_data(:,1:2), host_data(:,1:2), 'rows');
%     host_comb = host_comb_(target_comb(1,2) - host_comb_(:,2)< 15,:);
%     host_data = host_data();
    target_comb = intersect(target_data(:,1:2), target_data(:,1:2), 'rows');
    
    target_comb_bak = target_comb;
    fig_size = size(host_comb,1) + size(target_comb,1);
    
    fig_row =size(host_comb,1); floor(sqrt(fig_size));
    fig_col =size(target_comb,1)+1; ceil(fig_size/fig_row);
    
    
    if big_fig
        figure(base_fig + pair_num);
    end
    counter = 1;
    for j = 1 : size(host_comb,1)
        dirCam_host_use = dirCams{host_comb(j,1)};
        img_host = imread(fullfile(inputDir, camInfo(host_comb(j,1)).name, 'images',dirCam_host_use(host_comb(j,2)).name));
        
%         dirCam_host_use = dirCams{host_comb(j,1)};
%         img_host = imread(fullfile(inputDir, camInfo(host_comb(j,1)).name, 'images',dirCam_host_use(host_comb(j,2)).name));
        idx_host = find(host_data(:,1) == host_comb(j,1) & host_data(:,2) == host_comb(j,2) );
        if big_fig
        subplot(fig_row,fig_col,1+ fig_col*(j-1));imshow(img_host);hold on;plot(host_data(idx_host,3), host_data(idx_host,4),'or','MarkerSize',2,'LineWidth',2');title(sprintf('[host]: cam %d, frame %d',host_comb(j,1)-1, host_comb(j,2) ));
        else
            if check_cam_id < 0
                figure(base_fig + pair_num);
                subplot(1,fig_col,1);imshow(img_host);hold on;plot(host_data(idx_host,3), host_data(idx_host,4),'or','MarkerSize',2,'LineWidth',2');title(sprintf('[host]: cam %d, frame %d',host_comb(j,1)-1, host_comb(j,2) ));
            else
               if host_comb(j, 1) == check_cam_id
                   figure(base_fig + pair_num);
                subplot(1,fig_col,1);imshow(img_host);hold on;plot(host_data(idx_host,3), host_data(idx_host,4),'or','MarkerSize',2,'LineWidth',2');title(sprintf('[host]: cam %d, frame %d',host_comb(j,1)-1, host_comb(j,2) ));
               end
            end
        end
        %counter = counter + 1;
        
        if (host_comb(j,1) == 1)
            target_order = [3 1 2 0]+1;
        elseif (host_comb(j,1) == 2)
            target_order = [2 0 3 1]+1;
        elseif (host_comb(j,1) == 4)
            target_order = [1 3 0 2]+1;
        else
            target_order = [0 2 1 3]+1;
        end
        
        inde = [];
        for oo = 1 : length(target_order)
            [min_val, min_idx] = min(abs(target_comb_bak(:,1) - target_order(oo)));
            if(min_val == 0)
                inde = [inde; min_idx];
            else
                sadfhk = 1;
            end
        end
        
        target_comb = target_comb_bak(inde,:);
        for k = 1:size(target_comb,1)
            dirCam_target_use = dirCams{target_comb(k,1)};
            img_target = imread(fullfile(inputDir, camInfo(target_comb(k,1)).name, 'images',dirCam_target_use(target_comb(k,2)).name));
            idx = find(host_data(:,1) == host_comb(j,1) & host_data(:,2) == host_comb(j,2) & target_data(:,1) == target_comb(k,1) & target_data(:,2) == target_comb(k,2) );
            if big_fig
                subplot(fig_row,fig_col,1+ k + fig_col*(j-1));imshow(img_target);hold on;plot(target_data(idx,3), target_data(idx,4),'or','MarkerSize',2,'LineWidth',2');title(sprintf('[target]: cam %d, frame %d',target_comb(k,1)-1, target_comb(k,2) ));
            else
                if check_cam_id < 0
                    subplot(1,fig_col,1+ k);imshow(img_target);hold on;plot(target_data(idx,3), target_data(idx,4),'or','MarkerSize',2,'LineWidth',2');title(sprintf('[target]: cam %d, frame %d',target_comb(k,1)-1, target_comb(k,2) ));
                else
                    if host_comb(j, 1) == check_cam_id
                        subplot(1,fig_col,1+ k);imshow(img_target);hold on;plot(target_data(idx,3), target_data(idx,4),'or','MarkerSize',2,'LineWidth',2');title(sprintf('[target]: cam %d, frame %d',target_comb(k,1)-1, target_comb(k,2) ));
                    end
                end
            end
            drawnow;
            %             subplot(fig_row,fig_row,1);imshow(img_host);hold on;plot(host_data(idx,3), host_data(idx,4),'or','MarkerSize',2,'LineWidth',2');title(sprintf('cam %d, frame %d',host_comb(j,1)-1, host_comb(j,2) ));
            %             subplot(1,2,2);imshow(img_target);hold on;plot(target_data(idx,3), target_data(idx,4),'or','MarkerSize',2,'LineWidth',2');title(sprintf('cam %d, frame %d',target_comb(k,1)-1, target_comb(k,2) ));
            for m =  idx'
                
                temp_ = temp_data_{m};
                
                
                if (temp_.pid == check_pid)
                    draw2 = 1;
                    check_pid_count = check_pid_count+1;
                    prograss = [check_pid_count length(inds_)]
                else
                    draw2 = draw3;
                end
                
                timestamp_ = temp_.ts;
                timestamp_host_ = temp_.host_ts;
                [camIds_] = findClosestFrame(timestamp_, timestamp1, timestamp2, timestamp3, timestamp4);
                [camIds_host_] = findClosestFrame(timestamp_host_, timestamp1, timestamp2, timestamp3, timestamp4);
                
                
                dirCam_host_ = dirCams{temp_.host_cid};
                dirCam_target_ = dirCams{temp_.target_cid};
                
                
                
                img00 = imread(fullfile(inputDir, camInfo(temp_.host_cid).name, 'images',dirCam_host_(camIds_host_(temp_.host_cid)).name));
                img33 = imread(fullfile(inputDir, camInfo(temp_.target_cid).name, 'images',dirCam_target_(camIds_(temp_.target_cid)).name));
                if draw2
                    figure(1000),subplot(2,2,1);cla;imshow(img00);hold on;plot(temp_.host(1), temp_.host(2),'or', 'MarkerSize',3,'LineWidth',3);
                    plot(temp_.initial_guess(1), temp_.initial_guess(2),'og', 'MarkerSize',3,'LineWidth',3);
                    try
                        line([temp_.warp(1,1) temp_.warp(1,warp_mat_row)],[temp_.warp(1,warp_mat_row+1) temp_.warp(1,2*warp_mat_row)], 'Color',[0 1 0]);
                        line([temp_.warp(1,1) temp_.warp(warp_mat_row,1)],[temp_.warp(1,warp_mat_row+1) temp_.warp(warp_mat_row,warp_mat_row+1)], 'Color',[0 0 1]);
                        line([temp_.warp(warp_mat_row,warp_mat_row) temp_.warp(1,warp_mat_row)],[temp_.warp(warp_mat_row,2*warp_mat_row) temp_.warp(1,2*warp_mat_row)]);
                        line([temp_.warp(warp_mat_row,warp_mat_row) temp_.warp(warp_mat_row,1)],[temp_.warp(warp_mat_row,2*warp_mat_row) temp_.warp(warp_mat_row,warp_mat_row+1)]);
                    catch
                        sgkhjf = 1;
                    end
                    title(sprintf('host cam id: %d, frame id: %d, depth reproj err: %f, det: %f\n rejected: %d, res: %d, [fail search] count: [%d / %d]',temp_.host_cid-1,camIds_host_(temp_.host_cid), temp_.reproj_error(1), temp_.reproj_error(2), temp_.rejected, temp_.found, temp_.fail_count, temp_.search_count));
                    figure(1000),subplot(2,2,2);cla;imshow(img33);hold on;plot(temp_.pxA(1), temp_.pxA(2),'og', 'MarkerSize',3,'LineWidth',3);
                    plot(temp_.pxB(1), temp_.pxB(2),'ob', 'MarkerSize',3,'LineWidth',3);
                    figure(1000),subplot(2,2,[3 4]);cla;
                end
                if draw2
                   figure(1000),subplot(2,2,2);%title(sprintf('search failed, affine: [%0.3f %0.3f], angle: [%0.3f %0.3f]',  temp_.affine(1), temp_.affine(2), acosd(temp_.angle(1)), acosd(temp_.angle(2))));
                   title(sprintf('search failed, target cam id: %d, frame id: %d\n affine: [%0.3f %0.3f], angle: [%0.3f %0.3f %0.3f]',temp_.target_cid-1,camIds_(temp_.target_cid), temp_.affine(1), temp_.affine(2), acosd(temp_.angle(1)), acosd(temp_.angle(2)), temp_.angle(3)));
                end
                for step = 1 : size(temp_.px_step,1)
                    if draw2
                        figure(1000),subplot(2,2,2);
                        %plot(temp.warp(:,1:10),temp.warp(:,11:20),'.c');
                        plot(temp_.px_step(step,1), temp_.px_step(step,2),'oy','MarkerSize',2,'LineWidth',2');
                        
                        if(step == size(temp_.px_step,1)  - 0 )
                            opt_px = temp_.px_step(step,1:2);
                            try
                                line([opt_px(1)-warp_mat_row/2 opt_px(1)+warp_mat_row/2],[opt_px(2)-warp_mat_row/2 opt_px(2)-warp_mat_row/2], 'Color',[0 1 0]);
                                line([opt_px(1)-warp_mat_row/2 opt_px(1)-warp_mat_row/2],[opt_px(2)-warp_mat_row/2 opt_px(2)+warp_mat_row/2], 'Color',[0 0 1]);
                                line([opt_px(1)-warp_mat_row/2 opt_px(1)+warp_mat_row/2],[opt_px(2)+warp_mat_row/2 opt_px(2)+warp_mat_row/2]);
                                line([opt_px(1)+warp_mat_row/2 opt_px(1)+warp_mat_row/2],[opt_px(2)-warp_mat_row/2 opt_px(2)+warp_mat_row/2]);
                            catch
                                sdgku = 12;
                            end
                        end
                        
                        if(step == size(temp_.px_step,1))
                            assert(temp_.px_step(step,1) == temp_.px_opt(1));
                            subplot(2,2,2);plot(temp_.px_step(step,1), temp_.px_step(step,2),'om','MarkerSize',3,'LineWidth',3');
                            
                        end
                        title(sprintf('target cam id: %d, frame id: %d, step id: [%d / %d]\n affine: [%0.3f %0.3f], angle: [%0.3f %0.3f %0.3f]',temp_.target_cid-1,camIds_(temp_.target_cid), step,  size(temp_.px_step,1), temp_.affine(1), temp_.affine(2), acosd(temp_.angle(1)), acosd(temp_.angle(2)), temp_.angle(3)));
                        drawnow;
                    end
                    %             plot(featMatBX(feat_id,imgid), featMatBY(feat_id,imgid),'ob','MarkerSize',3,'LineWidth',3');
                    %             plot(featMatX(feat_id,imgid), featMatY(feat_id,imgid),'or','MarkerSize',3,'LineWidth',3');title(sprintf('trial id: %d\n', imgid));
                    %             figure(1),subplot(2,2,[3 4]);plot(z_buf);drawnow;
                    ashk = 1;
                    
                end
            end
            
            counter = counter + 1;
        end
        if ~big_fig
            if check_cam_id < 0
                pair_num = pair_num+1;
            else
                if host_comb(j, 1) == check_cam_id
                     pair_num = pair_num+1;
                end
            end
        end
    end
    if big_fig
        pair_num = pair_num+1;
    end
    if ~big_fig
%         close all;
    end
end

pointCloud  = [];
figure(2000);
draw3 = 1;0;
min_trace_len = 5;
for aa = 1 : size(depth_bounds,1)
    if(size(depth_bounds{aa,1},1) > min_trace_len)
        pointCloud = [pointCloud; depth_bounds{aa,1}(end, 3:5)];
        if (draw3)
            aaa = depth_bounds{aa};
            figure(2000),cla;plot(aaa(:,1:2));legend('triangulate','merged');
        end
    end
end

[~, depths] = NormalizeVector(pointCloud);
figure,pcshow(pointCloud(depths<2,:),'MarkerSize', 100);



test = 111; figure,subplot(2,1,1);imshow(zeros(480, 640));hold on;plot(depth_bounds{test}(:,6), depth_bounds{test}(:,7),'or');plot(depth_bounds{test}(:,8), depth_bounds{test}(:,9),'og');plot(depth_bounds{test}(:,10), depth_bounds{test}(:,11),'ob');subplot(2,1,2);plot(depth_bounds{test}(:,2))
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
function depth_filter_data = extractData(Data)



extr_data = 2;
delta = 10+8;

xyMatSize = 0;
offset  = 1;
depth_filter_data = cell(length(Data)-offset,1);
for i = 1+offset : length(Data)
    
   a = Data{i,1}; 

   data.ts = a(1);
   data.host_ts = a(2);
   data.frame_id = a(3)+1;
   data.host_cid = a(4)+1;
   data.target_cid = a(5)+1;
   data.found = a(6);
   data.host_level = a(7);
   data.target_search_level = a(8);
   data.ab = a(9:10)';
   data.initial_guess = a(11:12)'+1;
   
   data.angle = a(13:15);
   data.affine = a(16:17)';
   
   data.fail_count = a(18);
   data.search_count = a(19);
   data.rejected = a(20);
   
   data.batch_counter = a(3+delta);
   data.batch_id = a(4+delta);
   data.pid = a(5+delta)+1;
   data.host = a(6+delta:7+delta)'+1;
   data.pxA = a(8+delta:9+delta)'+1;
   data.pxB = a(10+delta:11+delta)'+1;
   data.px_opt = a(12+delta:13+delta)'+1;
   data.z = a(14+delta);
   data.z_fused = a(15+delta);
   data.xyz_w = a(16+delta:18+delta);
   if 1 % 有patch
       data.px_step = reshape(a(19+delta:end-200-extr_data),3,[])';
       data.warp = reshape(a(end-200-extr_data+1:end-extr_data),20,10)';
       data.px_step(:,1:2) = data.px_step(:,1:2) + 1;
       data.warp = data.warp+1;
   elseif 0
       data.px_step = reshape(a(19+delta:end-72-extr_data),3,[])';
       data.warp = reshape(a(end-72-extr_data+1:end-extr_data),12,6)';
       data.px_step(:,1:2) = data.px_step(:,1:2) + 1;
       data.warp = data.warp+1;
   else
       data.px_step = reshape(a(19+delta:end-xyMatSize-extr_data),3,[])';
       %        data.warp = reshape(a(end-72-extr_data+1:end-extr_data),12,6)';
       data.px_step(:,1:2) = data.px_step(:,1:2) + 1;
       data.warp = 0;
       %    data.warp = data.warp+1;
   end
   %    data.px_step(:,1:2) = data.px_step(:,1:2) + 1;
   %    data.warp = data.warp+1;
   data.reproj_error = a(end-extr_data+1:end);
   depth_filter_data{i-offset,1} = data;
end

end
function [Data] = readNames(fileName)


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
Data{cnt,1} = data;
cnt = cnt + 1;
end
fclose(fid);





end