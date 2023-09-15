function ConcatImg()

close all;

show_fig = true; false;
save_vid = true; false;



if ~show_fig
    save_vid = false;
end

inputDir = 'G:\matlab\data\direct\gt\D2_006\imgs_show';
total_duration = 120.7;
start_pixel = [237*1 687];
end_pixel = [1581*1 687];


inputDir = 'G:\matlab\data\direct\gt\D2_011\imgs_for_show_good\imgs';
total_duration = 89.7;
start_pixel = [167*1 650];
end_pixel = [1485*1 650];

inputDir = 'G:\matlab\data\direct\gt\D2_011\imgs_for_show_bad\imgs';
total_duration = 89.7;
start_pixel = [159*1 673];
end_pixel = [1503*1 673];

inputDir = 'G:\matlab\data\direct\gt\D2_001\imgs_for_show\imgs';
total_duration = 118.2;
start_pixel = [163*1 669];
end_pixel = [1508*1 669];

inputDir = 'G:\matlab\data\direct\gt\D2_011\hack_011_100\imgs';
total_duration = 89.7;
start_pixel = [158*1 646];
end_pixel = [1475*1 646];

% inputDir = 'G:\matlab\data\direct\gt\D2_011\hack_011_200\imgs';
% total_duration = 89.7;
% start_pixel = [158*1 667];
% end_pixel = [1476*1 667];


% inputDir = 'G:\matlab\data\direct\gt\D2_011\hack_011_2.0_orig\imgs';
% total_duration = 89.7;
% start_pixel = [157*1 668];
% end_pixel = [1502*1 668];


inputDir = 'G:\matlab\data\direct\gt\D2_011\hack_011_2.0_minus_5\imgs';
total_duration = 89.7;
start_pixel = [157*1 668];
end_pixel = [1502*1 668];




inputDir = 'G:\matlab\data\direct\gt\D2_011\img_opt_pose_df_bad_config_good_result\imgs';
total_duration = 89.7;
start_pixel = [162*1 670];
end_pixel = [1479*1 670];



inputDir = 'G:\matlab\data\direct\gt\D2_011\img_opt_pose_df_best_config_best_result\imgs';
total_duration = 89.7;
start_pixel = [164*1 653];
end_pixel = [1480*1 653];



inputDir = 'G:\matlab\data\direct\gt\D2_011\img_opt_pose_no_df_less_good_result_kf_err_increase\imgs';
total_duration = 89.7;
start_pixel = [161*1 647];
end_pixel = [1478*1 647];


inputDir = 'G:\matlab\data\direct\gt\D2_006\imgs_show_coverage\imgs';
total_duration = 120.8;
start_pixel = [151*1 650];
end_pixel = [1497*1 650];



inputDir = 'G:\matlab\data\direct\gt\D2_006\imgs_show_big_fov\imgs';
total_duration = 120.8;
start_pixel = [162*1 656];
end_pixel = [1508*1 656];




inputDir = 'G:\matlab\data\direct\gt\D2_006\imgs_show_big_fov_2\imgs';
total_duration = 120.8;
start_pixel = [151*1 659];
end_pixel = [1497*1 659];



inputDir = 'G:\matlab\data\direct\gt\D2_006\imgs_show_big_fov_3\imgs';
total_duration = 120.8;
start_pixel = [150*1 656];
end_pixel = [1495*1 656];


inputDir = 'G:\matlab\data\direct\gt\D2_006\imgs_show_big_fov_4\imgs';
total_duration = 120.8;
start_pixel = [152*1 647];
end_pixel = [1498*1 647];





inputDir = 'G:\matlab\data\direct\gt\D2_006\imgs_show_big_fov_5\imgs';
total_duration = 120.8;
start_pixel = [154*1 645];
end_pixel = [1499*1 645];


inputDir = 'G:\matlab\data\direct\gt\D2_006\imgs_show_big_fov_6\imgs';
total_duration = 120.8;
start_pixel = [152*1 654];
end_pixel = [1497*1 654];


inputDir = 'G:\matlab\data\direct\gt\D2_006\imgs_show_new_kf_1\imgs';
total_duration = 120.8;
start_pixel = [153*1 657];
end_pixel = [1499*1 657];


inputDir = 'G:\matlab\data\direct\gt\D2_006\imgs_show_new_kf_2\imgs';
total_duration = 120.8;
start_pixel = [151*1 649];
end_pixel = [1496*1 649];


inputDir = 'G:\matlab\data\direct\gt\D2_006\imgs_show_new_kf_3\imgs';
total_duration = 120.8;
start_pixel = [151*1 658];
end_pixel = [1496*1 658];



inputDir = 'G:\matlab\data\direct\gt\D2_006\imgs_show_new_kf_vio_1\imgs';
total_duration = 114.6;
start_pixel = [142*1 658];
end_pixel = [1476*1 658];


inputDir = 'G:\matlab\data\direct\gt\D2_006\imgs_show_new_kf_vo_1\imgs';
total_duration = 120.8;
start_pixel = [151*1 653];
end_pixel = [1497*1 653];


inputDir = 'G:\matlab\data\direct\gt\D2_006\imgs_show_new_kf_vo_2\imgs';
total_duration = 120.8;
start_pixel = [148*1 651];
end_pixel = [1494*1 651];


inputDir = 'G:\matlab\data\direct\gt\D2_006\imgs_show_new_kf_vo_3\imgs';
total_duration = 89.7;
start_pixel = [149*1 660];
end_pixel = [1494*1 660];


inputDir = 'G:\matlab\data\direct\gt\D2_011\imgs_0.25\imgs';
total_duration = 89.7;
start_pixel = [164*1 649];
end_pixel = [1481*1 649];


inputDir = 'G:\matlab\data\direct\gt\D2_004\bad_case\imgs';
total_duration = 115.4;
start_pixel = [149*1 657];
end_pixel = [1466*1 657];










% ------------------------------------------------------------------

dirInfo = dir(fullfile(inputDir,'0*'));

inputDir1 = fullfile(inputDir, dirInfo(1).name);
inputDir2 = fullfile(inputDir, dirInfo(2).name);


dirInfo1 = dir(fullfile(inputDir1, '*.png'));
dirInfo2 = dir(fullfile(inputDir2, '*.png'));



timestamp1 = zeros(length(dirInfo1),1);
for i = 1 : length(dirInfo1)
    idx = find(dirInfo1(i).name == '_');
    duration = str2double(dirInfo1(i).name(idx(1)+1 :idx(2)-1))*1e-9;
    timestamp1(i,1) = duration;
end
timestamp2 = zeros(length(dirInfo2),1);
coverage_ratio = zeros(length(dirInfo2),5);
for i = 1 : length(dirInfo2)
    idx = find(dirInfo2(i).name == '_');
    duration = str2double(dirInfo2(i).name(idx(1)+1 :idx(2)-1))*1e-9;
    rat0 = str2double(dirInfo2(i).name(idx(2)+1 :idx(3)-1))./100;
    rat1 = str2double(dirInfo2(i).name(idx(3)+1 :idx(4)-1))./100;
    rat2 = str2double(dirInfo2(i).name(idx(4)+1 :idx(5)-1))./100;
    rat3 = str2double(dirInfo2(i).name(idx(5)+1 :idx(6)-1))./100;
    rat_mean = str2double(dirInfo2(i).name(idx(6)+1 :end-4))./100;
    timestamp2(i,1) = duration;
    coverage_ratio(i,:) = [rat0 rat1 rat2 rat3 rat_mean];
end

[a,b,c] = intersect(timestamp1, timestamp2);


% assert(length(dirInfo1) == length(dirInfo2));

err_image = imread(fullfile(inputDir, 'error.png'));

% err_image = imresize(err_image, [size(err_image,1), 2*size(err_image,2)]);


pose = load(fullfile(inputDir, 'output.txt'));
Data = pose(:,2:8);
poseMat = [];
for i = 1 : size(Data,1)
    data = Data(i,:);
    xyzw = data(4:7);
    trans = data(1:3);
    
%     R2 = quatern2rotMat(xyzw([4 1 2 3]));
    R = quat2rotm(xyzw([4 1 2 3]));
    poseMat = [poseMat; [reshape(R,1,9), trans]];
%     rotm2quat(R)
%     rotMat2quatern(R2)
end

total_length = norm(start_pixel - end_pixel);

err_image_resized = 0;

if show_fig
    figure(1);
end

is_full_screen_set = false;

% figure(1);subplot(2,1,2);imshow(err_image);hold on;
cur_pix_list = [];
if save_vid
    v = VideoWriter(fullfile(inputDir, 'trace_error.avi'),'Uncompressed AVI');
%     v = VideoWriter(fullfile(inputDir, 'trace_error.avi'), 'Indexed AVI');
    v.FrameRate = 10;
    open(v);
end
close_window_interval = 100;
kf_flag = [];
Duration = [];
is_first_frame_set = false;
for i = 1 : length(a)
    
   idx = find(dirInfo1(b(i)).name == '_');
   is_kf = str2double(dirInfo1(b(i)).name(end-4));
   kf_flag = [kf_flag; is_kf];
   
   duration = str2double(dirInfo1(b(i)).name(idx(2)+1 : idx(3)-1))*1e-9;
   if ~is_first_frame_set
       delta_time = duration;
       is_first_frame_set = true;
   end
   duration = duration - delta_time;
   
   duration_ratio = duration / total_duration;
   
   
   Duration = [Duration; duration];
   img1 = imread(fullfile(inputDir1, dirInfo1(b(i)).name)); 
   img2 = imread(fullfile(inputDir2, dirInfo2(c(i)).name)); 
   coverage = coverage_ratio(c(i),:);
   img = [img1 img2];
   
   if ~err_image_resized
       ratio = size(img, 2)/size(err_image, 2);
       start_pixel(1) = ratio*start_pixel(1);
       end_pixel(1) = ratio*end_pixel(1);
       total_length = norm(start_pixel - end_pixel);
       
       err_image = imresize(err_image, [size(err_image,1),size(img,2)]);
      err_image_resized = 1; 
   end
   cur_pix = start_pixel + [total_length * duration_ratio 0];
   cur_pix_list = [cur_pix_list; [cur_pix is_kf coverage]];
   img = [img; err_image];
   vkf_id = find(cur_pix_list(:,3) == 1);
   tkf_id = find(cur_pix_list(:,3) == 2);
   
   if show_fig
       if( mod(i,close_window_interval) == 0)
           figure(1);
           close ;
           figure(1);
           is_full_screen_set = false;
       end
       
       imshow(img);hold on;
       
       
       
       plot(cur_pix_list(:,1), cur_pix_list(:,2)+size(img1,1),'xr','MarkerSize', 8,'LineWidth',8);
       plot(cur_pix_list(:,1), cur_pix_list(:,2)+size(img1,1)-100*(cur_pix_list(:,8)-1)-100,'.-m','MarkerSize', 1,'LineWidth',1);
       if(~isempty(vkf_id))
           plot(cur_pix_list(vkf_id,1), cur_pix_list(vkf_id,2)+size(img1,1),'xg','MarkerSize', 15,'LineWidth',15);
       end
       if(~isempty(tkf_id))
           plot(cur_pix_list(tkf_id,1), cur_pix_list(tkf_id,2)+size(img1,1),'xb','MarkerSize', 3,'LineWidth',3);
       end
       %    subplot(2,1,[1]);imshow(img);
       %    subplot(2,1,[2]);plot(cur_pix(:,1), cur_pix(:,2),'xr');
       if(~is_full_screen_set)
           set(gcf,'outerposition',get(0,'screensize'));
           is_full_screen_set = true;
       end
       drawnow;
       frame = getframe(gcf);
       image =  imresize(frame.cdata(28:922, 493:1406,:),0.5);
       if save_vid
           writeVideo(v,image);
       end
   end
   
   
end

if save_vid
    close(v);
end

kf_ids = find(kf_flag == 1);

timestamp_used = timestamp1(b,:);
poseMat_kf = [];
for k = 1 : length(kf_ids)
    time_ = timestamp_used(kf_ids(k));
    [min_delta, min_id] = min(abs(time_ - pose(:,1)));
    poseMat_kf = [poseMat_kf; poseMat(min_id,:)];
end


figure,subplot(1,2,1);plotPath(poseMat);plotPath(poseMat_kf); subplot(1,2,2);plot(Duration, coverage_ratio(:,5)); grid on;hold on;plot(Duration(kf_flag==1), coverage_ratio(kf_flag==1,5),'xr');
fprintf(sprintf('vkf num: %d\n',sum(kf_flag == 1)));


delta_pose = {};
for i = 1 : size(poseMat_kf,1)
   pose_i = [reshape(poseMat_kf(i,1:9),3,3) poseMat_kf(i,10:12)';0 0 0 1];
   delta_pose{i,1} = [];
    for j = 1: size(poseMat_kf,1)
        if(i == j)
           % continue; 
        end
           pose_j = [reshape(poseMat_kf(j,1:9),3,3) poseMat_kf(j,10:12)';0 0 0 1];
           delta_pose_ = inv(pose_j) * pose_i;
           delta_pose{i,1} = [delta_pose{i,1}; [rad2deg(norm(rodrigues(delta_pose_(1:3,1:3)))) norm(delta_pose_(1:3,4))]];
    end
    
    
end



end