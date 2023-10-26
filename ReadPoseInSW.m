function ReadPoseInSW()
inputDir = 'G:\matlab\data\direct\gt\D2_007';


inputDir = 'G:\matlab\data\direct\gt\D2_009';
total_duration = 89.6;
start_pixel = [162*1 662];
end_pixel = [1480*1 662];


inputDir = 'G:\matlab\data\direct\gt\D2_009\2';
total_duration = 89.6;
start_pixel = [150*1 603];
end_pixel = [1496*1 603];


inputDir = 'G:\matlab\data\direct\gt\D2_009\3';
total_duration = 89.6;
start_pixel = [153*1 585];
end_pixel = [1499*1 585];


inputDir = 'G:\matlab\data\direct\gt\D2_007\1';
total_duration = 100;
start_pixel = [142*1 647];
end_pixel = [1514*1 647];


inputDir = 'G:\matlab\data\direct\gt\D2_011\1';
total_duration = 90;
start_pixel = [141*1 617];
end_pixel = [1513*1 617];

inputDir = 'G:\matlab\data\direct\gt\D2_011\2';
total_duration = 90;
start_pixel = [141*1 617];
end_pixel = [1513*1 617];

inputDir = 'G:\matlab\data\direct\gt\D2_011\3';
total_duration = 90;
start_pixel = [141*1 617];
end_pixel = [1513*1 617];

inputDir = 'G:\matlab\data\direct\gt\D2_011\4';
total_duration = 90;
start_pixel = [141*1 617];
end_pixel = [1513*1 617];






err_image = imread(fullfile(inputDir, 'error.png'));
%-----------------------------------------------------



%% cur_time  |   sw_fid_iter_time  |   sw_fid_iter_framekind   | xyz xyzw  |   sw_front_time
%%    1      |         2           |             3             |  4 - 10   |       11


delta_fig = 10; 1;0;10;



data = load(fullfile(inputDir, 'output_sw.txt'));

data = data(data(:,3) > 0,:);


% id = find(data(:,3) == 2);
% vkf_timestamp = unique(sdata(id,2));


Data = data(:,4:10);
poseMat = [];
for i = 1 : size(Data,1)
    data1 = Data(i,:);
    xyzw = data1(4:7);
    trans = data1(1:3);
    
    %     R2 = quatern2rotMat(xyzw([4 1 2 3]));
    R = quat2rotm(xyzw([4 1 2 3]));
    poseMat = [poseMat; [reshape(R,1,9), trans]];
    %     rotm2quat(R)
    %     rotMat2quatern(R2)
end




% poseMat = poseMat(data(:,3) > 0,:);
% data = data(data(:,3) > 0,:);


timestamp_orig = data(:,1);

data(:,[1 2 11]) = data(:,[1 2 11]) - data(1,2);

id = find(data(:,3) == 2);
vkf_timestamp = unique(data(id,2));


front_timestamp = unique(data(id,11));


timestamps = unique(data(:,2));


figure(2+delta_fig);clf;plotPath(poseMat, 5);plotPath(poseMat(data(:,3) == 2,:),5);

 figure(4+delta_fig),imshow(err_image);

total_length = norm(start_pixel - end_pixel);
% figure(1),subplot(2,2,[1 2]);imshow(err_image);

Trans_diff = {};


base_offset = 0;

for i = 1 : length(timestamps)
    timestamp = timestamps(i);
    id = find(data(:,2) == timestamp);
    if data(id(1),3) ~= 2
        %        continue;
    end
    
    
    
    
    
    pose_opt = poseMat(id,:);
    duration = data(id,1);
    
    duration_ratio = duration / total_duration;
    cur_pix = start_pixel + [total_length * duration_ratio zeros(length(duration_ratio),1)];
    
    
    trans_diff = 1000.*(pose_opt(:,10:12) - pose_opt(1,10:12));
    base_offset = base_offset - min((trans_diff(:))); 
    Trans_diff{i,1} = [duration trans_diff + base_offset.*ones(size(trans_diff,1),3)];
    Trans_diff{i,2} = data(id(1),2);
    Trans_diff{i,3} = data(id, 11);
    Trans_diff{i,4} = [timestamp_orig(id,:) poseMat(id,:)];
    Trans_diff{i,5} = [length(id)];
    
    base_offset = base_offset + max((trans_diff(:))) + 2; 
    
    scale = 150./max(abs(trans_diff(:)));
    
    if 0
        figure(4),imshow(err_image);hold on;plot(cur_pix(:,1), [repmat(cur_pix(:,2),1,3) -  scale.*trans_diff]);title(sprintf('trans-diff, kf-kind: %d',data(id(1),3)));
        
        figure(1),subplot(2,3,[1 2 3]);cla;imshow(err_image);hold on;plot(cur_pix(:,1), [repmat(cur_pix(:,2),1,3) -  scale.*trans_diff]);title(sprintf('trans-diff, kf-kind: %d',data(id(1),3)));
        subplot(2,3,[4 5]);cla;plotPath(pose_opt);title(sprintf('first appearance: %fs', timestamp));
        subplot(2,3,[6]);plot(duration, trans_diff);title(sprintf('trans-diff(mm), kf-kind: %d',data(id(1),3)));
        %    subplot(2,2,4);plot(duration, trans_diff);title(sprintf('trans-diff, kf-kind: %d',data(id(1),3)));
    end
    
    
%     figure(3);clf;hold on;grid on;
%     for j = 1 : size(Trans_diff,1)
%         plot(Trans_diff{j,1}(:,1), Trans_diff{j,1}(:,2),'-r');
%         plot(Trans_diff{j,1}(:,1), Trans_diff{j,1}(:,3),'-g');
%         plot(Trans_diff{j,1}(:,1), Trans_diff{j,1}(:,4),'-b');
%     end
%     drawnow;
    
end





min_length = 1;%400;  %min(cell2mat(Trans_diff(:,5)));

output_init = [];
output_opt = [];
output_opt_2 = [];


figure(3+delta_fig);clf;hold on;grid on
for j = 1 : size(Trans_diff,1)
    plot(Trans_diff{j,1}(:,1), Trans_diff{j,1}(:,2),'-r');
    plot(Trans_diff{j,1}(:,1), Trans_diff{j,1}(:,3),'-g');
    plot(Trans_diff{j,1}(:,1), Trans_diff{j,1}(:,4),'-b');
    plot(Trans_diff{j,1}(:,1), 5.*Trans_diff{j,3},'-m');
    
    
    if min_length > size( Trans_diff{j,4},1)
        continue;
    end
        
    
    init_pose = Trans_diff{j,4}(1,:);
    opt_pose = [Trans_diff{j,4}(1,1) Trans_diff{j,4}(end,[2:end])];
    
    opt_pose_2 = [Trans_diff{j,4}(1,1) Trans_diff{j,4}(min_length,[2:end])];
    
    
    q_init_wxyz = rotm2quat(reshape(init_pose(2:10),3,3));
    q_opt_wxyz = rotm2quat(reshape(opt_pose(2:10),3,3));
    
    q_opt_wxyz_2 = rotm2quat(reshape(opt_pose_2(2:10),3,3));
    
    output_init = [output_init; [init_pose(1) init_pose(:,11 : 13) q_init_wxyz([2 3 4 1])]];
    output_opt = [output_opt; [opt_pose(1) opt_pose(:,11 : 13) q_opt_wxyz([2 3 4 1])]];
    output_opt_2 = [output_opt_2; [opt_pose_2(1) opt_pose_2(:,11 : 13) q_opt_wxyz_2([2 3 4 1])]];
    if ~isempty(find(vkf_timestamp == Trans_diff{j,2}))
        plot(Trans_diff{j,1}(end,1), Trans_diff{j,1}(end,2),'or');
        plot(Trans_diff{j,1}(end,1), Trans_diff{j,1}(end,3),'og');
        plot(Trans_diff{j,1}(end,1), Trans_diff{j,1}(end,4),'ob');
    end
end

 write2file(output_init, inputDir, 'output', 1);
 write2file(output_opt, inputDir, 'output', 2);
 write2file(output_opt_2, inputDir, 'output', 3);

end

function write2file(data, inputDir, name, id)
f_id=fopen(fullfile(inputDir,sprintf(strcat(name,'_%d.txt'),id)),'w');

for i = 1 : size(data,1)
    fprintf(f_id,'%f %f %f %f %f %f %f %f\n',data(i,1),data(i,2),data(i,3),data(i,4),data(i,5),data(i,6),data(i,7),data(i,8));
    
end
fclose(f_id);
end