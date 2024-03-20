function PlotThresholds()
inputDir = 'G:\matlab\data\direct\gt\D2_011\4';
a = load(fullfile(inputDir, 'thresholds.txt'));

%% cur_time | pid | level | iter | target_cid | is_timestamp_set | last_seen_timestamp | disparity6 | disparity_time6 | disparity_thr | disparity_actual | angle6 | angle_time6 | angle_thr | angle_actual | trace6 | trace_time6 | trace_thr | trace_actual | trace_point2x6

Data = {};


pids = unique(a(:,2));
timestamps = unique(a(:,1));
for i = 1 : length(pids)
    for j = 1 : 4
        for k = 1 : 12
            Data{i, j, k} = cell(length(timestamps), 1);
        end
    end
end
for i = 1 : size(a, 1)
    %Data{i,1}
    status = reshape(a(i,8:35), [], 2);
    trace = a(i,36:61);
    
    data.trace_inlier_flag = a(i,62);
    data.disp_inlier_flag = a(i,63);
    data.angle_inlier_flag = a(i,64);
    
    data.trace_2d = reshape(a(i,65:76), 6,[])';
    data.trace_time_2d = a(i,77:82);
    data.trace_predict_xy = a(i,83:84);
    
    
    data.cur_time = a(i,1);
    data.pid = a(i,2);
    data.level = a(i,3) + 1;
    data.iter = a(i,4) + 1;
    data.target_cid = a(i,5) + 1;
    data.is_timestamp_set = a(i,6);
    data.last_seen_timestamp = a(i,7);
    
    pid_index = find(pids == data.pid);
    time_index = find(timestamps == data.cur_time);
    
    data.disp_thr = status(13,1);
    data.disp_calc = status(14,1);
    data.disp = status(1:6,1);
    data.disp_time = status(7:12,1);
    
    data.angle = status(1:6,2);
    data.angle_time = status(7:12,2);
    data.angle_thr = status(13,2);
    data.angle_calc = status(14,2);
    
    data.trace_lens = trace(1:6);
    data.trace_times = trace(7:12);
    data.trace_thr = trace(13);
    data.trace_calc = trace(14);
    data.trace_points = reshape(trace(15:26), 6, 2)';
    
    stack_id = (data.level - 1) * 4 + data.iter;
    
    Data{pid_index,  data.target_cid, stack_id}{time_index} = data;
end



for pid = 1 : size(Data,1)
    data = Data(pid,:,:);
    for cid = 1 : size(data,2)
        data_time_stack = data(1,cid,:);
        traces_2d_stack = {};
        for stack_id = 1 : size(data_time_stack,3)
            data_time = data_time_stack{1,1,stack_id};
            disp_angle = [];
            traces = [];
            traces_2d = [];
            traces_predict_2d = [];
            trace_flag = [];
            disp_flag = [];
            angle_flag = [];
            level = [];
            cnt = 1;
            for k = 1 : size(data_time,1)
                if (~isempty(data_time{k,1}))
                    if (true || (data_time{k,1}.iter == 3 && data_time{k,1}.level == 0))
                        %                 if (1 && data_time{k,1}.iter == 3 && data_time{k,1}.level == 0)
                        disp_angle(cnt,:) = [data_time{k,1}.cur_time data_time{k,1}.disp_thr data_time{k,1}.angle_thr data_time{k,1}.disp_calc data_time{k,1}.angle_calc data_time{k,1}.trace_thr data_time{k,1}.trace_calc data_time{k,1}.is_timestamp_set];
                        traces(cnt, :) = data_time{k,1}.trace_points(:,5)';
                        traces(cnt, :) = data_time{k,1}.trace_points(:,3)';
                        traces_2d(cnt, :) = data_time{k,1}.trace_2d(:,3)';
                        traces_predict_2d(cnt, :) = data_time{k,1}.trace_predict_xy;
                        trace_flag(cnt,1) = data_time{k,1}.trace_inlier_flag;
                        disp_flag(cnt,1) = data_time{k,1}.disp_inlier_flag;
                        angle_flag(cnt,1) = data_time{k,1}.angle_inlier_flag;
                        level(cnt,1) = data_time{k,1}.level;
                        cnt = cnt + 1;
                    end
                end
            end
            if cnt > 2
                figure(stack_id),clf;subplot(3,2,1); %legend('thr', 'calc');title('disp');
                hold on;plot(find(disp_flag), disp_angle(disp_flag==1, 4),'og');plot(find(~disp_flag), disp_angle(disp_flag == 0, 4),'or');plot([disp_angle(:,[2 4])]);plot([disp_angle(:,[2])] - disp_angle(:,[4]));
                legend('inlier','outlier', 'thr', 'calc', 'diff');title('disp');
                subplot(3,2,3);%legend('thr', 'calc');title('angle');
                hold on;plot(find(angle_flag), disp_angle(angle_flag==1, 5),'og');plot(find(~angle_flag), disp_angle(angle_flag == 0, 5),'or');plot([disp_angle(:,[3 5])]);plot([disp_angle(:,[3])] - disp_angle(:,[5]));
                legend('inlier','outlier', 'thr', 'calc', 'diff');title('angle');
                subplot(3,2,5);%legend('thr', 'calc');title('trace');
                hold on;plot(find(trace_flag), disp_angle(trace_flag==1, 7),'og');plot(find(~trace_flag), disp_angle(trace_flag == 0, 7),'or');plot([disp_angle(:,[6 7])]);plot([disp_angle(:,[6])] - disp_angle(:,[7]));
                legend('inlier','outlier', 'thr', 'calc','diff');title('trace');
                %             hold on;plot(find(disp_angle(:,8)), disp_angle(disp_angle(:,8)==1, 7),'og');legend('thr', 'calc','set time');title('trace');
                subplot(3,2,[2]);imshow(ones(480, 640));hold on;
                plot(traces(logical(trace_flag),1), traces(logical(trace_flag),2),'og');
                plot(traces(~logical(trace_flag),1), traces(~logical(trace_flag),2),'or');
                if 0
                    plot(traces_2d(logical(trace_flag),1), traces_2d(logical(trace_flag),2),'oc');
                    plot(traces_2d(~logical(trace_flag),1), traces_2d(~logical(trace_flag),2),'om');
                end
                
                
                Angle = [];
                for m = 2 : size(traces,1)
                    pt1 = traces(m-1,:);
                    pt2 = traces(m,:);
                    dir = pt1 - pt2;
                    angle = rad2deg(atan2(dir(2), dir(1)));
                    angle = acosd(dot(traces(1,:)./norm(traces(1,:)), pt2./norm(pt2)));
                    Angle = [Angle;  angle];
                    
                    line([pt1(1),pt2(1)]',[pt1(2),pt2(2)]','Color', [0 0 1]);
                    if 0
                        pt3 = traces_2d(m-1,:);
                        pt4 = traces_2d(m,:);
                        line([pt3(1),pt4(1)]',[pt3(2),pt4(2)]','Color', [0 1 1]);
                    end
                end
                %             Angle(Angle < 0) = Angle(Angle < 0) + 360;
                %             Angle(Angle < 0) = 180 + Angle(Angle < 0);
                subplot(3,2,[5]); hold on; plot(Angle);
                subplot(3,2,[4]);plot([(disp_angle(:,[2]) - disp_angle(:,[4])) (disp_angle(:,[3]) - disp_angle(:,[5])) (disp_angle(:,[6]) - disp_angle(:,[7]))]);legend('disp','angle','trace');
                subplot(3,2,[6]);hist([(disp_angle(:,[2]) - disp_angle(:,[4])) (disp_angle(:,[3]) - disp_angle(:,[5])) (disp_angle(:,[6]) - disp_angle(:,[7]))], 50);legend('disp','angle','trace');
                if 1
                    figure(100 + stack_id),clf;subplot(2,2,1);plot([traces_2d(traces_predict_2d(:,1) < 1000,1) traces_predict_2d(traces_predict_2d(:,1) < 1000,1)]);title('x');
                    subplot(2,2,2);plot([traces_2d(traces_predict_2d(:,1) < 1000,2) traces_predict_2d(traces_predict_2d(:,1) < 1000,2)]);title('y');
                    subplot(2,2,3);plot([traces_2d(traces_predict_2d(:,1) < 1000,1) - traces_predict_2d(traces_predict_2d(:,1) < 1000,1)]);title('x');
                    subplot(2,2,4);plot([traces_2d(traces_predict_2d(:,1) < 1000,2) - traces_predict_2d(traces_predict_2d(:,1) < 1000,2)]);title('y');
                end
            end
            traces_2d_stack{stack_id,1} = traces_2d;
        end
%         trace_diff = [[traces_2d_stack{4,1} - traces_2d_stack{8,1}] [traces_2d_stack{4,1} - traces_2d_stack{12,1}]];
    end
end
end