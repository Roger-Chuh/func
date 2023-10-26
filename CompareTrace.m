function CompareTrace(TraceInfoOld_, TraceInfoNew_)

% TraceInfoNew = showTraceUpperStereo('point_trace_statistics_new_full.txt')
% point_trace = load(filename);
% save(mat_name,'point_trace');
% TraceInfoOld = showTraceUpperStereo('point_trace_statistics_old_full.txt')
% point_trace = load(filename);
% save(mat_name,'point_trace');
% CompareTrace(TraceInfoOld, TraceInfoNew)

time_old = cell2mat(TraceInfoOld_(:,5));
time_new = cell2mat(TraceInfoNew_(:,5));

[comm, id1, id2] = intersect(time_old, time_new);

TraceInfoOld = TraceInfoOld_(id1,:);
TraceInfoNew = TraceInfoNew_(id2,:);



draw = 0;



assert(size(TraceInfoOld,1) == size(TraceInfoNew,1));
rep_err = [];

old_vm_count = 0;
new_vm_count = 0;

old_pid_vec = 0;
new_pid_vec = 0;

for i = 1 : size(TraceInfoOld, 1)
    [pids_old, host_cids_old, host_pix_old] = GetPids(TraceInfoOld{i,3});
    [pids_new, host_cids_new, host_pix_new] = GetPids(TraceInfoNew{i,3});
    old_pid_vec = [old_pid_vec; pids_old];
    new_pid_vec = [new_pid_vec; pids_new];
    
    [a,b,c] = intersect(pids_old, pids_new);
    host_cids = host_cids_old(b);
    host_pixs = host_pix_old(b,:);
    
    img_big = uint8(zeros(480*2, 640*2));
    img_big(1:480,1:640) = TraceInfoOld{i,2}{2};
    img_big(1:480,641:end) = TraceInfoOld{i,2}{3};
    img_big(481:end,1:640) = TraceInfoOld{i,2}{1};
    img_big(481:end,641:end) = TraceInfoOld{i,2}{4};
    
    host_img_big = uint8(zeros(480*2, 640*2));
    host_img_big(1:480,1:640) = TraceInfoOld{i,4}{2};
    host_img_big(1:480,641:end) = TraceInfoOld{i,4}{3};
    host_img_big(481:end,1:640) = TraceInfoOld{i,4}{1};
    host_img_big(481:end,641:end) = TraceInfoOld{i,4}{4};
    
    
    figure(1); subplot(1,2,1);imshow([host_img_big]);hold on;title('host');
    subplot(1,2,2);imshow([img_big]);hold on;title('target');
    count = 0;
    for j = 1 : length(a)
        pid = a(j);
        host_pix = host_pixs(j,:);
        host_cid = host_cids(j,:);
        coord_old = zeros(4,3);
        for camId = 1 : size(TraceInfoOld{i,1},1)
            data = cell2mat(TraceInfoOld{i,1}{camId,1});
            if isempty(data)
               continue; 
            end
            id = find(data(:,3) == pid);
            if(~isempty(id))
                coord_old(camId,1:3) = data(id,[1:2 6]);
                asdgkh = 1;
            end  
        end
        coord_new = zeros(4,3);
        for camId = 1 : size(TraceInfoNew{i,1},1)
            data = cell2mat(TraceInfoNew{i,1}{camId,1});
            if isempty(data)
               continue; 
            end
            id = find(data(:,3) == pid);
            if(~isempty(id))
                coord_new(camId,1:3) = data(id,[1:2 6]);
                asdgkh = 1;
            end  
        end
        if draw %pid == 162
            PlotUV(pid, host_pix, coord_old, host_cid, false);
            PlotUV(pid, host_pix, coord_new, host_cid, true);
        end
        
        
        old_cids = find(coord_old(:,1) > 0);
        new_cids = find(coord_new(:,1) > 0);
        [comm_cid, comm1, comm2] = intersect(old_cids, new_cids);
        
        old_vm_count = old_vm_count + length(old_cids);
        new_vm_count = new_vm_count + length(new_cids);
        
        show = false;
        for jj = 1 : length(comm_cid)
            err = [coord_old(comm_cid(jj),:) - coord_new(comm_cid(jj),:) length(comm_cid)];
            rep_err = [rep_err; err];
            if(norm(err(1:2)) > 0.0001)
               show = true; 
            end
        end
        if show
            PlotUV(pid, host_pix, coord_old, host_cid, false);
            PlotUV(pid, host_pix, coord_new, host_cid, true);
            count = count + 1;
        end
        if count > 5
            close all ;
            figure(1); subplot(1,2,1);imshow([host_img_big]);hold on;title('host');
            subplot(1,2,2);imshow([img_big]);hold on;title('target');
            count = 0;
        end
    end
    
    
    
    
end


old_pid_vec1 = unique(old_pid_vec);
new_pid_vec1 = unique(new_pid_vec);

figure,subplot(1,2,1);plot(rep_err(rep_err(:,3) == 1, 1:2));subplot(1,2,2);hist(rep_err(rep_err(:,3) == 1, 1:2), 100);

[~,err_norm] = NormalizeVector(rep_err(:,1:2));

length(find(err_norm > 1))


figure,
end
function [pids, host_cids, host_pix] = GetPids(info)
pids = [];
host_cids = [];
host_pix = [];
for i = 1 : size(info,1)
    temp = cell2mat(info{i,1});
    if isempty(temp)
        continue;
    end
    pids = [pids; unique(temp(:,3))];
    host_cids = [host_cids; i*ones(length(unique(temp(:,3))),1)];
    host_pix = [host_pix; temp(:,1:2)];
end


assert(size(host_pix,1) == size(pids,1));



end
function PlotUV(pid, host_pix, coord, host_cid, is_new)
for kk = 1 : size(coord,1)
            if(coord(kk,1) > 0 && coord(kk,2) > 0 || kk == host_cid)
                if(kk == 1)
%                     subplot(1,2,2);title(sprintf('target, level: [%d %d]',))
                    if is_new
                        plot(coord(kk,1) + 0, 480 + coord(kk,2),'*r');
                    else
                        plot(coord(kk,1) + 0, 480 + coord(kk,2),'*b');
                    end
                    text(coord(kk,1)+2 +0,480 + coord(kk,2)+2, num2str(pid),'Color', [0 1 0]);
                   
                        if(kk == host_cid && is_new)
                            subplot(1,2,1);
                            plot(host_pix(:,1), 480 + host_pix(:,2),'*r');
                            text(host_pix(:,1)+2,480 + host_pix(:,2)+2, num2str(pid),'Color', [0 1 0]);
                        end
                   
                elseif kk == 2
                    subplot(1,2,2);
                    if is_new
                        plot(coord(kk,1) + 0, coord(kk,2),'*r');
                    else
                        plot(coord(kk,1) + 0, coord(kk,2),'*b');
                    end
                    text(coord(kk,1)+2 + 0,coord(kk,2)+2, num2str(pid),'Color', [0 1 0]);
                    if(kk == host_cid && is_new)
                        subplot(1,2,1);
                        plot(host_pix(:,1), host_pix(:,2),'*r');
                        text(host_pix(:,1)+2,host_pix(:,2)+2, num2str(pid),'Color', [0 1 0]);
                    end
                elseif kk == 3
                    subplot(1,2,2);
                    if is_new
                        plot(coord(kk,1) + 640, coord(kk,2),'*r');
                    else
                        plot(coord(kk,1) + 640, coord(kk,2),'*b');
                    end
                    text(coord(kk,1)+2 + 640,coord(kk,2)+2, num2str(pid),'Color', [0 1 0]);
                    if(kk == host_cid && is_new)
                        subplot(1,2,1);
                        plot(host_pix(:,1) + 640, host_pix(:,2),'*r');
                        text(host_pix(:,1)+2 + 640,host_pix(:,2)+2, num2str(pid),'Color', [0 1 0]);
                    end
                else
                    subplot(1,2,2);
                    if is_new
                        plot(coord(kk,1) + 640, coord(kk,2) + 480,'*r');
                    else
                       plot(coord(kk,1) + 640, coord(kk,2) + 480,'*b'); 
                    end
                    text(coord(kk,1)+2 + 640,coord(kk,2)+2 + 480, num2str(pid),'Color', [0 1 0]);
                    if(kk == host_cid && is_new)
                        subplot(1,2,1);
                        plot(host_pix(:,1) + 640, host_pix(:,2) + 480,'*r');
                        text(host_pix(:,1)+2 + 640,host_pix(:,2)+2 + 480, num2str(pid),'Color', [0 1 0]);
                    end
                end
                
            end
end
end