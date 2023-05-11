function AnalyseDepthFilter3()


inputDir = 'G:\matlab\data\direct\gt\1';
% inputDir = 'G:\matlab\data\direct\gt\self';
% inputDir = 'G:\matlab\data\direct\gt\ke\2';
inputDir = 'G:\matlab\data\direct\gt\ke\3';
data = load(fullfile(inputDir, 'trace_info_statistics.txt'));


ts = unique(data(:,1));
pids = unique(data(:,14));



%%   2    |    3   |   4   |    5      |    6   |  7  |  8  |  9  |  10 |  11  |    12    |    13    |   14     |    15     |     16      |      17        |      18        |      19        |   20      |    21
%% reproj | reproj | sigma | sigma_thr | idepth | cid | fid |  x  |   y |   z  | host_cid | host_fid | seed_id  | trace_ts  |  has_guess  | idepth_before  |  idepth_after  |   initial_pid  |  z_reset  |  new_thr
% for i = 1 : length(ts)
%     
%    ts_ = find(data(:,1) == ts(i)); 
%    data_ = data(find(data(:,1) == ts(i)),:);  
%    pid = unique(data_(:,14));
%    
%    
% end


offset = 0.02;

trace_len = -2;7;10;10; 5;

bad_count = 0;
draw = 0;1;0; 1;0;


pid_check = 206 - 1;



for i = 1 : length(pids)
    pid_ = pids(i);
    
   data_ = data(find(data(:,14) == pids(i)),:); 

   if (pid_check == pid_) % || (data_(end,20) == 1)
       draw = 1;
    else
        draw = 0;
    end
   a = (data_(:,1) - data_(1,1))./10^6;
   b = (data_(:,15) - data_(1,15))./10^6;

   
   id = find(b == 0);
   if length(unique(b)) < trace_len
%        if size(data_,1) < trace_len
       continue;
   end
    
   if(data_(end,5)<data_(end,4))
%       continue;
   end
   
   if(max(data_(:,3)) > 3)
      sdfhkl = 1; 
   end
   
   bad_count = bad_count+1;
   
   if draw
       figure(190),clf;subplot(2,3,1);hold on;plot3(data_(end,9),data_(end,10),data_(end,11),'og','MarkerSize',5,'LineWidth',5);plot3(data_(1,9),data_(1,10),data_(1,11),'ob','MarkerSize',5,'LineWidth',5);
       for(k = 1 : size(data_,1))
           text(data_(k,9)+offset,data_(k,10)+offset, data_(k,11)+offset,num2str(k));
       end
       axis equal;pcshow(data_(:,9:11));plot3(data_(:,9),data_(:,10),data_(:,11),'-xr');legend('final','initial','trace');
       subplot(2,3,2);plot(data_(:,[4 5 21]));legend('sigma', 'sigma thr init', 'sigma thr new'); grid on;
       subplot(2,3,3);plot(data_(:,2:3));legend('reproj single', 'reproj all');title(sprintf('pid: %d, traceLen: %d, has guess: %d\n host fid: %d, initial pid: %0.1f, sigma thr reset: %d', pid_+1, size(data_,1), data_(1,16), data_(1,13), data_(1,19)+1, data_(end,20)));grid on;
       bound1 = data_(:,18)-data_(:,4);
       bound1(bound1 < 0) = 0.05;
       bound1(bound1 < 0.05) = 0.05;
       subplot(2,3,[4 5]);cla;plot([  1./(data_(:,18)+data_(:,4)) 1./bound1 ]);hold on;plot(1./(data_(:,17:18)),'MarkerSize',5,'LineWidth',5);legend('bound','bound','idepth before', 'idepth after');title(sprintf('start id: %d',id(end)));grid on
       %    subplot(2,3,5);plot([data_(:,18)  data_(:,18)+data_(:,4) data_(:,18)-data_(:,4) ]);legend('idepth','bound','bound');
   end
end
[bad_count length(pids)]
pc = load(fullfile(inputDir, 'landmark_statistics.txt'));pc = pc(1:end,:);
idepth = pc(:,2);
landmark = pc(:,3:5);
sigma = pc(:,6);
goodness = pc(:,7);
[~,depths_] = NormalizeVector(landmark);
%  figure(4000),pcshow(landmark(goodness == 1 & depths_<2 & sigma < 1,:),'MarkerSize', 100);
figure,pcshow(landmark(goodness == 1 & depths_<15.5 & sigma < 0.1 & idepth > -0.8,:),'MarkerSize', 100);
end