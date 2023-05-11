function ProfilingChassis()
close all
inputDir = 'D:\Temp\20201123\dump1';
radL = 0.1295; radR = 0.1295; wheelBase = 0.318;
inputDir = 'D:\Temp\20201123\dump2';
inputDir = 'D:\Temp\20201124\dump2';

% inputDir = 'D:\Temp\20201124\dump1';
radL = 0.1285; radR = 0.1285; wheelBase = 0.330;
N = 2000;
cntRatio = 27;

dirInfo = dir(fullfile(inputDir, 'log*.txt'));

goldenVec = [1 0; 2 0; 3 0; 5 0; 0 1000; 0 180; 0 270; 0 90; 0 -1000];

goldenVec = [0 1000; 0 -1000; 0 -2000];
% goldenVec = [0 -1000];

goldenVec = load(fullfile(inputDir, 'config.txt'));
goldenVec0 = goldenVec;
goldenVec(:,2) = deg2rad(goldenVec(:,2));



dt = 0.01;

for i = 1 : length(dirInfo)
    
    log = load(fullfile(inputDir, dirInfo(i).name));
    log(:,6) = log(:,6)./cntRatio;
    log(:,7) = log(:,7)./cntRatio;
    %     id = find(log(:,6) == 0 & log(:,7) == 0);
    id = find(log(:,6) == 0);
    [indexCell] = splitIndex2(id);
    log = log(indexCell{1}(end):indexCell{end}(1)-1,:);
    
    
    v = log(:,4);
    w = log(:,5);
    theta = cumsum([0;w.*dt]);
    xt = cumsum(-v.*sin(theta(1:end-1)).*dt);
    zt = cumsum(v.*cos(theta(1:end-1)).*dt);
    
    ang = [rad2deg(theta)];
    xz = [xt zt];
    
    if goldenVec(i,1) ~= 0
        travel = goldenVec(i,1);
        travel = norm(xz(end,:));
        expL = travel*N/pi/radL;
        expR = travel*N/pi/radR;
        expErr = travel - goldenVec(i,1);
    else
        angle = abs(goldenVec(i,2));
        angle = abs(deg2rad(ang(end)));
        expL = wheelBase*angle*N/2/pi/radL;
        expR = wheelBase*angle*N/2/pi/radR;
        expErr = rad2deg(angle - abs(goldenVec(i,2)));
    end
    figure,plot(cumsum(abs(log(:,6))) - cumsum(abs(log(:,7))));legend('accum cnt L - accum cnt R');
    errL = cumsum(abs(log(:,6))) - expL;
    errR = cumsum(abs(log(:,7))) - expR;
    title(sprintf('%d\naccumCntL - expCntL: %0.1f\naccumCntR - expCntR: %0.1f\nintegratedData - expData: %0.4f (m / deg)',max(abs(goldenVec0(i,:))), errL(end), errR(end), expErr));
    
end


end
function [indexCell,idx] = splitIndex2(idxx)

idx = [1 find(diff(idxx') ~= 1)+1  numel(idxx')+1];
indexCell = mat2cell(idxx', 1, diff(idx))';     

end