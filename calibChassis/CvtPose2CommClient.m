function CvtPose2CommClient(inputDir)

dirInfo = dir(fullfile(inputDir, 'odom*.txt'));




for ii = 1 : length(dirInfo)
    
    pose1 = load(fullfile(inputDir, dirInfo(ii).name));
%     pose1 = pose1(1:3:end,:);
    if size(pose1,2) == 5
        pose1 = pose1(:,2:5);
    end
    
    [~,dist] = NormalizeVector(pose1(:,2:3));
    id = find(dist == 0);
    [indexCellRobo, indexRobo] = splitIndex2(id);
    pose = pose1(max(1,indexCellRobo{1}(end)-1:end),:);
    
    
    
%     pose(:,1) = 1000.*pose(:,1);
    pose(:,1) = pose(:,1) - pose(1,1) + 100;
    %     rmatSlam = [[slampath(7) 0 -slampath(5)]' [0 1 0]' slampath(5:7)']';
    pose = [pose(:,[1 2 2 3 4])];
    pose(:,3) = 0;
    pose0 = pose;
    dltAng1 = diff(pose(:,5));
    a = find(abs(dltAng1) > 2);
    a = a+1;
    for k = 1 : length(a)
        if dltAng1(a(k)-1) < 0
            pose(a(k):end,5) = pose(a(k):end,5) + 2*pi;
        else
            pose(a(k):end,5) = pose(a(k):end,5) - 2*pi;
        end
    end
    
    
    fid = fopen(fullfile(inputDir, 'comm_client_log.txt'),'w');
    for i = 1 : size(pose, 1)
        
        trans = pose(i, 2:4);
        
        rmat = roty(rad2deg(pose(i,5)));
        
        if i == 1
            initT = [roty(-90) [0 0 0]';0 0 0 1]*inv([rmat trans';0 0 0 1]);

        end
        transformed = initT*[rmat trans';0 0 0 1];
        
% %         rvec = rmat(:,3)';
        rvec = transformed(1:3,3)';
        transs = transformed(1:3,4)';
        poseMat(i,:) = [reshape(transformed(1:3,1:3),1,9) transformed(1:3,4)'];
        
        fprintf(fid, sprintf('Robot FID %d %d %d %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f\n',...
            i, round(pose(i,1)), round(pose(i,1)),transs(1),transs(2), transs(3), rvec(1), rvec(2), rvec(3)));
        
        fprintf(fid, sprintf('Slam FID %d %d %d %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f\n',...
            i, round(pose(i,1)), round(pose(i,1)),transs(1),transs(2), transs(3), rvec(1), rvec(2), rvec(3)));
    end
    fclose(fid);
    
end


end