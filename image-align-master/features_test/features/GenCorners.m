function GenCorners()
inputDir = 'G:\papers\direct_methods\image-align-master\features_test\features\save';

dirInfo = dir(fullfile(inputDir, '*.bmp'));

timestamps = [];
for i = 1 : length(dirInfo)
    timestamps = [timestamps; uint64(str2double(dirInfo(i).name(1:end-6)))];
    
end
timestamps = unique(timestamps);
for i = 1 : length(timestamps)
    
    fid1 = fopen(fullfile(inputDir,sprintf('%d.yaml', timestamps(i))),'w');
    fprintf(fid1,'%%YAML:1.0\n');
    fprintf(fid1,'---\n');
    for j = 1 : 4
        img = imread(fullfile(inputDir,sprintf('%d_%d.bmp', timestamps(i), j-1)));
        corners = Detect(img);
        corners = corners - 1;
        Corners{j,1} = corners;
        fprintf(fid1,'Cam_%d: !!opencv-matrix\n   rows: %d\n   cols: %d\n   dt: d\n', j-1, size(corners,1),2);
        for k = 1 : size(corners,1)
           if k == 1
               fprintf(fid1,sprintf('   data: [ %f, %0f, \n',corners(k,1), corners(k,2)));
           elseif k == size(corners,1)
               fprintf(fid1,sprintf('           %f, %0f ]\n',corners(k,1), corners(k,2)));
           else
               fprintf(fid1,sprintf('           %f, %0f, \n',corners(k,1), corners(k,2)));
           end
            
        end
    end
    fclose(fid1);
    
    
    
%     fprintf(fid1,sprintf('       %0.15f, %0.15f, %0.15f, %0.15f, %0.15f ]\n',intrMatL(2,2), intrMatL(2,3),intrMatL(3,1), intrMatL(3,2), intrMatL(3,3)));
%     fprintf(fid1,'D1: !!opencv-matrix\n   rows: 1\n   cols: 5\n   dt: d\n');
%     fprintf(fid1,sprintf('   data: [ %0.15f, %0.15f, \n',kcL(1), kcL(2)));
%     fprintf(fid1,sprintf('       %0.15f, %0.15f, %0.15f ]\n',kcL(3), kcL(4), kcL(5)));
%     
%     fprintf(fid1,'M2: !!opencv-matrix\n   rows: 3\n   cols: 3\n   dt: d\n');
%     fprintf(fid1,sprintf('   data: [ %0.15f, %0.15f,  %0.15f, %0.15f, \n',intrMatR(1,1), intrMatR(1,2), intrMatR(1,3), intrMatR(2,1)));
%     fprintf(fid1,sprintf('       %0.15f, %0.15f, %0.15f, %0.15f, %0.15f ]\n',intrMatR(2,2), intrMatR(2,3),intrMatR(3,1), intrMatR(3,2), intrMatR(3,3)));
%     fprintf(fid1,'D2: !!opencv-matrix\n   rows: 1\n   cols: 5\n   dt: d\n');
%     fprintf(fid1,sprintf('   data: [ %0.15f, %0.15f, \n',kcR(1), kcR(2)));
%     fprintf(fid1,sprintf('       %0.15f, %0.15f, %0.15f ]\n',kcR(3), kcR(4), kcR(5)));
%     
%     fclose(fid1);
    
    
    
    
end




end