function plotDatasetTimestamps()

inputDir = 'G:\matlab\data\direct\gt2\4';

camInfo = dir(fullfile(inputDir, 'Camera*'));
dirCam0 = dir(fullfile(inputDir, camInfo(1).name,'images','*.bmp'));
dirCam1 = dir(fullfile(inputDir, camInfo(2).name,'images','*.bmp'));
dirCam2 = dir(fullfile(inputDir, camInfo(3).name,'images','*.bmp'));
dirCam3 = dir(fullfile(inputDir, camInfo(4).name,'images','*.bmp'));
dirCams = {dirCam0, dirCam1, dirCam2, dirCam3};

timestamp1 = zeros(length(dirCam0),1);
for i = 1 : length(dirCam0)
    timestamp1(i,1) = str2double(dirCam0(i).name(1:end-4)) * 1e-6;
end
timestamp2 = zeros(length(dirCam1),1);
for i = 1 : length(dirCam1)
    timestamp2(i,1) = str2double(dirCam1(i).name(1:end-4)) * 1e-6;
end
timestamp3 = zeros(length(dirCam2),1);
for i = 1 : length(dirCam2)
    timestamp3(i,1) = str2double(dirCam2(i).name(1:end-4)) * 1e-6;
end
timestamp4 = zeros(length(dirCam3),1);
for i = 1 : length(dirCam3)
    timestamp4(i,1) = str2double(dirCam3(i).name(1:end-4)) * 1e-6;
end



inputDir = 'G:\matlab\data\direct\gt\D2_011\4\imgs4';

depthDir = dir(fullfile(inputDir, '0*'));
timestamp_depth = zeros(length(depthDir),1);
for i = 1 : length(depthDir)
    timestamp_depth(i,1) = str2double(depthDir(i).name(7:end-19)) * 1e-6;
end

traceDir = dir(fullfile(inputDir, 'Trace*'));
timestamp_trace = zeros(length(traceDir),1);
for i = 1 : length(traceDir)
    id = find(traceDir(i).name == '_');
    timestamp_trace(i,1) = str2double(traceDir(i).name(id(1)+1:id(2)-1)) * 1e-6;
end

end