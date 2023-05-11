function CropImagePair()
inputDir ='G:\matlab\data\omni';
MakeDirIfMissing(fullfile(inputDir,'01'));
dirInfo = dir(fullfile(inputDir,'*.png'));
% for i = 1 : length(dirInfo)
%    img = rgb2gray(imread(fullfile(inputDir,dirInfo(i).name))); 
%    imgL = img(:,1:size(img,2)/2);
%    imgR = img(:,size(img,2)/2+1:end);
% %    imwrite(imgL,fullfile(inputDir,'1',sprintf('imgL_%03d.png',i)));
% %    imwrite(imgR,fullfile(inputDir,'1',sprintf('imgR_%03d.png',i)));
% end

names = {'Camera0','Camera1'};
srcs = {'src0','src1'};
for i = 1 : length(dirInfo)
   img = rgb2gray(imread(fullfile(inputDir,dirInfo(i).name))); 
   imgL = img(:,1:size(img,2)/2);
   imgR = img(:,size(img,2)/2+1:end);
%    imwrite(imgL,fullfile(inputDir,'1',sprintf('imgL_%03d.png',i)));
%    imwrite(imgR,fullfile(inputDir,'1',sprintf('imgR_%03d.png',i)));
end


camInfo = dir(fullfile(inputDir,'0*'));

tsInfo = cell(4,2);
frameNum = [];
for i = 1 : 2%length(camInfo)
    tsInfo{i,1} = loadjson(fullfile(inputDir, camInfo(i).name,'data.json'));
    tsInfo{i,1}.Sequence.CameraInfo.width = 960;512;
    tsInfo{i,1}.Sequence.CameraInfo.height = 960;512;
    tsInfo{i,1}.Sequence.CameraInfo.Camera.width = 960;512;
    tsInfo{i,1}.Sequence.CameraInfo.Camera.height = 960;512;
    srcInfo = dir(fullfile(inputDir,srcs{i},'*.png'));
    for j = 1:length(tsInfo{i,1}.Sequence.Frameset.Frame)
        tsInfo{i,2} = [tsInfo{i,2};tsInfo{i,1}.Sequence.Frameset.Frame(j).timestamp];
        if j <= length(dirInfo)
             img = (imread(fullfile(inputDir,srcs{i},srcInfo(j).name))); 
        else
            img = (imread(fullfile(inputDir,srcs{i},srcInfo(length(dirInfo)).name))); 
        end
        imwrite(img,fullfile(inputDir, camInfo(i).name,tsInfo{i,1}.Sequence.Frameset.Frame(j).filename));
        frameNum = [frameNum;length(tsInfo{i,1}.Sequence.Frameset.Frame)];
    end
    
end

err = [tsInfo{1,2}(1:min(frameNum),1) tsInfo{2,2}(1:min(frameNum),1) tsInfo{3,2}(1:min(frameNum),1)] - tsInfo{4,2}(1:min(frameNum),1);
figure, subplot(1,2,1);plot(diff(tsInfo{1,2})./1e6);subplot(1,2,2), plot(err./1e6)


end