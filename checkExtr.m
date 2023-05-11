function checkExtr()
inputDir = 'G:\matlab\data\extr\1';
dirInfo = dir(fullfile(inputDir,'0*'));

for i = 1 : length(dirInfo)
    
    info = dir(fullfile(inputDir, dirInfo(i).name, '*.bmp'));
    if length(info) == 0
        info = dir(fullfile(inputDir, dirInfo(i).name, '*.png'));
    end
    
    frameNum = floor(length(info)/2);
    
    for j = 1 :1% length(info)
        img1 = imread(fullfile(inputDir, dirInfo(i).name, info(j).name));
        for jj = j+1:length(info)
            img2 = imread(fullfile(inputDir, dirInfo(i).name, info(jj).name));
            figure(1),imshowpair(img1, img2);title(sprintf('camera %d, image [%d %d]',i, j, jj));drawnow;
        end
    end
    
end

end