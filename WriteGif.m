function WriteGif(inputDir, gifName,delay)
dirInfo = dir(fullfile(inputDir, '*.png'));

if length(dirInfo) == 0
    dirInfo = dir(fullfile(inputDir, '*.bmp'));
end

nImages = length(dirInfo);
filename = fullfile(inputDir, gifName); % Specify the output file name
[a,b,c] = size(imread(fullfile(inputDir, dirInfo(end).name)));

for idx = 1 : nImages
%     img = imresize(imread(fullfile(inputDir, dirInfo(idx).name)),[360,640]);
%     img = imresize(imread(fullfile(inputDir, dirInfo(idx).name)),[922,1914]);
    img = imresize(imread(fullfile(inputDir, dirInfo(idx).name)),[1040,1920]);
    
    img = imresize(imread(fullfile(inputDir, dirInfo(idx).name)),[480,640]);
    img = imresize(imread(fullfile(inputDir, dirInfo(idx).name)),[480,1280]);
    img = imresize(imread(fullfile(inputDir, dirInfo(idx).name)),round([640,1706])./1);
    img = imresize(imread(fullfile(inputDir, dirInfo(idx).name)),round([960,1280])./1);
    img = imresize(imread(fullfile(inputDir, dirInfo(idx).name)),round([740,920])./1);
%     img = img(:,1:1706/2,:);
    if(ndims(img) ~= 3)
        img = cat(3, img, img, img);
    end
% %     img = imresize(img(21:386,143:583,:), 2);
    [A,map] = rgb2ind(img,256);
    if idx == 1
        for j = 1
            imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',delay);
        end
    else
%         imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',2);
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',delay);
    end
end
end