function WriteVid(inputDir, gifName, varargin)

% WriteVideo('G:\matlab\data\direct\gt\D2_011\imgs_bad', 'test.avi')
% WriteVid('G:\matlab\data\direct\gt\D2_011\imgs_bad_full', 'test_comp.avi','G:\matlab\data\direct\gt\D2_011\imgs_good_full');

inputDir2 = [];
if nargin == 3
    inputDir2 = varargin{1};
else
    
end



dirInfo = dir(fullfile(inputDir, '*.png'));

nImages = length(dirInfo);
if(~isempty(inputDir2))
    dirInfo2 = dir(fullfile(inputDir2, '*.png'));    
    nImages = min([nImages length(dirInfo2)]);
end




% filename = fullfile(inputDir, gifName); % Specify the output file name
% [a,b,c] = size(imread(fullfile(inputDir, dirInfo(end).name)));



v = VideoWriter(fullfile(inputDir, gifName),'Uncompressed AVI');
% v = VideoWriter(fullfile(inputDir, gifName),'Motion JPEG AVI');
% v = VideoWriter(fullfile(inputDir, gifName),'Motion JPEG 2000');
% v.Quality = 95;
v.FrameRate = 10;
% v = VideoWriter(fullfile(inputDir, gifName));
open(v);
repeat = 1;
for idx = 1 : nImages
%     img = imresize(imread(fullfile(inputDir, dirInfo(idx).name)),[360,640]);
%     img = imresize(imread(fullfile(inputDir, dirInfo(idx).name)),[922,1914]);
%     img = imresize(imread(fullfile(inputDir, dirInfo(idx).name)),[1040,1920]);
%     
%     img = imresize(imread(fullfile(inputDir, dirInfo(idx).name)),[480,640]);
%     img = imresize(imread(fullfile(inputDir, dirInfo(idx).name)),[480,1280]);
%     img = imresize(imread(fullfile(inputDir, dirInfo(idx).name)),round([640,1706])./1);
    img = imread(fullfile(inputDir, dirInfo(idx).name));
    img  = img(29:1431, 567:2409,:);
    if(~isempty(inputDir2))
        img2 = imread(fullfile(inputDir2, dirInfo(idx).name));
        img2  = img2(29:1431, 567:2409,:);
        img = [img img2];
    end
%     for j = 1:repeat
        writeVideo(v,img);
%     end
%     img = img(:,1:1706/2,:);
%     if(ndims(img) ~= 3)
%         img = cat(3, img, img, img);
%     end
% % %     img = imresize(img(21:386,143:583,:), 2);
%     [A,map] = rgb2ind(img,256);
%     if idx == 1
%         for j = 1
%             imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',delay);
%         end
%     else
% %         imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',2);
%         imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',delay);
%     end
end
close(v);
end