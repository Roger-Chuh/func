function frame2video_temp()
inputDir = 'D:\Temp\20201221\1\video\video';

inputDir = 'D:\Work\marker_mapper1.0.15\20201228_11633\video\L';
dirInfo = dir(fullfile(inputDir, '*.png'));

if 0
    [~, ind1] = sort([dirInfo(1:length(dirInfo)).datenum], 'ascend');
    
    I = [ind1];
    
    tempInfo = cell(length(I),1);
    for fg = 1:length(I)
        tempInfo{fg,1} = dirInfo(I(fg)).name;
    end
    for fgg = 1:length(I)
        dirInfo(fgg).name = tempInfo{fgg};
    end
end

if 0
    for k = 1 : length(dirInfo)
        id = find(dirInfo(k).name == '_');
        Id(k,1) = str2double(dirInfo(k).name(id+1:end-4)) ;
        
    end
    
    [~, sortId] = sort(Id, 'ascend');
    
    I = [sortId];
    
    tempInfo = cell(length(I),1);
    for fg = 1:length(I)
        tempInfo{fg,1} = dirInfo(I(fg)).name;
    end
    for fgg = 1:length(I)
        dirInfo(fgg).name = tempInfo{fgg};
    end
end

% dirInfo = dir(fullfile(inputDir, '*.png'));
img = imresize(imread(fullfile(inputDir, dirInfo(1).name)), [540 960]);
img = double(imread(fullfile(inputDir, dirInfo(1).name)));
img(img > 180) = 1;
img = img./max(img(:));
img = uint8(255*img);
imgg = cat(3, img,img,img);

aviobj=VideoWriter(fullfile(inputDir, 'demo.avi'));
aviobj.FrameRate = 15;  25;
open(aviobj);


for i = 1 : length(dirInfo)
    imgg = imread(fullfile(inputDir, dirInfo(i).name));
    imgg(1:550,:) = 0;
     writeVideo(aviobj, imgg);
end
close(aviobj);

end