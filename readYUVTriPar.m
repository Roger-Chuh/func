function readYUVTriPar( inputDir , inputDirRGB, width, height,scale)



global cfg
if length(scale) == 1
    scale = [scale,scale];
end
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

width_rgb = cfg.img_width_rgb;
height_rgb = cfg.img_height_rgb;

try
    cfg.check_depth;
catch
    cfg.check_depth = 0;
end

nonYUV = 0;
% close all
if ~cfg.check_depth
    dirInfo = dir(fullfile(inputDir,'*.yuv'));
    dirInfoRGB = dir(fullfile(inputDirRGB,'*.yuv'));
    if length(dirInfo) == 0
        nonYUV = 1;
        dirInfo = dir(fullfile(inputDir,'*.bmp'));
        dirInfoRGB = dir(fullfile(inputDirRGB,'*.bmp'));
    end
    if length(dirInfo) == 0
        dirInfo = dir(fullfile(inputDir,'*.jpg'));
        dirInfoRGB = dir(fullfile(inputDirRGB,'*.jpg'));
    end
    if length(dirInfo) == 0
        dirInfo = dir(fullfile(inputDir,'*.png'));
        dirInfoRGB = dir(fullfile(inputDirRGB,'*.png'));
        return;
    end
    if 0
        dirInfoRGB = dir(fullfile(inputDirRGB,'*.yuv'));
    end
    
else
    read_board_depth(inputDir, width, height);
    return;
end


if 0
    try
        for i = 1 : length(dirInfo)/2
            id_ = find(dirInfo(i).name == '_');
            timeNum1(i,1) = str2double(dirInfo(i).name(id_(1)+1:id_(2)-1));
            
        end
        
        for i = length(dirInfo)/2 + 1 : length(dirInfo)
            id_ = find(dirInfo(i).name == '_');
            timeNum2(i - length(dirInfo)/2,1) = str2double(dirInfo(i).name(id_(1)+1:id_(2)-1));
            
        end
    catch
        skgfj = 1;
    end
end

dircell=struct2cell(dirInfo);
aaAa = dircell(2,:)';
%     sss = datetime(aaAa,'Format','HH:mm:ss','Locale','zh_CN');

if 1
    %% 20201112 big bug,(change ## if 0 -> if 1 ##) why wait till just now!!!
    [~, ind1] = sort([dirInfo(1:length(dirInfo)/2).datenum], 'ascend');
    [~, ind2] = sort([dirInfo(length(dirInfo)/2+1:end).datenum], 'ascend');
    ind2 = ind2 + length(dirInfo)/2;
    I = [ind1 ind2];
    
    tempInfo = cell(length(I),1);
    for fg = 1:length(I)
        tempInfo{fg,1} = dirInfo(I(fg)).name;
    end
    for fgg = 1:length(I)
        dirInfo(fgg).name = tempInfo{fgg};
    end
    
    %%
    if 1
        
        
        %% 20201112 big bug,(change ## if 0 -> if 1 ##) why wait till just now!!!
        [~, ind1RGB] = sort([dirInfoRGB(1:length(dirInfoRGB)/1).datenum], 'ascend');
        %             [~, ind2] = sort([dirInfo(length(dirInfo)/2+1:end).datenum], 'ascend');
        %             ind2 = ind2 + length(dirInfo)/2;
        %             I = [ind1 ind2];
        IRGB = [ind1RGB];
        
        tempInfoRGB = cell(length(IRGB),1);
        for fg = 1:length(IRGB)
            tempInfoRGB{fg,1} = dirInfoRGB(IRGB(fg)).name;
        end
        for fgg = 1:length(IRGB)
            dirInfoRGB(fgg).name = tempInfoRGB{fgg};
        end
    end
    
    
    
end

outFlag = [];
% % outFlag = [1 2 17 19 20 21 22 37 39 40 41 42  43 60];
% % outFlag = [3 4 6 7 8 11 13 14 15 17 19 20 23 24 27 28 35 36 37 38 40 41 45 48 49 52 53 56 57 58 61 67 68 70 71 73 75 76];
outFlag = []; %setdiff([1:40]',[1;2;3;5;6;8;13;14;15;16;17;18;19;20;21;22;23;24;25;26;27;28;29;30;31;32;33;34;35;36;37;38;39;40]);
frameNum = length(dirInfo)/2;



thr = 130;
H = fspecial('gaussian',21);
cnt = 1;


if cfg.patternNum == 2
    colRng1 = round([width- width/(cfg.patternNum+1) : width]);
    colRng2 = round([1 : width/(cfg.patternNum+1)]);
else
    colRng1 = [];
    colRng2 = [];
end

patternNum = cfg.patternNum;
blurLevel = cfg.blurLevel;


idList = [1:frameNum 1:frameNum 1:frameNum]';

dirInfoCat = [dirInfo; dirInfoRGB];
inputDirCat = [repmat({inputDir}, 2*frameNum,1); repmat({inputDirRGB}, frameNum,1)];

whichCat = [repmat({'L'}, frameNum,1); repmat({'R'}, frameNum,1); repmat({'R'}, frameNum,1)];
widthList = [repmat(width, 2*frameNum, 1);repmat(width_rgb, frameNum, 1)];
heightList = [repmat(height, 2*frameNum, 1);repmat(height_rgb, frameNum, 1)];
% for i 
if ~nonYUV
    if cfg.coreNum > 0
        parfor i =    1 : 3*frameNum  %[1 2 4 6 8 13 14 16 18 20 25 26 28 30 35 36 37 39 41 47 48 49 51 53 55 58 60 61 63 65 67]
            
            readFile(inputDirCat{i}, dirInfoCat(i).name,  idList(i), widthList(i), heightList(i),colRng1,colRng2, patternNum, blurLevel, scale, whichCat{i});
            
            
        end
    else
        for i =    1 : 3*frameNum  %[1 2 4 6 8 13 14 16 18 20 25 26 28 30 35 36 37 39 41 47 48 49 51 53 55 58 60 61 63 65 67]
            
            readFile(inputDirCat{i}, dirInfoCat(i).name,  idList(i),widthList(i), heightList(i),colRng1,colRng2, patternNum, blurLevel, scale, whichCat{i});
            
            
        end
    end
    
else
    if cfg.coreNum > 0
        parfor i =    1 : 3*frameNum  %[1 2 4 6 8 13 14 16 18 20 25 26 28 30 35 36 37 39 41 47 48 49 51 53 55 58 60 61 63 65 67]
            
            renameFile(inputDirCat{i}, dirInfoCat(i).name,  idList(i), widthList(i), heightList(i),colRng1,colRng2, patternNum, blurLevel, scale, whichCat{i});
            
            
        end
    else
        for i =    1 : 3*frameNum  %[1 2 4 6 8 13 14 16 18 20 25 26 28 30 35 36 37 39 41 47 48 49 51 53 55 58 60 61 63 65 67]
            
            renameFile(inputDirCat{i}, dirInfoCat(i).name,  idList(i), widthList(i), heightList(i),colRng1,colRng2, patternNum, blurLevel, scale, whichCat{i});
            
            
        end
    end
end

% % for i = 1 : length(dirInfo)/2
% %     id(i,1) = str2double(dirInfo(i).name(6:end-4));
% % end
end
function renameFile(inputDir, filename, i, width, height,colRng1,colRng2, patternNum, blurLevel, scale, which)
 CC = imread(fullfile(inputDir,filename));
 if 1
     delete(fullfile(inputDir,filename));
 end
%  DD = dispI420(fullfile(inputDir,dirInfo(i + frameNum).name),width,height);
%  EE = dispI420(fullfile(inputDir,dirInfo(i + frameNum).name),width,height);
    
    if 0
        CC = imfilter(uint8(CC), H);
        DD = imfilter(uint8(DD), H);
    elseif 0
        CC = imgaussfilt(uint8(CC), 1.4);
%         DD = imgaussfilt(uint8(DD), 1.4);
    elseif 1
        if blurLevel ~= 0
            CC = imgaussfilt(uint8(CC), blurLevel);
%             DD = imgaussfilt(uint8(DD), blurLevel);
        else
            
        end
    else
        sadghk = 1;
    end
   
    if  patternNum == 1 % mod(i,2) == 0
        imwrite(imresize(uint8(CC),[height/scale(1),width/scale(2)]), fullfile(inputDir,sprintf('img%s_%05d.png',which,i)));
%         imwrite(imresize(uint8(DD),[height/scale(1),width/scale(2)]), fullfile(inputDir,sprintf('imgR_%05d.png',i)));
    else
        CC_ = CC; %DD_ = DD;
        CC_(:,colRng1,:) = 0;
%         DD_(:,colRng1,:) = 0;
        imwrite(imresize(uint8(CC_),[height/scale(1),width/scale(2)]), fullfile(inputDir,sprintf('img%s_%05d.png',which,2*i-1)));
%         imwrite(imresize(uint8(DD_),[height/scale(1),width/scale(2)]), fullfile(inputDir,sprintf('imgR_%05d.png',2*i-1)));
%         cnt = cnt + 1;
        CC_ = CC; %DD_ = DD;
        CC_(:,colRng2,:) = 0;
%         DD_(:,colRng2,:) = 0;
        imwrite(imresize(uint8(CC_),[height/scale(1),width/scale(2)]), fullfile(inputDir,sprintf('img%s_%05d.png',which,2*i)));
%         imwrite(imresize(uint8(DD_),[height/scale(1),width/scale(2)]), fullfile(inputDir,sprintf('imgR_%05d.png',2*i)));
%         cnt = cnt + 1;
    end
end

function readFile(inputDir, filename, i, width, height,colRng1,colRng2, patternNum, blurLevel, scale, which)
 CC = dispI420(fullfile(inputDir,filename),width,height);
%  DD = dispI420(fullfile(inputDir,dirInfo(i + frameNum).name),width,height);
%  EE = dispI420(fullfile(inputDir,dirInfo(i + frameNum).name),width,height);
    
    if 0
        CC = imfilter(uint8(CC), H);
        DD = imfilter(uint8(DD), H);
    elseif 0
        CC = imgaussfilt(uint8(CC), 1.4);
%         DD = imgaussfilt(uint8(DD), 1.4);
    elseif 1
        if blurLevel ~= 0
            CC = imgaussfilt(uint8(CC), blurLevel);
%             DD = imgaussfilt(uint8(DD), blurLevel);
        else
            
        end
    else
        sadghk = 1;
    end
   
    if  patternNum == 1 % mod(i,2) == 0
        imwrite(imresize(uint8(CC),[height/scale(1),width/scale(2)]), fullfile(inputDir,sprintf('img%s_%05d.png',which,i)));
%         imwrite(imresize(uint8(DD),[height/scale(1),width/scale(2)]), fullfile(inputDir,sprintf('imgR_%05d.png',i)));
    else
        CC_ = CC; %DD_ = DD;
        CC_(:,colRng1,:) = 0;
%         DD_(:,colRng1,:) = 0;
        imwrite(imresize(uint8(CC_),[height/scale(1),width/scale(2)]), fullfile(inputDir,sprintf('img%s_%05d.png',which,2*i-1)));
%         imwrite(imresize(uint8(DD_),[height/scale(1),width/scale(2)]), fullfile(inputDir,sprintf('imgR_%05d.png',2*i-1)));
%         cnt = cnt + 1;
        CC_ = CC; %DD_ = DD;
        CC_(:,colRng2,:) = 0;
%         DD_(:,colRng2,:) = 0;
        imwrite(imresize(uint8(CC_),[height/scale(1),width/scale(2)]), fullfile(inputDir,sprintf('img%s_%05d.png',which,2*i)));
%         imwrite(imresize(uint8(DD_),[height/scale(1),width/scale(2)]), fullfile(inputDir,sprintf('imgR_%05d.png',2*i)));
%         cnt = cnt + 1;
    end
end


function rgbImg = dispI420(inputfile,width,height)
fd = fopen(inputfile, 'rb');
if fd == -1
    error('open file error');
end
% width = 1920;
% height = 1080;
y = fread(fd, width*height,'uint8');
u = fread(fd, width*height/4,'uint8');
v = fread(fd, width*height/4,'uint8');
fclose(fd);
yimg = reshape(y, [width,height]);
uimg = reshape(u, [width/2,height/2]);
vimg = reshape(v, [width/2,height/2]);
if 0
    rgbImg = Yuv420P2Rgb(yimg', uimg', vimg');
else
    rgbImg = cat(3, yimg', yimg', yimg');
end
% rgbImg = im2frame(uint8(rgbImg));
% rgbImg = (uint8(rgbImg));
% % figure(1),movie(fg);
end
function rgbImg = Yuv420P2Rgb(y, u, v)
% R = Y + 1.402 * (V - 128);
% G = Y - 0.34414 * (U - 128) - 0.71414 * (V - 128);
% B = Y + 1.772 * (U - 128);
% I420: YYYYYYYY UU VV =>YUV420P
% YV12: YYYYYYYY VV UU =>YUV420P
% 4:1:1

rImg = zeros(size(y));
gImg = zeros(size(y));
bImg = zeros(size(y));
rImg(1:2:end,1:2:end) = y(1:2:end,1:2:end) + 1.0402 * (v - 128);
rImg(2:2:end,1:2:end) = y(2:2:end,1:2:end) + 1.0402 * (v - 128);
rImg(2:2:end,2:2:end) = y(2:2:end,2:2:end) + 1.0402 * (v - 128);
rImg(1:2:end,2:2:end) = y(1:2:end,2:2:end) + 1.0402 * (v - 128);

gImg(1:2:end,1:2:end) = y(1:2:end,1:2:end) + 0.34414 * (u - 128) - 0.71414 * (v - 128);
gImg(2:2:end,1:2:end) = y(2:2:end,1:2:end) + 0.34414 * (u - 128) - 0.71414 * (v - 128);
gImg(2:2:end,2:2:end) = y(2:2:end,2:2:end) + 0.34414 * (u - 128) - 0.71414 * (v - 128);
gImg(1:2:end,2:2:end) = y(1:2:end,2:2:end) + 0.34414 * (u - 128) - 0.71414 * (v - 128);

bImg(1:2:end,1:2:end) = y(1:2:end,1:2:end) + 1.772 * (u - 128);
bImg(2:2:end,1:2:end) = y(2:2:end,1:2:end) + 1.772 * (u - 128);
bImg(2:2:end,2:2:end) = y(2:2:end,2:2:end) + 1.772 * (u - 128);
bImg(1:2:end,2:2:end) = y(1:2:end,2:2:end) + 1.772 * (u - 128);

rgbImg = zeros(size(y,1), size(y,2), 3);
rgbImg(:,:,1) = rImg;
rgbImg(:,:,2) = gImg;
rgbImg(:,:,3) = bImg;
end
function CC = Precess(CC, thr)

 CC = rgb2gray(CC);
 CC(CC > thr) = 0;
    CC = im2double(CC);
    CC = CC./max(CC(:));
    CC= uint8(255.*(CC));
    CC = cat(3,CC,CC,CC);

end