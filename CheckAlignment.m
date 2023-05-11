function CheckAlignment()
% function CheckAlignment(inputDir, paraDir, varargin)

inputDir = 'G:\matlab\data\alignment\check_align';
paraDir = 'G:\matlab\data\alignment\data';

global cfg

if 0
    if (nargin == 2)
        draww = '0';
    elseif (nargin == 3)
        draww = varargin{1};
    else
        error('Too many input arguments');
    end
end

draww = '1';


draw = str2double(draww);

MakeDirIfMissing(fullfile(inputDir,'align'));
delete(fullfile(inputDir,'align','*'));

close all
load(fullfile(paraDir,'rectMap.mat'));




readYUVTriPar( inputDir ,fullfile(inputDir,'RGB'), cfg.img_width, cfg.img_height, cfg.img_scale);


dirInfo = dir(fullfile(inputDir,'img*.png'));
frameNum = length(dirInfo)/2;

dirInfoRGB = dir(fullfile(inputDir,'RGB','img*.png'));

dispRange = 256;
areaRatio = 0.1; 0.2;
downScaleRatio = 1; 0.5;

intrMatL = cfg.LRIntrMat;
baseline = cfg.LRBaseline;
intrMatRGB = cfg.RGBIntrMat;
L2RGB = cfg.L2RGB;

alignDiffList = [];
% tic;
for i = 1 : frameNum
    
    if cfg.switch_lr
        imgR = imread(fullfile(inputDir, dirInfo(i).name));
        imgL = imread(fullfile(inputDir, dirInfo(i + frameNum).name));
    else
        imgL = imread(fullfile(inputDir, dirInfo(i).name));
        imgR = imread(fullfile(inputDir, dirInfo(i + frameNum).name));
    end
    imgRGB = imread(fullfile(inputDir,'RGB', dirInfoRGB(i).name));
    if ndims(imgL) == 3
        rectImgL(:,:,1) = uint8(interp2(matX, matY, double(imgL(:,:,1)),reMapLX,reMapLY));
        rectImgR(:,:,1) = uint8(interp2(matX, matY, double(imgR(:,:,1)),reMapRX,reMapRY));
        rectImgRGB(:,:,1) = uint8(interp2(matRGBX, matRGBY, double(imgRGB(:,:,1)),reMapRGBX,reMapRGBY));
        
        rectImgL(:,:,2) = uint8(interp2(matX, matY, double(imgL(:,:,2)),reMapLX,reMapLY));
        rectImgR(:,:,2) = uint8(interp2(matX, matY, double(imgR(:,:,2)),reMapRX,reMapRY));
        rectImgRGB(:,:,2) = uint8(interp2(matRGBX, matRGBY, double(imgRGB(:,:,2)),reMapRGBX,reMapRGBY));
        
        rectImgL(:,:,3) = uint8(interp2(matX, matY, double(imgL(:,:,3)),reMapLX,reMapLY));
        rectImgR(:,:,3) = uint8(interp2(matX, matY, double(imgR(:,:,3)),reMapRX,reMapRY));
        rectImgRGB(:,:,3) = uint8(interp2(matRGBX, matRGBY, double(imgRGB(:,:,3)),reMapRGBX,reMapRGBY));
        
        dispMap = disparity(rgb2gray(imresize(rectImgL,downScaleRatio)), rgb2gray(imresize(rectImgR,downScaleRatio)),'DisparityRange', [0 dispRange*downScaleRatio]);
        dispMap(dispMap<5) = nan;
        dispMap(dispMap>dispRange*downScaleRatio) = nan;
        dispMap = imresize(dispMap,1/downScaleRatio)./downScaleRatio;
        
        rectImgL = rgb2gray(rectImgL);
        rectImgR = rgb2gray(rectImgR);
        rectImgRGB = rgb2gray(rectImgRGB);
    else
        rectImgL(:,:,1) = uint8(interp2(matX, matY, double(imgL(:,:,1)),reMapLX,reMapLY));
        rectImgR(:,:,1) = uint8(interp2(matX, matY, double(imgR(:,:,1)),reMapRX,reMapRY));
        rectImgRGB(:,:,1) = uint8(interp2(matRGBX, matRGBY, double(imgRGB(:,:,1)),reMapRGBX,reMapRGBY));
        
        dispMap = disparity(imresize(rectImgL,downScaleRatio), imresize(rectImgR,downScaleRatio),'DisparityRange', [0 dispRange*downScaleRatio]);
        
        dispMap(dispMap<5) = nan;
        dispMap(dispMap>dispRange*downScaleRatio) = nan;
        dispMap = imresize(dispMap,1/downScaleRatio)./downScaleRatio;
    end
    dispMap = RemoveDepthHoles(dispMap, areaRatio);
%     imwrite(rectImgL, fullfile(inputDir, sprintf('rectL_%05d.png',i)));
%     imwrite(rectImgR, fullfile(inputDir, sprintf('rectR_%05d.png',i)));
    [imgDiffMean, maskPrv, maskWarp, valid_mask_in_warp_prv, inValidMask_in_cur] = ImageAlignmentError(rectImgL, rectImgRGB, dispMap, dispMap, intrMatL,intrMatRGB,baseline, L2RGB, draw);
    if draw
       saveas(gcf, fullfile(inputDir,'align',sprintf('aligned_img_%04d.png', i))); 
    end
    alignDiffList = [alignDiffList; imgDiffMean];
end
% toc;
% figure,imshow([rectImgL rectImgR]);



fid1 = fopen(fullfile(inputDir,'align_error.txt'),'w');%????

for v = 1 : length(alignDiffList)
    fprintf(fid1,sprintf('%0.5f\n ',alignDiffList(v)));
end

fclose(fid1);

% exit(0);

end