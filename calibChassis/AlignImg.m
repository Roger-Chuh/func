function AlignImg()

close all

inputDir = 'C:\Users\rongjiezhu\Documents\WXWork\1688852647983342\Cache\File\2021-07\usbcam�궨ͼ\ir\test03';
load(fullfile(inputDir, 'calib.mat'));

dirInfo = dir(fullfile(inputDir, 'img*.jpg'));
frameNum = length(dirInfo)/2 - 1;
imgSize = [640 480];
[xGrid, yGrid] = meshgrid(1 : imgSize(2), 1 : imgSize(1));
pix = [xGrid(:) yGrid(:)];

intrMatL = [stereoParam.focLeft(1) 0 stereoParam.cenLeft(1); 0 stereoParam.focLeft(2) stereoParam.cenLeft(2); 0 0 1];
intrMatR = [stereoParam.focRight(1) 0 stereoParam.cenRight(1); 0 stereoParam.focRight(2) stereoParam.cenRight(2); 0 0 1];

kcL = stereoParam.kcLeft;
kcR = stereoParam.kcRight;

rotMat = rodrigues(stereoParam.rotVecRef);

intrMatNew = intrMatL;
% intrMatNew(1,1) = 400;
% intrMatNew(2,2) = 400;

pixUndistL = remapRect(pix', intrMatNew, intrMatL, kcL, eye(3));
pixUndistR = remapRect(pix', intrMatNew, intrMatR, kcR, eye(3));
reMapLX = reshape(pixUndistL(:,1), imgSize);
reMapLY = reshape(pixUndistL(:,2), imgSize);
reMapRX = reshape(pixUndistR(:,1), imgSize);
reMapRY = reshape(pixUndistR(:,2), imgSize);

for i = frameNum : frameNum
    imgL = rgb2gray(imread(fullfile(inputDir, dirInfo(i).name)));
    imgR = rgb2gray(imread(fullfile(inputDir, dirInfo(i+frameNum).name)));
    imgL_undist = uint8(interp2(xGrid, yGrid, double(imgL(:,:,1)),reMapLX,reMapLY));
    imgR_undist = uint8(interp2(xGrid, yGrid, double(imgR(:,:,1)),reMapRX,reMapRY));
    rgb2ir = ImgTransform(imgR_undist, intrMatNew, rotMat');
    
    figure,subplot(1,2,1);imshowpair(rgb2ir, imgL_undist);subplot(1,2,2),imshowpair(rgb2ir, imgR_undist)
end

end