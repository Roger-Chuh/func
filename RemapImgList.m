function RemapImgList(inputDir, paraDir)

% example: RemapImgList('./test', './')



dirInfo = dir(fullfile(inputDir, 'img*.png'));
if length(dirInfo) == 0
    dirInfo = dir(fullfile(inputDir, '*.jpg'));
end
load(fullfile(paraDir,'remap.mat'));

frameNum = length(dirInfo)/2;

for i = 1 : frameNum
    
    imgR = imread(fullfile(inputDir, dirInfo(i).name));
    imgL = imread(fullfile(inputDir, dirInfo(i + frameNum).name));
    
    if ndims(imgL) == 3
        rectImgL(:,:,1) = uint8(interp2(xGrid, yGrid, double(imgL(:,:,1)),xRect2OrigL,yRect2OrigL));
        rectImgR(:,:,1) = uint8(interp2(xGrid, yGrid, double(imgR(:,:,1)),xRect2OrigR,yRect2OrigR));
        
        rectImgL(:,:,2) = uint8(interp2(xGrid, yGrid, double(imgL(:,:,2)),xRect2OrigL,yRect2OrigL));
        rectImgR(:,:,2) = uint8(interp2(xGrid, yGrid, double(imgR(:,:,2)),xRect2OrigR,yRect2OrigR));
        
        rectImgL(:,:,3) = uint8(interp2(xGrid, yGrid, double(imgL(:,:,3)),xRect2OrigL,yRect2OrigL));
        rectImgR(:,:,3) = uint8(interp2(xGrid, yGrid, double(imgR(:,:,3)),xRect2OrigR,yRect2OrigR));
        
    else
        rectImgL(:,:,1) = uint8(interp2(xGrid, yGrid, double(imgL(:,:,1)),xRect2OrigL,yRect2OrigL));
        rectImgR(:,:,1) = uint8(interp2(xGrid, yGrid, double(imgR(:,:,1)),xRect2OrigR,yRect2OrigR));
    end
    
    imwrite(rectImgL, fullfile(inputDir, sprintf('rectL_%05d.png',i)));
    imwrite(rectImgR, fullfile(inputDir, sprintf('rectR_%05d.png',i)));
    
end
figure,imshow([rectImgL rectImgR]);
end