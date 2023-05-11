function [rectImgL, rectImgR, xRangeL, xRangeR,maskL,maskR] = RectifyImagePair(stereoParam, imgLeft, imgRight)
% RectifyImagePair performs rectification on image pair specified by
% imgFileLeft and imgFileRight.
% By Ji Zhou

if ischar(imgLeft)
    imgL = imread(imgLeft);
else
    imgL = imgLeft;
end
if ischar(imgRight)
    imgR = imread(imgRight);
else
    imgR = imgRight;
end
imgSize = size(imgLeft);
imgSize = imgSize(1:2);

if 0 %% 20200702
    [rectParamL, rectParamR] = GetRectifyParam(stereoParam, size(imgL));
else
    [rectParamL, rectParamR] = GetRectifyParam2(stereoParam, size(imgL));
end

[rectImgL, xRangeL0, maskL] = RectifyOneImage(imgL, rectParamL);
[rectImgR, xRangeR0, maskR] = RectifyOneImage(imgR, rectParamR);

dltL = imgSize(2) - xRangeL0(2);
dltR = xRangeR0(1);
dlt = min(dltL,dltR);
if 1
    if abs(diff(xRangeL0)) >  abs(diff(xRangeR0))
        xRangeR0(1) = xRangeR0(1) - abs((abs(diff(xRangeL0)) -  abs(diff(xRangeR0))));
    else
        xRangeL0(2) = xRangeL0(2) + abs((abs(diff(xRangeL0)) -  abs(diff(xRangeR0))));
    end
    xRangeR0 = xRangeR0 + size(imgL,2);
    xRangeL = xRangeL0;
    xRangeR = xRangeR0;
else
    xRangeL = [1 imgSize(2) + dlt];
    xRangeR = [imgSize(2)-dlt+1 imgSize(2)*2];
end


% % % % % % % % % % % % %  ShowStereoPair(rectImgL, rectImgR);
% % % % imshow([rectImgL rectImgR]);
% % % % figure, imshowpair(rectImgL,rectImgR);
end


function ShowStereoPair(imgLeft, imgRight)

[imgH,imgW, numChnl] = size(imgLeft);
splitW = 10;
imgDisp = zeros(imgH, 2*imgW+splitW, numChnl, 'uint8');
imgDisp(:,1:imgW, :) = imgLeft;
imgDisp(:,(1:imgW)+splitW+imgW, :) = imgRight;
figure, imshow(imgDisp);


end