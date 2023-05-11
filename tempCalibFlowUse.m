function [imgL, imgR, rectImgL, rectImgR, meanError, KK_newL, imgLRGB, imgRRGB] = tempCalibFlowUse(inputDir,scalee, yuvv,widthh,heightt,calibFuncDir,cbSizee,switchLRR)
global cfg inputDir0 paraDir0

width = str2double(widthh);
height = str2double(heightt);
scale = str2double(scalee);
cbSize = str2double(cbSizee);
  
readYUVDualPar( inputDir , width, height,scale)


%����һ������������Ϣ�Ľṹ��config
config.dX = cbSize; 
config.dY = cbSize; 
config.estDistortion = [1;1;1;1;cfg.est_k3];

% ���aruco�ǵ㲢������������
[imgL, imgR, rectImgL, rectImgR, meanError, KK_newL, imgLRGB, imgRRGB] = CallArucoDetection(inputDir, fullfile(calibFuncDir,'Aruco'), paraDir0, config);

end

