function readCalibData(inputDir)

dirInfo = dir(fullfile(inputDir, '*.y*'));

calib = ReadYaml(fullfile(inputDir, dirInfo(1).name));

stereoParam.focLeft = [calib.CameraMatrixL.data(1,1); calib.CameraMatrixL.data(2,2)];
stereoParam.cenLeft = [calib.CameraMatrixL.data(1,3); calib.CameraMatrixL.data(2,3)];
stereoParam.alphaLeft = 0;
stereoParam.kcLeft = calib.DistCoeffsL.data';

stereoParam.focRight = [calib.CameraMatrixR.data(1,1); calib.CameraMatrixR.data(2,2)];
stereoParam.cenRight = [calib.CameraMatrixR.data(1,3); calib.CameraMatrixR.data(2,3)];
stereoParam.alphaRight = 0;
stereoParam.kcRight = calib.DistCoeffsR.data';

stereoParam.rotVecRef = rodrigues(calib.R_stereo.data);
stereoParam.transVecRef = calib.T_stereo.data';

stereoParam.rotMatLeft = calib.R1.data;
stereoParam.rotMatRight = calib.R2.data;
stereoParam.intrMatLeftNew = calib.P1.data(1:3,1:3);
stereoParam.intrMatRightNew = calib.P2.data(1:3,1:3);

stereoParam.CameraMatrixRGB = calib.CameraMatrixRGB.data;
stereoParam.CameraMatrixRGB_new = calib.CameraMatrixRGB_new.data;
stereoParam.DistCoeffsRGB = calib.DistCoeffsRGB.data';
stereoParam.R_RGB = rodrigues(calib.R_RGB.data);
stereoParam.T_RGB = calib.T_RGB.data';

stereoParam_ = stereoParam;
save(fullfile(inputDir,'stereoParam_.mat'),'stereoParam_');

end