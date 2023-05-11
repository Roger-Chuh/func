function testOmniRectification()
if 0
    inputDir = 'D:\Temp\20190125\sensor210\Left';
    load('D:\Temp\20190125\sensor210\Left\oCamModelL.mat');
    oCamModel = oCamModelL;
    dirInfo = dir(fullfile(inputDir, '*.png'));
else
    inputDir = 'D:\Temp\20210118\fisheye';
    load('D:\Temp\20210118\fisheye\oCamModel.mat')
    dirInfo = dir(fullfile(inputDir, '*.jpg'));
end



img = imread(fullfile(inputDir, dirInfo(1).name));
U_same = ocam_undistort_map_use(oCamModel,'OutputView','same');
intrMat_same = U_same.K';
nim_same = ocam_undistort(img, U_same);

end