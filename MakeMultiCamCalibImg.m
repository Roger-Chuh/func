function MakeMultiCamCalibImg()

inputDir = 'D:\Matlab\calibration-toolbox-0.9\scripts\demo\images';
outputDir = 'D:/Matlab/calibration-toolbox-0.9/scripts/demo/images2';

dirInfo = dir(fullfile(inputDir, '*.png'));

% %YAML:1.0
% ---
% images:
%    - "D:/Matlab/calibration-toolbox-0.9/scripts/demo/images2/"

fid1 = fopen('D:\Matlab\calibration-toolbox-0.9\scripts\demo\imagelist.yaml','w');
fprintf(fid1,'%%YAML:1.0\n');
fprintf(fid1,'---\n');
fprintf(fid1,'images:\n');


for i = 1 : length(dirInfo)
    imgName = dirInfo(i).name;
    imgNameNew = [num2str(str2double(imgName(1))-1) imgName(2:end)];
    fprintf(fid1,strcat('   - "',outputDir,'/',imgNameNew,'"','\n'));
    img = imread(fullfile(inputDir, dirInfo(i).name));
    imwrite(img, fullfile(outputDir, imgNameNew));
end
fclose(fid1);

end