function UndistortImgList()
inputDir = 'D:\Temp\20181026\2';
paraDir = 'D:\Temp\20181026\calibFishEyeParam2';


dirInfo = dir(fullfile(inputDir, '*.jpg'));


[~, ind1] = sort([dirInfo(1:length(dirInfo)).datenum], 'ascend');

I = [ind1];

tempInfo = cell(length(I),1);
for fg = 1:length(I)
    tempInfo{fg,1} = dirInfo(I(fg)).name;
end
for fgg = 1:length(I)
    dirInfo(fgg).name = tempInfo{fgg};
end





load(fullfile(paraDir,'oCamModel.mat'));

U_same = ocam_undistort_map(oCamModel,'OutputView','same');
intrMat_same = U_same.K';

for i = 1 : length(dirInfo)
    
    nim = ocam_undistort(imread(fullfile(inputDir, dirInfo(i).name)),U_same);
    imwrite(nim, fullfile(inputDir,sprintf('img_%03d.png', i)));% ±£¥Ê÷°
end


end