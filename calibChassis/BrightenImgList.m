function BrightenImgList()
inputDir = 'D:\Temp\20201221\20201221_161029\video\R';


dirInfo = dir(fullfile(inputDir, '*.png'));





for i = 1 : length(dirInfo)
    
    cbImg = (imread(fullfile(inputDir, dirInfo(i).name)));
    cbImg0 = double(cbImg);
    cbImg0(cbImg>110) = 1;
    cbImg0 = 255.*cbImg0./max(cbImg0(:));
    cbImg = uint8(cbImg0);
    
    cbImg = cat(3, cbImg,cbImg,cbImg);
    imwrite(cbImg, fullfile(inputDir, dirInfo(i).name));
end



end