function ChangeByEval(inputDir)

dirL = fullfile(inputDir, 'img_l');
dirR = fullfile(inputDir, 'img_r');

changeExtensionByEval(dirL, 'imgL');
changeExtensionByEval(dirR, 'imgR');
movefile(fullfile(dirL, '*.png'), inputDir);
movefile(fullfile(dirR, '*.png'), inputDir);


end