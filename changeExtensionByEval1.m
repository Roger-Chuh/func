% inputDir = 'D:\temp5\20171010\1\1';
function changeExtensionByEval1(inputDir)
dirInfo = dir(fullfile(inputDir,'*.png'));
frameNum = length(dirInfo); 1/2;  %/2;
workDir = pwd;
cd (inputDir)
for i = 1:    frameNum  %length(dirInfo)
%    img = imread(fullfile(inputDir,dirInfo(i).name)); 
%    name = str2double(dirInfo(i).name(5:end-4)); 
%    name = name.*100; name = round(name);
%    name = name*10000000 + 1000000000000000000;
%    imwrite(img,strcat(inputDir,'\',sprintf('%19.0f',name),'.png'));



   eval(['!rename' 32 dirInfo(i).name 32 sprintf('imgL_%04d.jpg', i)]);
%    eval(['!rename' 32 dirInfo(i + frameNum).name 32 sprintf('imgR_%04d.png', i)]);

%     eval(['!rename' 32 dirInfo(i).name 32 sprintf('rectR_%04d.png', i)]);


end
cd (workDir)
end