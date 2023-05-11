function genDataset()
inputDir = 'D:\Auto\simu\output';

dirInfo = dir(fullfile(inputDir, 'Replay*.mat'));

intrMat = [329.115520046, 0,             320.0
 0,             329.115520046, 240.0
 0,             0,                 1]; 

baseline = 0.5;
[xGrid,yGrid] = meshgrid(1:640, 1:480);
Pix = [xGrid(:) yGrid(:)];
metricPrevPtCcsGT = intrMat\HomoCoord(Pix',1);
metricPrevPtCcsGT = normc(metricPrevPtCcsGT);
scaleMat = reshape(metricPrevPtCcsGT(3,:), 480, 640);
for i = 1 : length(dirInfo)
    clear 'data';
    load(fullfile(inputDir, dirInfo(i).name))
    img = rgb2gray(data{1});
    depth = scaleMat.*data{3};
    disp = intrMat(1,1)*baseline./depth.*2^16;
    imwrite(img, fullfile(inputDir, sprintf('tsukuba_fluorescent_L_%05d.png',i-1)));
    imwrite(uint16(disp),fullfile(inputDir, sprintf('tsukuba_disparity_L_%05d.png',i-1)),'png','bitdepth',16);
%     imwrite(uint16(disp),fullfile(inputDir, sprintf('tsukuba_disparity_L_%05d.png',i-1)),'png','bitdepth',16);
end


end