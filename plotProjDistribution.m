function plotProjDistribution()

inputDir = 'G:\matlab\data\direct\gt\D2_011\4\intr\temp';
dirInfo = dir(fullfile(inputDir, '*.txt'));
image_height = 480;
image_width = 640;
Mat3 = [];
for i = 1 : length(dirInfo)
    data = load(fullfile(inputDir, dirInfo(i).name));
    pix = round(data(:,1:2));
    ind = sub2ind([image_height image_width], pix(:,2),pix(:,1));
    mat = zeros(image_height, image_width);
    mat(ind) = data(:,3);
    mat2 = imgaussfilt(mat,2);
    scale = median(mat(mat~=0))/median(mat2(mat2~=0));
    [xMat, yMat] = meshgrid(1:image_width,1:image_height);
    
    mat3 = scale.*mat2;
    Mat3 = [Mat3; mat3'];
    
    figure,subplot(1,2,1);surf(xMat, yMat, mat3);subplot(1,2,2);contour(mat3, 50);axis equal;colorbar;
end

end