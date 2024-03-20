function Mat3 = showErrorDistribution()


% Mat3 = showErrorDistribution();
% Mat4 = showErrorDistribution();
% Mat = [Mat3 Mat4];
% figure,contour(Mat', 100);axis equal; colorbar;





inputDir = 'G:\matlab\data\distrib\1';
inputDir = 'G:\matlab\data\distrib\1_1';
% inputDir = 'G:\matlab\data\distrib\2';
% inputDir = 'G:\matlab\data\distrib\2_1';
% inputDir = 'G:\matlab\data\distrib\3';
% inputDir = 'G:\matlab\data\distrib\x1';
% inputDir = 'G:\matlab\data\distrib\4';
% inputDir = 'G:\matlab\data\distrib\ce';
inputDir = 'G:\matlab\data\distrib\ju';

% inputDir = 'G:\matlab\data\distrib\d1';
inputDir = 'G:\matlab\data\distrib\hao_1';%eq
inputDir = 'G:\matlab\data\distrib\hao_2';%no eq
inputDir = 'G:\matlab\data\distrib\d1_1';%eq
inputDir = 'G:\matlab\data\distrib\d1_2';%no eq
inputDir = 'G:\matlab\data\distrib\fail1_1';
inputDir = 'G:\matlab\data\distrib\hao_3';
inputDir = 'G:\matlab\data\distrib\close_1';
inputDir = 'G:\matlab\data\close\ning\Device_D2HD232328D9000957_2022_08_01_17_56_51\controller_2022-08-01-07-21-35\distrib';
inputDir = 'G:\matlab\data\close\hao\Device_D2HD232328D9000859_2022_08_01_19_10_57\distrib';
inputDir = 'G:\matlab\data\distrib\revisit_1';
inputDir = 'G:\matlab\data\distrib\revisit_1_2';
inputDir = 'G:\matlab\data\distrib\revisit_2_1';
inputDir = 'G:\matlab\data\distrib\revisit_2_2';

inputDir = 'G:\matlab\data\distrib\463';
inputDir = 'G:\matlab\data\distrib\exp\board';
inputDir = 'G:\matlab\data\distrib\exp\916';
inputDir = 'G:\matlab\data\distrib\exp\793';
inputDir = 'G:\matlab\data\distrib\exp\cur';
inputDir = 'G:\matlab\data\distrib\d1_ae';
inputDir = 'G:\matlab\data\direct\gt\D2_011\4\error_distri';

board_num = 3;5;


dirInfo = dir(fullfile(inputDir, '*.txt'));
if board_num == 5 % 20240205
    assert(length(dirInfo) == 20);
    camId = [1:5:20];
else
    assert(length(dirInfo) == 12);%6);
    camId = [1:board_num:12];%6];
end
% camId = [1:5:20];




if 0
    image_width = 1280;640;
    image_height = 1024;480;
else
    image_width = 640;
    image_height = 480;
end

is_rgb_stereo = 1;

if board_num == 3
%     camId = [1:board_num:20];
%    camId = camId(1:2); 
end

Mat3 = [];
for i = 1 : length(camId)
%     a = load(fullfile(inputDir, dirInfo(camId(i)).name));
    a = [];
%     for j = camId(i)+0:camId(i)+4
    for j = camId(i)+0:camId(i)+(board_num-1)
        a = [a; load(fullfile(inputDir, dirInfo(j).name))];
        
    end
    [~, err] = NormalizeVector(a(:,3:4));  
    pix = round(a(:,1:2));
    ind = sub2ind([image_height image_width], pix(:,2),pix(:,1));
    pix = round(a(:,1:2));
    ind = sub2ind([image_height image_width], pix(:,2),pix(:,1));
    mat = zeros(image_height, image_width);
    mat(ind) = err;
    mat2 = imgaussfilt(mat,2);
    scale = median(mat(mat~=0))/median(mat2(mat2~=0));
    [xMat, yMat] = meshgrid(1:image_width,1:image_height);
    
    mat3 = scale.*mat2;
    
    Mat3 = [Mat3; mat3'];
    
    figure,subplot(1,2,1);surf(xMat, yMat, mat3);subplot(1,2,2);contour(mat3, 50);axis equal;colorbar;
    
%     figure,subplot(1,2,1);surface(xMat, yMat, mat);
%     subplot(1,2,2),surf(imresize(xMat,0.5), imresize(yMat,0.5), 2.*imresize(mat,0.5));
end








end