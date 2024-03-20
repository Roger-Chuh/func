function testAprilTagWarp()

img = imread('G:\matlab\data\direct\gt\D2_011\4\1.bmp');
% img = rgb2gray(imread('G:\matlab\data\direct\gt\D2_011\4\2.png'));

pix = [56 119];
source = [56 119; 73 68; 79 147; 97 90];
% source = [323 793; 704 183; 866 1138; 1325 424];
target = [1 1; 50 1; 1 50; 50 50];

tform = estimateGeometricTransform( source,target,'projective');


a = [[378.486, 353.959]
[398.827, 348.18]
[389.853, 333.531]
[371.427, 340.637]];

b = [[0, 10]; [10, 10]; [10, 0]; [0, 0]];


H = tform.T';


check = H * pextend(source');
check = pflat(check);
check = check(1:2,:)';


[xMat, yMat] = meshgrid(1:50,1:50);

pixList = [xMat(:) yMat(:)];

check2 = inv(H) * pextend(pixList');
check2 = pflat(check2);
check2 = check2(1:2,:)';

warped_intensities = interp2(double(img), check2(:, 1), check2(:, 2), 'cubic', 0);

warped = reshape(warped_intensities, 50, 50);
figure,subplot(1,2,1);imshow(img);hold on; plot(source(:,1), source(:,2),'.r');subplot(1,2,2);imshow(uint8(warped));
check - target
end