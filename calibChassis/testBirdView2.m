function testBirdView2(K,I0)

% close all


% K = eye(3);
% K(1,1) = 2.0003201001638317e+03;
% K(1,3) = 9.6889274597167969e+02;
% K(2,2) = 2.0003201001638317e+03;
% K(2,3) = 5.4882026672363281e+02;


% K = 1000.*[1.6755         0    0.9676
%     0    1.6755    0.5475
%     0         0    0.0010];

focalLength = [K(1,1) K(2,2)];
principalPoint = [K(1,3) K(2,3)];

imageSize = [800 1280];
camIntrinsics = cameraIntrinsics(focalLength,principalPoint,imageSize);
height = 2; 2.1798;
pitch =0; -0.21; 0;
sensor = monoCamera(camIntrinsics,height,'Pitch',pitch);
distAheadOfSensor = 5;130; 300; 130; 300;
spaceToOneSide = 3; 10; 30;
bottomOffset = 1;
% outView = [bottomOffset,distAheadOfSensor,-spaceToOneSide,spaceToOneSide];

% outView = [xmin,xmax,ymin.ymax];
outView = [0,2,-2,2];

outImageSize = [5*800 1280];imageSize; % [NaN,250];
% outImageSize = [1080 1920];imageSize; % [NaN,250];
birdsEye = birdsEyeView(sensor,outView,outImageSize);

% 
pix = [266 329;1210 132;1797 822;1277 584];

if 0
%     I0 = imread('D:\Temp\20200829\imagel5400_081500_rectify.jpg');
    I = ImgTransform(rgb2gray(I0), K, rotx(0.21));
    pixRot = Orig2Rect(pix, K, K, rotx(0.21),zeros(5,1));
else
%     I0 = imread('E:\bk_20180627\SLAM\slam_algo\L123.png');
    I = ImgTransform((I0), K, rotx(0.0));
    pixRot = Orig2Rect(pix, K, K, rotx(0.0),zeros(5,1));
end

% figure,subplot(1,2,1),imshow(I0);hold on;plot(pix(:,1),pix(:,2),'or');
% subplot(1,2,2),imshow(I);hold on;plot(pixRot(:,1),pixRot(:,2),'or');


% figure
% imshow(I)
% title('Original image')
B = transformImage(birdsEye,I);
% imagePoint = vehicleToImage(birdsEye,[2000,0]);
% annotatedB = insertMarker(B,imagePoint - 5);
% annotatedB = insertText(annotatedB,imagePoint,'20 meters');
figure,subplot(1,2,1),imshow(I);
% subplot(1,2,2);imshow(annotatedB)
subplot(1,2,2);imshow(B)
title('Bird''s eye view')







end