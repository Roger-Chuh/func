function CalcDisp()

close all;

img1 = imread('C:\Users\Administrator\Desktop\left.pgm');
img2 = imread('C:\Users\Administrator\Desktop\right.pgm');





disp = disparity(img1, img2, 'DisparityRange', [80 128]);
disp(disp < 83) = nan;
disp(disp > 125) = nan;
disp_new = imrotate(disp, 90);
[xMat, yMat] = meshgrid(1:640 , 1:480);
pix = [xMat(:) yMat(:)];
Q =  [ 1., 0., 0., -3.0759966080713497e+02; 0., 1., 0., -2.3564654886450342e+02; 0., 0., 0., 2.8904126210912881e+02;0.       0., -8.6211332447092115e+00, 0. ];
dispList = disp_new(:);
Q(4,3) = -Q(4,3);
xx = Q*pextend([pix dispList]');
X = [xx(1,:)./xx(4,:); xx(2,:)./xx(4,:); xx(3,:)./xx(4,:); ]';
depthMap = reshape(X(:,3), size(disp_new));
figure,pcshow(X)
figure,PcShow(X)
figure,subplot(2,2,1);imshow(imrotate(img1,90));title('left image');
subplot(2,2,2);imshow([img1 img2]);title('left - right');
subplot(2,2,3);imshow(disp_new, []);title('disparity');
subplot(2,2,4);imshow(depthMap, []);title('depth');


end