function debugDepth()

close all;



if 0
    intrMat = [192.0982055664062 0 100; 0 192.0982055664062 100; 0 0 1];
    imgL = imread('C:\Users\Roger\Desktop\facTest3\fisheye\up_top_img.png');
    imgR = imread('C:\Users\Roger\Desktop\facTest3\fisheye\down_top_img.png');
    depth = load('C:\Users\Roger\Desktop\facTest3\fisheye\depth.txt');
    [xMat, yMat] = meshgrid(1:200, 1:200);
else
    intrMat = [960.4910888671875  0 500; 0 960.4910888671875  500; 0 0 1];
    imgL = imread('C:\Users\Roger\Desktop\facTest3\fisheye\1\up_top_img.png');
    imgR = imread('C:\Users\Roger\Desktop\facTest3\fisheye\1\down_top_img.png');
    depth = load('C:\Users\Roger\Desktop\facTest3\fisheye\1\depth.txt');
    [xMat, yMat] = meshgrid(1:1000, 1:1000);
end
baseline = 0.064;

disp = intrMat(1,1)*baseline./depth;

dispMat = disparity(imgL, imgR,'DisparityRange',[0 128]);
dispMat(dispMat > 32) = nan;
dispMat(dispMat < 1) = nan;
figure,imshow([disp; dispMat], []);



pix_un = (inv(intrMat)*pextend([xMat(:) yMat(:)]'))';

depthList = depth(:);

valid = find(depthList>0 & depthList < 3);
xyz = repmat(depth(:),1,3).*pix_un;

id = 25287; [xyz(id,1)/xyz(id,3) xyz(id,2)/xyz(id,3)] - pix_un(id,1:2)
figure,pcshow(xyz(valid,:));
end