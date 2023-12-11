function testCircleDetection()
img = imread('G:\matlab\data\direct\gt\D2_011\4\t\23613886382394.bmp');
px = load('G:\matlab\data\direct\gt\D2_011\4\t\0_points.txt');
px = px + 1;
figure,imshow(img);hold on;plot(px(:,1), px(:,2),'.r')
px_round = round(px);
ind = sub2ind(size(img), px_round(:,2), px_round(:,1));
a = zeros(size(img));
a(ind) = 1;
img_dt = bwdist(a,'euclidean');
figure,imshow(img_dt,[])
mask = img_dt < 15;
figure,imshow(mask)
[L, NUM] = bwlabel(mask, 8);
figure,imshow(L, []);
figure,contour(img_dt,200)



end