function debugVINS()

close all


str = '1\';
% str = '2\';
% str = '3\';
% str = '4\';
% str = '5\';
% str = '7\';
str = '9\';
str = '11\';
str = '12\';

str = '10\';


img11 = imread(['C:\Users\Roger\Desktop\virtual-machine-data\vins-fisheye\',str,'up_orig_img.png']);
img22 = imread(['C:\Users\Roger\Desktop\virtual-machine-data\vins-fisheye\',str,'down_orig_img.png']);

img1 = imread(['C:\Users\Roger\Desktop\virtual-machine-data\vins-fisheye\',str,'up_top_img.png']);
n = size(img1,1);
img2 = imread(['C:\Users\Roger\Desktop\virtual-machine-data\vins-fisheye\',str,'down_top_img.png']);
img3 = imread(['C:\Users\Roger\Desktop\virtual-machine-data\vins-fisheye\',str,'up_side_img.png']);
img4 = imread(['C:\Users\Roger\Desktop\virtual-machine-data\vins-fisheye\',str,'down_side_img.png']);

img3_1 = img3(:,1:n);
img3_2 = img3(:,n+1:2*n);
img3_3 = img3(:,2*n+1:3*n);
img3_4 = img3(:,3*n+1:4*n);

figure,imshow(img11);

% figure,imshow([imrotate(img3_3,180);img1;img3_1]);
% figure,imshow([imrotate(img3_4,-90) img1 imrotate(img3_2,90)]);
pix = [220 361];
pix = [1 1];

pix = [0 0];
pix = [0 1];
pix = [220 361];

topCol = size(img1,2);
topRow = size(img1,1);
sideCol = size(img3_1, 2);
sideRow = size(img3_1, 1);

if 0
    figure,subplot(1,5,1);imshow(img1);hold on;plot(pix(1)+1, pix(2)+1,'.r');title('0');
    subplot(1,5,2);imshow(img3_1);hold on;plot(pix(1)+1, pix(2)+1,'.r');title('1');
    subplot(1,5,3);imshow(img3_2);hold on;plot(pix(1)+1, pix(2)+1,'.r');title('2');
    subplot(1,5,4);imshow(img3_3);hold on;plot(pix(1)+1, pix(2)+1,'.r');title('3');
    subplot(1,5,5);imshow(img3_4);hold on;plot(pix(1)+1, pix(2)+1,'.r');title('4');
end

img = uint8(zeros(topCol + sideRow*2, topCol + sideRow*2));

img(1:sideRow, sideRow + 1 : sideRow + topCol) = imrotate(img3_3,180);
img(sideRow+1:sideRow+topCol,1:sideRow) = imrotate(img3_4,-90);
img(sideRow+1:sideRow+topCol, sideRow+1:sideRow+topCol) = img1;
img(sideRow+1:sideRow+topCol, sideRow+topCol+1:end) = imrotate(img3_2,90);
img(sideRow + topCol+1:end, sideRow+1:sideRow+topCol) = img3_1;



if 0
    % 1 based
    pix0 = pix + [sideRow sideRow];
    pix1 = [pix]                                       + [sideRow (sideRow + topCol)];
    pix2 = [pix(2) (topCol+1 - pix(1))]                + [(sideRow + topCol) sideRow];
    pix3 = [(topCol+1 - pix(1)) (sideRow+1 - pix(2))]  + [sideRow 0];
    pix4 = [(sideRow+1 - pix(2)) pix(1)]               + [0 sideRow];
else
    % 0 based
    pix0 = pix                                         + [sideRow sideRow];
    pix1 = [pix]                                       + [sideRow (sideRow + topCol)];
    pix2 = [pix(2) (topCol-1 - pix(1))]                + [(sideRow + topCol) sideRow];
    pix3 = [(topCol-1 - pix(1)) (sideRow-1 - pix(2))]  + [sideRow 0];
    pix4 = [(sideRow-1 - pix(2)) pix(1)]               + [0 sideRow];
end

pix00 = pix0 - [sideRow sideRow];
pix11 = pix1 - [sideRow (sideRow + topCol)];
pix22 = [(topCol-1 -(pix2(2) - sideRow))       (pix2(1) - (sideRow + topCol))];
pix33 = [(topCol-1 - (pix3(1) - sideRow))       (sideRow-1 - pix3(2))];
pix44 = [(pix4(2) -sideRow)      (sideRow-1 - pix4(1))];


figure,subplot(1,5,1);imshow(img1);hold on;plot(pix(1)+1, pix(2)+1,'or');plot(pix00(1)+1, pix00(2)+1,'.g');title('0');
subplot(1,5,2);imshow(img3_1);hold on;plot(pix(1)+1, pix(2)+1,'or');plot(pix11(1)+1, pix11(2)+1,'.g');title('1');
subplot(1,5,3);imshow(img3_2);hold on;plot(pix(1)+1, pix(2)+1,'or');plot(pix22(1)+1, pix22(2)+1,'.g');title('2');
subplot(1,5,4);imshow(img3_3);hold on;plot(pix(1)+1, pix(2)+1,'or');plot(pix33(1)+1, pix33(2)+1,'.g');title('3');
subplot(1,5,5);imshow(img3_4);hold on;plot(pix(1)+1, pix(2)+1,'or');plot(pix44(1)+1, pix44(2)+1,'.g');title('4');

figure,imshow(img);hold on;
plot(pix0(1)+1, pix0(2)+1,'.r');text(pix0(1)+3,pix0(2)+3,num2str(0), 'Color',[1 0 1],'FontSize',15);
plot(pix1(1)+1, pix1(2)+1,'.r');text(pix1(1)+3,pix1(2)+3,num2str(1), 'Color',[1 0 1],'FontSize',15);
plot(pix2(1)+1, pix2(2)+1,'.r');text(pix2(1)+3,pix2(2)+3,num2str(2), 'Color',[1 0 1],'FontSize',15);
plot(pix3(1)+1, pix3(2)+1,'.r');text(pix3(1)+3,pix3(2)+3,num2str(3), 'Color',[1 0 1],'FontSize',15);
plot(pix4(1)+1, pix4(2)+1,'.r');text(pix4(1)+3,pix4(2)+3,num2str(4), 'Color',[1 0 1],'FontSize',15);


end