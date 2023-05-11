function MarkerHist()

close all

img = imread('D:\Temp\20201214\marker54.png');


polygonList = [93 16 190 78 195 173 100 132;...
    426 241 443 252 443 291 427 284;...
    501 296 509 301 509 327 501 324;...
    582 347 597 347 598 362 582 362;...
    723 334 736 334 738 349 723 349;...
    906 280 917 276 917 297 906 300;...
    978 252 996 245 996 273 978 278;...
    1204 267 1254 254 1258 300 1206 309];

figure,imshow(img);hold on;plot([polygonList(:,1);polygonList(:,3);polygonList(:,5);polygonList(:,7)], [polygonList(:,2);polygonList(:,4);polygonList(:,6);polygonList(:,8)],'.r')


gain = 1;

for i = 1 : size(polygonList, 1)
    tempImg = uint8(zeros(size(img(:,:,1))));
    tempCnr = reshape(polygonList(i,:), 2,[])';
    
    
    
    tempCnr = [tempCnr; tempCnr(1,:)];
    [~, len] = NormalizeVector(tempCnr(1:end-1,:) - tempCnr(2:end,:));
    if 0
        targetLen = floor(max(len));
    else
        targetLen = 40;
    end
    targetCoord = [1 1; targetLen 1; targetLen targetLen; 1 targetLen];
    if 0
        [tform,inlierpoints1,inlierpoints2] = estimateGeometricTransform( tempCnr(1:end-1,:), targetCoord,'projective');
    else
        [H,Hnorm,inv_Hnorm] = compute_homography(tempCnr(1:end-1,:)', targetCoord');
    end
    a = pflat(H*pextend(targetCoord'));
    [xMat, yMat] = meshgrid(1:targetLen, 1:targetLen);
    xy = [xMat(:) yMat(:)];
    srcCoord = pflat(H*pextend(xy'));
    warped_intensities = interp2(double(img), srcCoord(1,:)', srcCoord(2,:)', 'linear', 0);
    targetImg = uint8(reshape(warped_intensities, targetLen, targetLen));
%         targetImg([1 end],:) = min(min(targetImg(2:end-1, 2:end-1)));
%         targetImg(:, [1 end]) = min(min(targetImg(2:end-1, 2:end-1)));
    [counts,x] = imhist(targetImg,16);
    T = otsuthresh(counts);
    BW = imbinarize(targetImg,T);
    BW([1:2 end-1:end],:) = false; % min(min(targetImg(2:end-1, 2:end-1)));
    BW(:, [1:2 end-1:end]) = false; % min(min(targetImg(2:end-1, 2:end-1)));
    % figure,imshow(BW)
    figure,subplot(1,2,1),imshow([(gain.*targetImg) uint8(255.*BW)]);subplot(1,2,2);hist(gain.*double(targetImg(:)),round(length(warped_intensities)/10));
    
    if 0
        J = imwarp(img,tform);
        figure(1),clf;imshow(tempImg); hold on;
        for j = 1 : 4
            plot(tempCnr(j:j+1,1), tempCnr(j:j+1,2),'-w','LineWidth',5); drawnow;
        end
        saveas(gcf,'tempMask.jpg');
        tempMask = imread('tempMask.jpg');
    end
end

end