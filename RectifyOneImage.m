function [rectImg,xRange,imgg] = RectifyOneImage(img, rectParam)

rectImg = 255*ones(size(img), 'uint8');
chnOfst = size(img,1)*size(img,2);
for iChn = 1:size(img, 3)
    rectImgChn = 255*ones(size(img,1), size(img,2), 'uint8');
    rectImgChn(rectParam.dstInd) = uint8(sum(bsxfun(@times, rectParam.coef, double(img(rectParam.srcInd + (iChn-1)*chnOfst))), 3));
    rectImg(:,:,iChn) = rectImgChn;
end
IM = zeros(size(img,1),size(img,2));
[y1,x1] = ind2sub(size(IM),rectParam.srcInd(:,:,1));
[y2,x2] = ind2sub(size(IM),rectParam.srcInd(:,:,2));
[y3,x3] = ind2sub(size(IM),rectParam.srcInd(:,:,3));
[y,x] = ind2sub(size(IM),rectParam.dstInd);
xy = [[x1 x2 x3]' [y1 y2 y3]'];
xy = intersect(xy,xy,'rows');
IM = zeros(size(img,1),size(img,2));
ind = sub2ind(size(IM),xy(:,2),xy(:,1));
IM(ind) = 1;
minX = min(xy(:,1));
maxX = max(xy(:,1));
xRange = [minX maxX];
% end


% % % % % % % % [ppx,ppy]=ind2sub(size(imgg),ind)
% % % % % % % % 
% % % % % % % % 
if 1
imgg = zeros(size(img,1),size(img,2));
indd = rectParam.dstInd;
imgg(indd) = 1;
% % figure,imshow(imgg);
end
% % % % % % % % 
% % % % % % % % 
% % % % % % % % for i = 1 : size(img,1)
% % % % % % % %     temp = imgg(i,:)
% % % % % % % %     find
% % % % % % % %     
% % % % % % % %     
% % % % % % % % end
% % % % % % % % 
% % % % % % % % % end
% % % % % % % % 
% % % % % % % % 
% % % % % % % % ind = [];
% % % % % % % % % % % % imgg = zeros(size(img,1),size(img,2));
% % % % % % % % % % % % 
% % % % % % % % % % % % a=rectParam.dstInd; 
% % % % % % % % % % % % b = diff(a); b = [median(b);b];
% % % % % % % % % % % % aa = a(end:-1:1);
% % % % % % % % % % % % bb = diff(aa); bb = [median(bb); bb];
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % c = sort(a,'ascend'); 
% % % % % % % % % % % % d = diff(c);d = [1;d];
% % % % % % % % % % % % cc = sort(a, 'descend');
% % % % % % % % % % % % dd = diff(cc);dd = [1;dd];
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % iind = [a(b~=median(b));c(d~=1);aa(bb~=median(bb));cc(dd~=-1) ];
% % % % % % % % % % % % 
% % % % % % % % % % % % imgg(iind) = 1;
% % % % % % % % % % % % figure,imshow(imgg);
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % [ppx,ppy]=ind2sub(size(imgg),iind);
% % % % % % % % % % % % figure,imshow(imgg);hold on;plot(ppy,ppx,'*');
% % % % % % % % % % % % pix = [ppy ppx];
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % [ia, ib, ic] = unique(pix(:,1));
% % % % % % % % % % % % 
% % % % % % % % % % % % % % % bw = edge((imgg),'canny',[0 0.00000001]);
% % % % % % % % % % % % imggg = zeros(size(img,1),size(img,2));
% % % % % % % % % % % % imggg(iind) = 1;
% % % % % % % % % % % % figure,imshow(imggg);
% % % % % % % % % % % % % % % % % storeDist = [];
% % % % % % % % % % % % % % % % % for i = 1 : length(ia)
% % % % % % % % % % % % % % % % %     tt = ia(i,:);
% % % % % % % % % % % % % % % % % % % % for ii = 1 : size(imggg,2)
% % % % % % % % % % % % % % % % % %     temp = imgg(:,i);
% % % % % % % % % % % % % % % % %     idx = find(pix(:,1) == tt);
% % % % % % % % % % % % % % % % %     if length(idx) == 2
% % % % % % % % % % % % % % % % %         yDist = abs(pix(idx(1),2)-pix((2),2));
% % % % % % % % % % % % % % % % %         storeDist = [storeDist; [yDist pix(idx(1),:) pix(idx(2),:)]];
% % % % % % % % % % % % % % % % %     end
% % % % % % % % % % % % % % % % %         
% % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % end
% % % % % % % % % % % % 
% % % % % % % % % % % % storeDist = [];
% % % % % % % % % % % % for i = 1 : length(ia)
% % % % % % % % % % % %     tt = ia(i,:);
% % % % % % % % % % % % % % % for ii = 1 : size(imggg,2)
% % % % % % % % % % % % %     temp = imgg(:,i);
% % % % % % % % % % % %     idx = find(pix(:,1) == tt);
% % % % % % % % % % % %     if length(idx) == 2
% % % % % % % % % % % %         yDist = abs(pix(idx(1),2)-pix(idx(2),2));
% % % % % % % % % % % %         if yDist > 0.5*size(imggg,1)
% % % % % % % % % % % %         storeDist = [storeDist; [yDist pix(idx(1),:) pix(idx(2),:)]];
% % % % % % % % % % % %         end
% % % % % % % % % % % %     end
% % % % % % % % % % % %        idx  =[]; 
% % % % % % % % % % % % 
% % % % % % % % % % % % end
% % % % % % % % % % % % 
% % % % % % % % % % % % idddxx = find(storeDist(:,1) == min(storeDist(:,1)));
% % % % % % % % % % % % 
% % % % % % % % % % % % closePix = storeDist(idddxx,2:end);
% % % % % % % % % % % % 
% % % % % % % % % % % % yPix = [min(median(closePix(:,2)), median(closePix(:,4))) max(median(closePix(:,2)), median(closePix(:,4)))];
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % [ia2, ib2, ic2] = unique(pix(:,2));
% % % % % % % % % % % % 
% % % % % % % % % % % % storeDist2 = [];
% % % % % % % % % % % % for i = 1 : length(ia2)
% % % % % % % % % % % %     tt = ia2(i,:);
% % % % % % % % % % % % % % % for ii = 1 : size(imggg,2)
% % % % % % % % % % % % %     temp = imgg(:,i);
% % % % % % % % % % % %     idx = find(pix(:,1) == tt);
% % % % % % % % % % % %     if length(idx) == 2
% % % % % % % % % % % %         xDist = abs(pix(idx(1),1)-pix(idx(2),1));
% % % % % % % % % % % %         if xDist>0.5*size(imggg,2)
% % % % % % % % % % % %         
% % % % % % % % % % % %         storeDist2 = [storeDist2; [xDist pix(idx(1),:) pix(idx(2),:)]];
% % % % % % % % % % % %         end
% % % % % % % % % % % %     end
% % % % % % % % % % % %        idx  =[]; 
% % % % % % % % % % % %  
% % % % % % % % % % % % 
% % % % % % % % % % % % end
% % % % % % % % % % % % 
% % % % % % % % % % % % idddxx = find(storeDist2(:,1) == min(storeDist2(:,1)));
% % % % % % % % % % % % 
% % % % % % % % % % % % closePix2 = storeDist2(idddxx,2:end);
% % % % % % % % % % % % 
% % % % % % % % % % % % xPix = [min(median(closePix2(:,2)), median(closePix2(:,4))) max(median(closePix2(:,2)), median(closePix2(:,4)))];
% % % % % % % % % % % % 
% % % % % % % % % % % % corner = [xPix(1) yPix(1);xPix(2) yPix(1);xPix(2) yPix(2);xPix(1) yPix(2)]  %;xPix(1) yPix(1)];
% % % % % % % % % % % % 
% % % % % % % % % % % % figure,imshow(imgg);hold on;plot(corner(:,1),corner(:,2),'-*r');
% % % % % % % % % % % % crop = [xPix; yPix];
% % % % % % % % % % % % 
% % % % % % % % % % % % imgggg = imgg(min(corner(:,2)):max(corner(:,2)), min(corner(:,1)):max(corner(:,1)));
end


% % % end
