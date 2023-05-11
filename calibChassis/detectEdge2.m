function [pix,imggggg] = detectEdge2(img,mask,vec, mask_)



if exist('vec') == 2
    vec = [];
end


mask = ones(size(mask));
if ndims(img) == 3;
    img = rgb2gray(img);
end
if max(img(:))~=1;
    try
        imgg = immultiply(img,uint8(mask));
    catch
        imgg = immultiply(img,(mask));
    end
else
    imgg = immultiply(img,mask);
end
if 0
    bw = edge(imgg,'canny',vec);
else
    bw = edge(imgg,'Sobel',[],'vertical');
end

h1 = fspecial('average',[1 1]);
imgfs = imfilter(double(bw),h1);


imgfs(imgfs(:)>0)=1;
imgfs(imgfs(:)<=0)=0;


% leftover = immultiply(img,imgfs);
indd=find(imgfs(:)~=0);
[yy,xx] = ind2sub(size(img),indd);
pix = [xx yy];
thr = 5;

pix = pix(pix(:,1)>thr&pix(:,1)<size(img,2)-thr,:);
pix = pix(pix(:,2)>thr&pix(:,2)<size(img,1)-thr,:);


[validY, validX] = ind2sub(size(mask_), find(mask_ > 0));

pix = intersect(pix, [validX validY], 'rows');
imggggg = zeros(size(img));
ing = sub2ind(size(img),pix(:,2),pix(:,1));
imggggg(ing) = 1;


% % figure,imshow(imgfs);


end