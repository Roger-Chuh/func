function [Irec,ind_new,ind_1,ind_2,ind_3,ind_4,a1,a2,a3,a4] = rect_index(I,R,f,c,k,alpha,KK_new);
marg = 10;
skipX = 75; skipY = 75;
skipX = 85; skipY = 85;
% % % % skipX = 10; skipY = 10;
if nargin < 5,
    k = [0;0;0;0;0];
    if nargin < 4,
        c = [0;0];
        if nargin < 3,
            f = [1;1];
            if nargin < 2,
                R = eye(3);
                if nargin < 1,
                    error('ERROR: Need an image to rectify');
                    %break;
                end;
            end;
        end;
    end;
end;


if nargin < 7,
    if nargin < 6,
        KK_new = [f(1) 0 c(1);0 f(2) c(2);0 0 1];
    else
        KK_new = alpha; % the 6th argument is actually KK_new
    end;
    alpha = 0;
end;



% Note: R is the motion of the points in space
% So: X2 = R*X where X: coord in the old reference frame, X2: coord in the new ref frame.


if ~exist('KK_new'),
    KK_new = [f(1) alpha_c*fc(1) c(1);0 f(2) c(2);0 0 1];
end;


[nr,nc] = size(I);

Irec = 255*ones(nr,nc);

[mx,my] = meshgrid(1:nc, 1:nr);
px = reshape(mx',nc*nr,1);
py = reshape(my',nc*nr,1);

if 0
    rays = inv(KK_new)*[(px - 1)';(py - 1)';ones(1,length(px))];
else
    rays = inv(KK_new)*[(px - 0)';(py - 0)';ones(1,length(px))];
end

% Rotation: (or affine transformation):

rays2 = R'*rays;

x = [rays2(1,:)./rays2(3,:);rays2(2,:)./rays2(3,:)];


% Add distortion:
xd = apply_distortion(x,k);


% Reconvert in pixels:

px2_ = f(1)*(xd(1,:)+alpha*xd(2,:))+c(1);
py2_ = f(2)*xd(2,:)+c(2);

px2 = f(1)*(xd(1,:)+alpha*xd(2,:))+c(1);
py2 = f(2)*xd(2,:)+c(2);


% Interpolate between the closest pixels:

px_0 = floor(px2);
py_0 = floor(py2);


% % good_points = find((px_0 >= 0) & (px_0 <= (nc-2)) & (py_0 >= 0) & (py_0 <= (nr-2)));
good_points = find((px_0 > 0) & (px_0 <= (nc-2)) & (py_0 > 0) & (py_0 <= (nr-2)));

px2 = px2(good_points);
py2 = py2(good_points);
px_0 = px_0(good_points);
py_0 = py_0(good_points);

alpha_x = px2 - px_0;
alpha_y = py2 - py_0;

a1 = (1 - alpha_y).*(1 - alpha_x);
a2 = (1 - alpha_y).*alpha_x;
a3 = alpha_y .* (1 - alpha_x);
a4 = alpha_y .* alpha_x;

if 0
    ind_1 = px_0 * nr + py_0 + 1;
    ind_2 = (px_0 + 1) * nr + py_0 + 1;
    ind_3 = px_0 * nr + (py_0 + 1) + 1;
    ind_4 = (px_0 + 1) * nr + (py_0 + 1) + 1;
else
    ind_1 = (px_0-1) * nr + py_0 + 0;
    ind_2 = (px_0 + 0) * nr + py_0 + 0;
    ind_3 = (px_0-1) * nr + (py_0 + 1) + 0;
    ind_4 = (px_0 + 0) * nr + (py_0 + 1) + 0;
end
ind_new = (px(good_points)-1)*nr + py(good_points);


Irec(ind_new) = a1 .* I(ind_1) + a2 .* I(ind_2) + a3 .* I(ind_3) + a4 .* I(ind_4);


return;
[py11,px11] = ind2sub(size(I),ind_1);
[py22,px22] = ind2sub(size(I),ind_2);
[py33,px33] = ind2sub(size(I),ind_3);
[py44,px44] = ind2sub(size(I),ind_4);
[xu,yu] = meshgrid(1:nc, 1:nr);

% % imgL = imread(strcat('D:\bk_20180627_\nextvpu\algorithm\output\calibration\cbccalib\stereo_pairs\iEye\hua_01_01','\cbpair_left001.png'));
% % nim = interp2(xu,yu,double(rgb2gray(imgL)),reshape(px2_,nc,nr)',reshape(py2_,nc,nr)');

% % imgL = imread('D:\bk_20180627_\nextvpu\algorithm\output\calibration\cbccalib\stereo_pairs\iEye\i2c15_6\cbpair_left001.bmp');
% % nim = interp2(xu,yu,double(rgb2gray(imgL)),reshape(px2_,nc,nr)',reshape(py2_,nc,nr)');
% % figure,imshow(ones(480,640));hold on;plot(cbcXYLR{1}(1,:),cbcXYLR{1}(2,:),'or');plot(pixDist{1}(:,1),pixDist{1}(:,2),'.g');plot(pixDist{2}(:,1),pixDist{2}(:,2),'.b');legend('orig','rect-mod','rect-orig');




% xd=apply_distortion(x,kc);
% px2__=f(1)*px2_+c(1);
% py2__=f(2)*xd(2,:)+c(2);
Kinit = [f(1) 0 c(1);0 f(2) c(2);0 0 1];
skipNum = 50;

figure(11),plotUndistortion(nc, nr, Kinit, Kinit,k,eye(3),skipNum);
figure(12),pixNewAll = plotUndistortion(nc, nr, Kinit, KK_new,k,R,skipNum);
% % [xuu,yuu] = meshgrid(1:nc/skipNum:nc, 1:nr/skipNum:nr);
[xuu,yuu] = meshgrid(1:skipX:nc, 1:skipY:nr);
[xuuAll,yuuAll] = meshgrid(1:1:nc, 1:1:nr);

pixDist1 = remapRect([xuu(:) yuu(:)]', Kinit, Kinit,k, eye(3));
pixOrig = remapRect([xuu(:) yuu(:)]', KK_new, Kinit,k, eye(3));

pixOrigAll = remapRect([xuuAll(:) yuuAll(:)]', KK_new, Kinit,k, R);
pixOrig = remapRect([xuu(:) yuu(:)]', KK_new, Kinit,k, R);


xMat = reshape(pixOrig(:,1),size(xuu,1),size(xuu,2));
yMat = reshape(pixOrig(:,2),size(xuu,1),size(xuu,2));



dx=pixOrig(:,1)-xuu(:);dx_ = reshape(dx,size(xuu,1),size(xuu,2))';dx_ = dx_(:)';
dy=pixOrig(:,2)-yuu(:);dy_ = reshape(dy,size(xuu,1),size(xuu,2))';dy_ = dy_(:)';
% Q=quiver(px+1,py+1,dx,dy);
% figure,
figure(1),clf; imshow(zeros(nr,nc));hold on;plot(pixOrig(:,1),pixOrig(:,2),'xg');hold on;plot(xuu(:),yuu(:),'.r'); legend('orig','rect');
Q=quiver(xuu(:)+0,yuu(:)+0,dx,dy,0.3);
hold on;
plot(Kinit(1,3),Kinit(2,3),'o');
% plot((nc-1)/2+1,(nr-1)/2+1,'x');
plot(KK_new(1,3),KK_new(2,3),'x');
dr=reshape(sqrt((dx_.*dx_)+(dy_.*dy_)),size(xuu,2),size(xuu,1))';
[C,h]=contour(xuu,yuu,dr,'k');
clabel(C,h);

figure(2),clf;imshow(zeros(nr,nc));hold on;
for i = 1 : size(xMat,1)
    p1 = polyfit(xMat(i,:),yMat(i,:),2);
    if p1(1) > 0
        if min(yMat(i,:)) > marg && min(yMat(i,:)) < nr - marg
            yHeight(i,:) = repmat(min(yMat(i,:)),1,size(xMat,2));
        elseif min(yMat(i,:)) <=marg
            yHeight(i,:) = repmat(0,1,size(xMat,2));
        else
            yHeight(i,:) = repmat(nr,1,size(xMat,2));
        end
    else
        if max(yMat(i,:)) > marg && max(yMat(i,:)) < nr - marg
            yHeight(i,:) = repmat(max(yMat(i,:)),1,size(xMat,2));
        elseif max(yMat(i,:)) < marg
            yHeight(i,:) = repmat(0,1,size(xMat,2));
        else
            yHeight(i,:) = repmat(nr,1,size(xMat,2));
        end
    end
    [~,~,VV] = svd([xMat(i,:);yHeight(i,:);ones(1,size(xMat,2))]');
    initLine = VV(:,end)';
    if initLine(end) ~= 0
        linek = [initLine(1)/initLine(3) initLine(2)/initLine(3) 1];
    else
        linek = initLine;
    end
    LineHori(i,:) = linek;
    if 1
        points2 = lineToBorderPoints(linek, size(I));
        figure(2),line(points2(:, [1,3])', points2(:, [2,4])','Color',[0 1 0]);
        
        figure(1),plot(xMat(i,:),yMat(i,:),'-b');
    end
    %     figure(2),hold on;plot(xMat(i,:),yHeight(i,:),'-c');
end
for i = 1 : size(xMat,2)
    p2 = polyfit(yMat(:,i),xMat(:,i),2);
    if p2(1) > 0
        if min(xMat(:,i)) > marg && min(xMat(:,i)) < nc -marg
            xWidth(:,i) = repmat(min(xMat(:,i)),size(xMat,1),1);
        elseif min(xMat(:,i)) <= marg
            xWidth(:,i) = repmat(0,size(xMat,1),1);
        else
            xWidth(:,i) = repmat(nc,size(xMat,1),1);
        end
    else
        if max(xMat(:,i)) > marg && max(xMat(:,i)) < nc-marg
            xWidth(:,i) = repmat(max(xMat(:,i)),size(xMat,1),1);
        elseif max(xMat(:,i)) <= marg
            xWidth(:,i) = repmat(0,size(xMat,1),1);
        else
            xWidth(:,i) = repmat(nc,size(xMat,1),1);
        end
    end
    
    [~,~,VV] = svd([xWidth(:,i)';yMat(:,i)';ones(1,size(xMat,1))]');
    initLine = VV(:,end)';
    if initLine(end) ~= 0
        linek = [initLine(1)/initLine(3) initLine(2)/initLine(3) 1];
    else
        linek = initLine;
    end
    LineVerti(i,:) = linek;
    if 1
        points2 = lineToBorderPoints(linek, size(I));
        figure(2),line(points2(:, [1,3])', points2(:, [2,4])','Color',[0 1 0]);
        
        
        figure(1),plot(xMat(:,i),yMat(:,i),'-b');
    end
    %     figure(2),hold on;plot(xWidth(:,i),yMat(:,i),'-c');
end



% [C,h]=contour(xuu,yuu,dr,'k');
% clabel(C,h);
% % % figure,quiver(xuu,yuu,reshape(pixDist(:,1),size(xuu,1),size(xuu,2))-xuu, reshape(pixDist(:,2),size(xuu,1),size(xuu,2))-yuu,'r');
pixDist2 = remapRect(pixOrig', Kinit, Kinit,k, eye(3));
% % % % % % figure,imshow(zeros(nr,nc));hold on;plot(pixOrig(:,1),pixOrig(:,2),'.g');hold on;plot(xuu(:),yuu(:),'.r');legend('orig','rect');
% % % % % % figure,imshow(zeros(nr,nc));hold on;plot(pixDist2(:,1),pixDist2(:,2),'.g');hold on;plot(pixOrig(:,1),pixOrig(:,2),'.r');legend('orig','rect');


xPoint = []; yPoint = [];
for i = 1 : size(LineHori,1)
    LineHoriTmp = repmat(LineHori(i,:),size(LineVerti,1),1);
    rowPoint = cross(LineVerti,LineHoriTmp);
    rowPoint = pflat(rowPoint');
    rowPoint = rowPoint(1:2,:)';
    if  1 %mean(rowPoint(:,2)) > 10 && mean(rowPoint(:,2)) < nr-10
        rowPoint = rowPoint(rowPoint(:,1) >= 0-0.0001 & rowPoint(:,1) <= nc + 0.0001,:);
        xPoint = [xPoint;rowPoint(:,1)'];
        yPoint = [yPoint;rowPoint(:,2)'];
% %         figure(2),plot(rowPoint(:,1),rowPoint(:,2),'or','LineWidth',3);
    else
    end
end
xPoint(:,1) = 0;
xPoint(:,end) = nc;
xPoint = round(xPoint);
yPoint(1,:) = 0;
yPoint(end,:) = nr;
yPoint = round(yPoint);
id1 = find(yPoint(:,1) > 0 & yPoint(:,1) < nr);
yPoint = yPoint(id1(1)-1:id1(end)+1,:);
xPoint = xPoint(id1(1)-1:id1(end)+1,:);
id2 = find(xPoint(1,:) > 0 & xPoint(1,:) < nc);
xPoint = xPoint(:,id2(1)-1:id2(end)+1);
yPoint = yPoint(:,id2(1)-1:id2(end)+1);

figure(2),plot(xPoint(:),yPoint(:),'or','LineWidth',3);

pixOrigNew = remapRect([round(xPoint(:)) round(yPoint(:))]', KK_new, Kinit,k, R);
figure,imshow(zeros(nr,nc));hold on;plot(round(xPoint(:)),round(yPoint(:)),'.r');plot(pixOrigNew(:,1),pixOrigNew(:,2),'.g');legend('rect','orig');
xPointOrig = reshape(pixOrigNew(:,1),size(xPoint));
yPointOrig = reshape(pixOrigNew(:,2),size(xPoint));



for i = 1 : size(xPoint,1)
    xPointTmp = xPoint(i,:);
    xPointOrigTmp = xPointOrig(i,:);
    y_ = yPoint(i,:);
    
    params = splineParams(xPointTmp', xPointOrigTmp');
    Params1{i,1} = params;
    [Y, xx] = spline_(xPointTmp', params);
    if 1
        v = splineval(xPointTmp, params, xx);
    else
        coeff = params(:,4:-1:1);
        [Y, xx] = spline_(xPointTmp', params);
        [~,index] = histc(xx',[-inf,xPointTmp(2:end-1),inf]);
        xxNew = xx-xPointTmp(index)';
        
        v = coeff(index,1);
        for ii=2:4
            v = xxNew(:).*v + coeff(index,ii);
        end
        v = reshape(v,[1,length(xx)]);
    end
    x__ = xx; %xPointTmp(1):xPointTmp(end)-1;
    y__ = repmat(y_(1),length(x__),1);
    pixOrigNew__ = remapRect([x__ y__]', KK_new, Kinit,k, R);
    figure(50),clf,plot(xPointTmp', xPointOrigTmp','-x');hold on;plot(xx,Y,'-');plot(x__,pixOrigNew__(:,1));legend('control point','fitted','golden')
    
    
end

for i = 1 : size(xPoint,2)
    yPointTmp = yPoint(:,i);
    xPointTmp_ = xPoint(:,i);
    pixOrigTmp_ = remapRect([xPointTmp_ yPointTmp]', KK_new, Kinit,k, R);
    x_ = xPoint(:,i);
    
    
    yPointOrigTmp = yPointOrig(:,i);
    params = splineParams(yPointTmp, yPointOrigTmp);
    [Y, xx] = spline_(yPointTmp, params);
    v = splineval(yPointTmp', params, xx);
    Params2{i,1} = params;
    %     x__ = xx; %xPointTmp(1):xPointTmp(end)-1;
    y__ = xx;
    x__ = repmat(x_(1),length(y__),1);
    pixOrigNew__ = remapRect([x__ y__]', KK_new, Kinit,k, R);
    
    figure(51),clf;plot(yPointTmp', yPointOrigTmp','-x');hold on;plot(xx,Y,'-');plot(y__,pixOrigNew__(:,2));legend('control point','fitted','golden')
    
    
    % %     figure,plot(yPointTmp, yPointOrigTmp,'-x');hold on; plot(xx,Y,'-');
    
    
end

% [xuuAll,yuuAll] = meshgrid(1:1:nc, 1:1:nr);
% % xuuAll = xPoint;
% % yuuAll = yPoint;

pixOrigAll = remapRect([xuuAll(:) yuuAll(:)]', KK_new, Kinit,k, R);


PixOrigTmp = [];PixNew3 = [];PixRectTmp = [];

for i = 1 : size(xPoint,1) - 1
    
    xTmp = xuuAll(yuuAll(:,1) >= yPoint(i,1) & yuuAll(:,1) < yPoint(i + 1,1),:);
    yTmp = yuuAll(yuuAll(:,1) >= yPoint(i,1) & yuuAll(:,1) < yPoint(i + 1,1),:);
    for j = 1 : size(xPoint,2)-1
        xxTmp = xTmp(:,xTmp(1,:) >= xPoint(1,j) & xTmp(1,:) < xPoint(1,j+1));
        yyTmp = yTmp(:,xTmp(1,:) >= xPoint(1,j) & xTmp(1,:) < xPoint(1,j+1));
        
        
        pixOrigTmp = remapRect([xxTmp(:) yyTmp(:)]', KK_new, Kinit,k, R);
        pixOrigSplineX1 = splineval(xPointTmp, Params1{i}, xxTmp(1,:)');
        pixOrigSplineX2 = splineval(xPointTmp, Params1{i+1}, xxTmp(1,:)');
        ratioX = (xxTmp - xPoint(1,j))./(xPoint(1,j+1) - xPoint(1,j));
        ratioX = ratioX(1,:)';
        
        xxyyTmp = [xxTmp(:) yyTmp(:)];
        
        try
            pixOrigSplineY1 = splineval(yPointTmp', Params2{j}, yyTmp(:,1));
        catch
            asvlk = 1;
        end
        pixOrigSplineY2 = splineval(yPointTmp', Params2{j+1}, yyTmp(:,1));
        ratioY = (yyTmp - yPoint(i,1))./(yPoint(i+1,1) - yPoint(i,1));
        ratioY = ratioY(:,1);
        if 0
            pixOrigSplineX = (1 - ratioY).*pixOrigSplineX1 + ratioY.*pixOrigSplineX2;
            pixOrigSplineY = (1 - ratioX).*pixOrigSplineY1 + ratioX.*pixOrigSplineY2;
        end
        pixNew = [];pixNew2 = [];
        
        
        [pixOrigSplineX1_1, index_x11] = splineval(xPointTmp, Params1{i}, xxyyTmp(:,1));
        [pixOrigSplineX2_1, index_x22]= splineval(xPointTmp, Params1{i+1}, xxyyTmp(:,1));
        [pixOrigSplineY1_1, index_y11] = splineval(yPointTmp', Params2{j}, xxyyTmp(:,2));
        [pixOrigSplineY2_1, index_y22] = splineval(yPointTmp', Params2{j+1}, xxyyTmp(:,2));
        ratioXX_ = (xxyyTmp(:,1) - xPointTmp(index_x11))./(xPointTmp(index_x11+1)-xPointTmp(index_x11));
        ratioYY_ = (xxyyTmp(:,2) - yPointTmp(index_y11))./(yPointTmp(index_y11+1)-yPointTmp(index_y11));
        
        xInterp_ = (1 - ratioYY_).*pixOrigSplineX1_1 + ratioYY_.*pixOrigSplineX2_1;
        yInterp_ = (1 - ratioXX_).*pixOrigSplineY1_1 + ratioXX_.*pixOrigSplineY2_1;
        pixNew3 = [xInterp_ yInterp_];
        
        
        if 0
            for kk = 1 : size(xxyyTmp,1)
                xxyyTmp_ = xxyyTmp(kk,:);
                [pixOrigSplineX1_, index_x1] = splineval(xPointTmp, Params1{i}, xxyyTmp_(1)');
                [pixOrigSplineX2_, index_x2]= splineval(xPointTmp, Params1{i+1}, xxyyTmp_(1)');
                [pixOrigSplineY1_, index_y1] = splineval(yPointTmp', Params2{j}, xxyyTmp_(2));
                [pixOrigSplineY2_, index_y2] = splineval(yPointTmp', Params2{j+1}, xxyyTmp_(2));
                pt1 = [pixOrigSplineX1_ pixOrigSplineY1_];
                pt2 = [pixOrigSplineX1_ pixOrigSplineY2_];
                pt3 = [pixOrigSplineX2_ pixOrigSplineY1_];
                pt4 = [pixOrigSplineX2_ pixOrigSplineY2_];
                ratio1 = [];
                ratioXX = (xxyyTmp_(1) - xPointTmp(index_x1))./(xPointTmp(index_x1+1)-xPointTmp(index_x1));
                ratioYY = (xxyyTmp_(2) - yPointTmp(index_y1))./(yPointTmp(index_y1+1)-yPointTmp(index_y1));
                
                xInterp = (1 - ratioYY).*pixOrigSplineX1_ + ratioYY.*pixOrigSplineX2_;
                yInterp = (1 - ratioXX).*pixOrigSplineY1_ + ratioXX.*pixOrigSplineY2_;
                
                XX = (xxyyTmp_(1)-xPointTmp(index_x1))/(xPointTmp(index_x1+1)-xPointTmp(index_x1));
                YY = (xxyyTmp_(2)-yPointTmp(index_y1))/(yPointTmp(index_y1+1)-yPointTmp(index_y1));
                XNew = [1-XX XX]*[pt1(1) pt2(1);pt3(1) pt4(1)]*[1-YY;YY];
                YNew = [1-XX XX]*[pt1(2) pt2(2);pt3(2) pt4(2)]*[1-YY;YY];
                
                pixNew(kk,:) = [xInterp yInterp];  %[XNew YNew];
                pixNew2(kk,:) = [XNew YNew];
                %            plot(xInterp,yInterp,'.g');
                %            drawnow;
            end
        end
%         figure(100),clf;imshow(zeros(1080,1920));hold on;plot(pixOrigTmp(:,1),pixOrigTmp(:,2),'or'); plot(pixNew3(:,1),pixNew3(:,2),'.g');
%         drawnow;
%         figure(101),plot(pixOrigTmp(:,1)-pixNew3(:,1),pixOrigTmp(:,2)-pixNew3(:,2),'+r');drawnow
        %         pixOrigSpline = [pixOrigSplineX,pixOrigSplineY];
%         [pixOrigSplineGridX, pixOrigSplineGridY] = meshgrid(pixOrigSplineX,pixOrigSplineY);
        
        PixOrigTmp = [PixOrigTmp;pixOrigTmp];
        PixNew3 = [PixNew3;pixNew3];
        PixRectTmp = [PixRectTmp; [xxTmp(:) yyTmp(:)]];
        
        if 0
            
            aa = intersect(pixOrigTmp(:,1),pixOrigTmp(:,1));
            figure,imshow(zeros(1080,1920));hold on;plot(aa,repmat(1,length(aa),1),'or');plot(pixOrigSplineX,repmat(1,length(pixOrigSplineX),1),'.g')
            bb = intersect(pixOrigTmp(:,2),pixOrigTmp(:,2));
            figure,imshow(zeros(1080,1920));hold on;plot(repmat(1,length(bb),1),bb,'or');plot(repmat(1,length(pixOrigSplineY),1),pixOrigSplineY,'.g')
            figure,imshow(zeros(1080,1920));hold on;plot(pixOrigTmp(:,1),pixOrigTmp(:,2),'or');plot(pixOrigSplineGridX(:),pixOrigSplineGridY(:),'.g');
        end
        
    end
end
    idShow = randperm(size(PixNew3,1));
    showNum = 10000;
    figure(100),clf;imshow(zeros(1080,1920));hold on;plot(PixOrigTmp(idShow(1:showNum),1),PixOrigTmp(idShow(1:showNum),2),'or'); plot(PixNew3(idShow(1:showNum),1),PixNew3(idShow(1:showNum),2),'.g');
    figure(101),plot(PixOrigTmp(idShow(1:showNum),1)-PixNew3(idShow(1:showNum),1),PixOrigTmp(idShow(1:showNum),2)-PixNew3(idShow(1:showNum),2),'+r');
    [~, projErr] = NormalizeVector(PixOrigTmp-PixNew3);
    idErr = find(projErr > 0.25);
    figure(2),plot(PixOrigTmp(idErr,1),PixOrigTmp(idErr,2),'or'); plot(PixNew3(idErr,1),PixNew3(idErr,2),'.g');
%     [~,~,idDraw] = intersect(PixOrigTmp,pixOrigAll,'rows');
    imgL = rgb2gray(imread('D:\bk_20180627_\nextvpu\algorithm\output\calibration\cbccalib\stereo_pairs\iEye\zed2\imgL__left0001.png'));
    imgL1 = zeros(nr,nc);
    PixNew3_ = round(PixNew3);
    idIn = PixNew3_(:,1) >= 1 & PixNew3_(:,1) <= nc & PixNew3_(:,2) >= 1 & PixNew3_(:,2) <= nr;
%     idDraw_ = idDraw(idIn);
    PixNew3_ = PixNew3_(idIn,:);
    PixRectTmp_ = PixRectTmp(idIn,:);
    ind1 = sub2ind(size(imgL1),PixNew3_(:,2),PixNew3_(:,1));
    ind2 = sub2ind(size(imgL1),PixRectTmp_(:,2),PixRectTmp_(:,1));
    imgL1(ind2) = imgL(ind1);
return;


% Convert in indices:

fact = 3;

[XX,YY]= meshgrid(1:nc,1:nr);
[XXi,YYi]= meshgrid(1:1/fact:nc,1:1/fact:nr);

%tic;
Iinterp = interp2(XX,YY,I,XXi,YYi);
%toc

[nri,nci] = size(Iinterp);


ind_col = round(fact*(f(1)*xd(1,:)+c(1)))+1;
ind_row = round(fact*(f(2)*xd(2,:)+c(2)))+1;

good_points = find((ind_col >=1)&(ind_col<=nci)&(ind_row >=1)& (ind_row <=nri));
end

function pixNew = plotUndistortion(nc, nr, Kinit, kk_new,kc,R,skipNum)


[mx,my] = meshgrid(1:nc/skipNum:(nc-0),1:nr/skipNum:(nr-0));
[nnx,nny]=size(mx);
px=reshape(mx',nnx*nny,1);
py=reshape(my',nnx*nny,1);
% % kk_new=[fc(1) alpha_c*fc(1) cc(1);0 fc(2) cc(2);0 0 1];
rays=inv(kk_new)*[px';py';ones(1,length(px))];
rays = R'*rays;


x=[rays(1,:)./rays(3,:);rays(2,:)./rays(3,:)];


title2=strcat('Complete Distortion Model');

% % fh1 = 2;

%if ishandle(fh1),
%    close(fh1);
%end;
% % % % % % figure; clf;
xd=apply_distortion(x,kc);
px2=Kinit(1,1)*xd(1,:)+Kinit(1,3);
py2=Kinit(2,2)*xd(2,:)+Kinit(2,3);
pixNew = [px2 py2];
dx=px2'-px;
dy=py2'-py;
% Q=quiver(px+1,py+1,dx,dy);
Q=quiver(px+0,py+0,dx,dy);
hold on;
plot(Kinit(1,3),Kinit(2,3),'o');
plot((nc-1)/2+1,(nr-1)/2+1,'x');
dr=reshape(sqrt((dx.*dx)+(dy.*dy)),nny,nnx)';
[C,h]=contour(mx,my,dr,'k');
clabel(C,h);
Mean=mean(mean(dr));
Max=max(max(dr));
title(title2);

end
function pixDist = remapRect(pixRect, KRect, KUndist,distCoeff, RL)

alpha = 0;

rays = inv(KRect)*pextend(pixRect);


% Rotation: (or affine transformation):

rays2 = RL'*rays;

x = [rays2(1,:)./rays2(3,:);rays2(2,:)./rays2(3,:)];


% Add distortion:
xd = apply_distortion(x,distCoeff);


% Reconvert in pixels:

px2_ = KUndist(1,1)*(xd(1,:) + alpha*xd(2,:)) + KUndist(1,3);
py2_ = KUndist(2,2)*xd(2,:) + KUndist(2,3);
pixDist = [px2_;py2_]';

end


function [v, indexx] = splineval(xPointTmp, params, xx)
coeff = params(:,4:-1:1);
% %     [Y, xx] = spline_(xPointTmp', params);
[~,index] = histc(xx',[-inf,xPointTmp(2:end-1),inf]);
xxNew = xx-xPointTmp(index)';

v = coeff(index,1);

indexx = unique(index);
for ii=2:4
    v = xxNew(:).*v + coeff(index,ii);
end
v = reshape(v,[1,length(xx)])';
end
