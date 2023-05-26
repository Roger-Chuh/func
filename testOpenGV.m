function testOpenGV()
close all
K = [500 0 320; 0 500 240; 0 0 1];
markerSize = 150; 100; 200;
markerRow = 4;
markerCol = 5;
[xMat, yMat] = meshgrid(markerSize : markerSize : (markerCol - 0) * markerSize, markerSize : markerSize : (markerRow - 0) * markerSize);
marker1 = [xMat(:) yMat(:) zeros(markerRow*markerCol,1)];

Cam2MarkerRot = rodrigues([0.03 0.02 0.01]);% + 0.0010*(rand(1,3)-0.5));
Cam2MarkerRot = rodrigues([0.3 0.2 0.01]+ 0.0010*(rand(1,3)-0.5));
Cam2MarkerTrans = [505 1000 -3000]';

Cam2MarkerRot2 = rodrigues([0.33 0.21 0.013]+ 0.0010*(rand(1,3)-0.5));
Cam2MarkerTrans2 = [500 1001 -3200]';


Cam2Marker = [Cam2MarkerRot Cam2MarkerTrans; 0 0 0 1];
Cam2Marker2 = [Cam2MarkerRot2 Cam2MarkerTrans2; 0 0 0 1];
Marker2Cam = inv(Cam2Marker);
Marker2Cam2 = inv(Cam2Marker2);
[pt2d,pt3d] = TransformAndProject(marker1, K, Marker2Cam(1:3,1:3), Marker2Cam(1:3,4));
[pt2d2,pt3d2] = TransformAndProject(marker1, K, Marker2Cam2(1:3,1:3), Marker2Cam2(1:3,4));
pt2dX = reshape(pt2d(:,1), size(xMat));
pt2dY = reshape(pt2d(:,2), size(yMat));
pt2dX2 = reshape(pt2d2(:,1), size(xMat));
pt2dY2 = reshape(pt2d2(:,2), size(yMat));
for i = 2 : markerCol
    ab1 = norm([pt2dX(1,i-1) - pt2dX(2,i-1), pt2dY(1,i-1) - pt2dY(2,i-1)]);
    ab2 = norm([pt2dX(1,i) - pt2dX(2,i), pt2dY(1,i-1) - pt2dY(2,i)]);
    
    
    ac1 = norm([pt2dX(1,i-1) - pt2dX(3,i-1), pt2dY(1,i-1) - pt2dY(3,i-1)]);
    ac2 = norm([pt2dX(1,i) - pt2dX(3,i), pt2dY(1,i) - pt2dY(3,i)]);
    
    ad1 = norm([pt2dX(1,i-1) - pt2dX(4,i-1), pt2dY(1,i-1) - pt2dY(4,i-1)]);
    ad2 = norm([pt2dX(1,i) - pt2dX(4,i), pt2dY(1,i) - pt2dY(4,i)]);
    
    ratio1 = (ab1/ac1)/((ab1/ad1));
    ratio2 = (ab2/ac2)/((ab2/ad2));
    ratio1 - ratio2
end

for i = 1 : markerCol
    ab1 = norm([pt2dX(1,i) - pt2dX(2,i), pt2dY(1,i) - pt2dY(2,i)]);
    ab2 = norm([pt2dX2(1,i) - pt2dX2(2,i), pt2dY2(1,i) - pt2dY2(2,i)]);
    
    ac1 = norm([pt2dX(1,i) - pt2dX(3,i), pt2dY(1,i) - pt2dY(3,i)]);
    ac2 = norm([pt2dX2(1,i) - pt2dX2(3,i), pt2dY2(1,i) - pt2dY2(3,i)]);
    
    ad1 = norm([pt2dX(1,i) - pt2dX(4,i), pt2dY(1,i) - pt2dY(4,i)]);
    ad2 = norm([pt2dX2(1,i) - pt2dX2(4,i), pt2dY2(1,i) - pt2dY2(4,i)]);
    
    ratio1 = (ab1/ac1)/((ab1/ad1));
    ratio2 = (ab2/ac2)/((ab2/ad2));
    ratio1 = (ab1/ac1)/((ab1/ad1));
    ratio2 = (ab2/ac2)/((ab2/ad2))
    ratio1 - ratio2
end
for i = 2:markerCol
    
   cross_ratio1 = getCrossRation([pt2dX(:,i) pt2dY(:,i)]) ;
   cross_ratio2 = getCrossRation([pt2dX2(:,i) pt2dY2(:,i)]);
   cross_ratio3 = getCrossRation([pt2dX(:,i-1) pt2dY(:,i-1)]);
   cross_ratio4 = getCrossRation([xMat(:,i-1) yMat(:,i-1)]);
  [ (cross_ratio1 - cross_ratio2)  (cross_ratio1 - cross_ratio3) (cross_ratio1 - cross_ratio4)]
  
end
figure,imshow(zeros(480, 640));hold on;plot(pt2d(:,1), pt2d(:,2), 'or');
pt2d_norm = inv(K)*pextend(pt2d');
pt2d_norm2 = NormalizeVector(pt2d_norm');
if 1
    X = opengv('p3p_kneip_ransac',marker1',pt2d_norm2');
    a = TransformAndProject(marker1, K, X(1:3,1:3)', -X(1:3,1:3)'*X(1:3,4));
    R = X(:,1:3);
    t = X(:,4);
else
    temp = K \ [I; ones(1,size(I,2))];
    I_norms = sqrt(sum(temp.*temp));
    I_normalized = temp ./ repmat(I_norms,3,1);
    X = opengv('epnp',P,I_normalized);
    R = X(:,1:3);
    t = X(:,4);
    
end

end
function cross_ratio = getCrossRation(pt)
x1x2 = det([pt(1,1) pt(2,1); pt(1,2) pt(2,2)]);
x3x4 = det([pt(3,1) pt(4,1); pt(3,2) pt(4,2)]);
x1x3 = det([pt(1,1) pt(3,1); pt(1,2) pt(3,2)]);
x2x4 = det([pt(2,1) pt(4,1); pt(2,2) pt(4,2)]);
cross_ratio = x1x2*x3x4/x1x3/x2x4;
end