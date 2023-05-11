function testHomo()

close all

K = [499 0 321; 0 501 239; 0 0 1];

if 0
    K_err = [499 0 320; 0 501 239; 0 0 1] + [10 0 -10; 0 -10 10; 0 0 0];
elseif 0
        K_err = K./1.5;
        K_err(1,3) = K(1,3);
        K_err(2,3) = K(2,3);
        K_err(3,3) = 1;
else
    K_err = K;
end

kc = [0.21;-0.1;-0.03;0.02;0];
% kc = [0 0 0 0 0]';
if 0
    kc_err = [0.01 -0.02 0.04 -0.001 0]';
else
    kc_err = zeros(5,1);
end
% kc = [0 0 0 0 0]';

[xMat, yMat] = meshgrid(50*(1:10),50*(1:8));
xyz = [xMat(:) yMat(:)];
xyz(:,3) = 1;


pose = [rodrigues([0.3 0.2 -0.3]) [-300 -150 500]';0 0 0 1];

xyz2 = (pose*[pextend(xyz')])';


pt2d = pflat(K*xyz2(:,1:3)');

%完全无噪声的角点
pt2d22 = pt2d(1:2,:)';

% 对内参和畸变施加噪声(认为是opencv直接检测到的角点)
% pt2d2_noise = remapRect(pt2d22', K_err, K,kc, eye(3));
% 上面写错了，应该只对畸变施加噪声，内参仍是golden，(pt2d2_obs认为是opencv直接检测到的角点)
pt2d2_obs = remapRect(pt2d22', K, K, kc, eye(3));

%用错误的内参和正确的畸变做去畸变
% pt2d2 = Orig2Rect(pt2d2_noise, K, K_err, eye(3), kc);
pt2d2 = Orig2Rect(pt2d2_obs, K_err, K_err, eye(3), kc + kc_err);

xMat2d = reshape(pt2d2(:,1),size(xMat));
yMat2d = reshape(pt2d2(:,2),size(xMat));

imageSize = [480 640];

figure,imshow(zeros(imageSize));hold on;plot(pt2d2(:,1), pt2d2(:,2),'-og');
hori = [];
for i = 1 : size(xMat,1)
    hori_line = [xMat2d(i,:)' yMat2d(i,:)'];
    hori_line(:,3) = 1;
    [~,~,V] = svd(hori_line);
    line_para = V(:,3)';
    line_para = sign(line_para(3)).*line_para;
    line_para = line_para./norm(line_para(1:2));
    line_para1 = line_para./line_para(3);
    hori = [hori; line_para];
    hori_line_fit_err(i,:) = dot(repmat(line_para', 1 ,size(hori_line,1)), hori_line');
    points = lineToBorderPoints(line_para,imageSize) ;
    line(points([1,3]),points([2,4]),'Color',[1 0 0]);
end

vert = [];
for i = 1 : size(xMat,2)
    vert_line = [xMat2d(:,i) yMat2d(:,i)];
    vert_line(:,3) = 1;
    [~,~,V] = svd(vert_line);
    line_para = V(:,3)';
    
    line_para = sign(line_para(3)).*line_para;
    line_para = line_para./norm(line_para(1:2));
    
    
    vert_line_fit_err(i,:) = dot(repmat(line_para', 1 ,size(vert_line,1)), vert_line');
    line_para1 = line_para./line_para(3);
    vert = [vert; line_para];
    
    points = lineToBorderPoints(line_para,imageSize) ;
    line(points([1,3]),points([2,4]),'Color',[0 0 1]);
    
    
end

[~,~,hori_van] = svd(hori);
hori_van = hori_van(:,3);
hori_van = hori_van./norm(hori_van(1:2));
hori_van1 = hori_van./hori_van(3);
hori_van_dir_err = sum(abs(dot(repmat(hori_van1,1,size(hori,1)), hori')));
% hori_van_fir_err = ((dot(repmat(hori_van1,1,size(hori,1)), hori')))

[~,~,vert_van] = svd(vert);
vert_van = vert_van(:,3);
vert_van = vert_van./norm(vert_van(1:2));
vert_van1 = vert_van./vert_van(3);
vert_van_dir_err = sum(abs(dot(repmat(vert_van1,1,size(vert,1)), vert')));
% vert_van_fir_err = ((dot(repmat(vert_van1,1,size(vert,1)), vert')))


hori_vert_van_dir_err = [hori_van_dir_err vert_van_dir_err];

hori_dir = inv(K)*hori_van1;
hori_dir = hori_dir./norm(hori_dir);
vert_dir = inv(K)*vert_van1;
vert_dir = vert_dir./norm(vert_dir);


hori_dir_err = inv(K_err)*hori_van1;
hori_dir_err = hori_dir_err./norm(hori_dir_err);
vert_dir_err = inv(K_err)*vert_van1;
vert_dir_err = vert_dir_err./norm(vert_dir_err);



err1 = dot(hori_dir, vert_dir);
err2 = dot(hori_dir_err, vert_dir_err);
err2 = CalcDegree(hori_dir_err', vert_dir_err');


hori_vert_van_dir__ang_err = [hori_vert_van_dir_err err2]















pt2d_un = inv(K)*pt2d;







[tform,inlierPtsDistorted,inlierPtsOriginal] = estimateGeometricTransform(pt2d2,xyz(:,1:2),'projective');
[tform2,inlierPtsDistorted2,inlierPtsOriginal2] = estimateGeometricTransform(pt2d_un(1:2,:)',xyz(:,1:2),'projective');


T = tform.T';

T2 = tform2.T';


mapped = pflat(T*pextend(pt2d2'));
mapped2 = pflat(T2*pt2d_un);

end