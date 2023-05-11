function testGestureUndistortion()

close all


deOutlierThr = 1; 2;  1;
num_levels = 1; 2; 4; 2; 3; 2; 4; 3; 4; 2; 3; 1; 2; 4;  6; 4; 2;
areaRatio = 0.08; 0.001; 0.08; 0.25; 0.08; 0.005;
initial_sigma = 5;
default_dof = 5;
para = 5;


inputDir0 = 'G:\matlab\data\gesture\sr305';
inputDepthDir = 'G:\matlab\data\gesture\sr305\Depth_0';
dirInfo0 = dir(fullfile(inputDir0, '0*'));
dirInfo = dir(fullfile(inputDir0,dirInfo0(1).name,'*.bmp'));
dirDepthInfo = dir(fullfile(inputDepthDir, '*.png'));

useRS = false; true;


camParam = [473.5599472 472.74264 326.1020203 250.0337029 -0.1483133073 -0.008184294911 -0.001698424582 -0.0004613224404 0];
camParam = [475.399, 475.399 318.445, 245.336 [-0.162634913921236 0.012270673412312 -0.002178680927890 -0.002448749761278 0.003953800470367]];
if useRS
    camParam = [475.399, 475.399 318.445, 245.336 [0.153221,0.101447,0.00371338,0.00427771,0.002554]];
%     camParam = [475.399, 475.399 318.445, 245.336 [0.153221,0.101447,0.00427771,0.002554, 0.00371338]];
end
if 0
    camParams = [235.9157615 235.8859066 319.1806366 239.8611012 0.2354646331 -0.2694915482 0.2146688075 -0.08427136 0 ;
        235.9157615  235.8859066 319.1806366 239.8611012 0.2354646331 -0.2694915482 0.2146688075 -0.08427136 0;
        235.9157615 235.8859066 319.1806366 239.8611012 0.2354646331 -0.2694915482 0.2146688075 -0.08427136 0;
        237.8438076 237.376296 321.4993816  238.6343759 0.2190979312 -0.1758694847 0.06978008689 -0.0123121619 0;
        235.9157615 235.8859066 319.1806366  239.8611012 0.2354646331 -0.2694915482 0.2146688075 -0.08427136 0 ;
        235.9157615 235.8859066 319.1806366 239.8611012 0.2354646331 -0.2694915482 0.2146688075 -0.08427136 0  ];
    
    camExs = {[0.8647022242 0.07272610075 0.4969919293 -359.7388473
        -0.2105544808 0.9508131599 0.227202873 -55.5349488
        -0.4560228877 -0.3011067073 0.8374830605 322.3314639
        0 -0 -0 1 ];
        [0.769118293 0.01756819909 0.6388649386 -391.5309403
        -0.2473821671 0.9298833411 0.2722481137 -75.2212005
        -0.5892869545 -0.3674347975 0.7195363471 330.0818662
        0 -0 -0 1 ];
        [0.769118293 0.01756819909 0.6388649386 -391.5309403
        -0.2473821671 0.9298833411 0.2722481137 -75.2212005
        -0.5892869545 -0.3674347975 0.7195363471 330.0818662
        0 -0 -0 1 ];
        [-0.5263033604 -0.8492068116 0.04304142163 356.9529387
        0.8439264581 -0.5155060037 0.1484644519 -212.2669178
        -0.1038889126 0.1144611345 0.9879806388 3.793192368
        -0 -0 -0 1 ];
        [0.769118293 0.01756819909 0.6388649386 -391.5309403
        -0.2473821671 0.9298833411 0.2722481137 -75.2212005
        -0.5892869545 -0.3674347975 0.7195363471 330.0818662
        0 -0 -0 1 ];
        [0.769118293 0.01756819909 0.6388649386 -391.5309403
        -0.2473821671 0.9298833411 0.2722481137 -75.2212005
        -0.5892869545 -0.3674347975 0.7195363471 330.0818662
        0 -0 -0 1 ]};
else
    camParams = [233.2686658 232.1117792 324.6003619  238.8964722  0.2328129411 -0.1823894791 0.06185659214 -0.007427229484 0 ;
        234.9885021 235.0691697 318.2176619  239.8823665  0.2174298565 -0.1784846445 0.07801622566 -0.01845956545 0;
        235.7166187 235.883154 321.7146205  240.0173367 0.2070347798 -0.1817326374 0.07896338097 -0.01245081646 0;
        237.9928022 236.7035253 321.0940689  238.9517572 0.2049304802 -0.1494268931 0.05150551133 -0.007911991223 0;
        234.0699835 233.6192187 320.3657268 239.6227658 0.2171037556 -0.1632187169 0.0674911462 -0.01382833219 0 ;
        235.0028935 235.237877 319.6477546 238.0883737 0.2109787705 -0.1694305837 0.0774062931 -0.0184859387 0  ];
    
    camExs = {[-0.8710770679 0.1059687321 0.4795783248 -110.391727
        -0.06570125205 0.9425296005 -0.3275992943 66.76587196
        -0.4867320487 -0.3168731291 -0.8140536424 633.8147144
        0 0 -0 1  ];
        [0.905291581 -0.2593245833 -0.3364489766 -91.2387976
        -0.08693893084 -0.8883780898 0.4508059359 -173.1113887
        -0.4157989606 -0.3788603042 -0.8267866075 692.0050067
        -0 0 -0 1 ];
        [0.7827079685 0.02782273135 0.6217669432 -396.1236851
        -0.2449329143 0.9321595666 0.2666203477 -74.36695669
        -0.572167898 -0.3609770602 0.7364234234 323.9959838
        0 -0 -0 1 ];
        [ -0.5248414536 -0.8491459491 0.05909827129 355.4439856
        0.8463125496 -0.5131391613 0.1429799618 -214.0782889
        -0.091085218 0.1250574196 0.9879595765 0.1829828146
        -0 -0 -0 1  ];
        [-0.3744915123 -0.05056446071 0.9258506049 -114.954204
        -0.5737758808 0.7970120895 -0.1885549461 281.2524225
        -0.728379946 -0.6018429732 -0.3274869307 578.737508
        0 0 -0 1 ];
        [0.2314769276 -0.9301545684 -0.2850103697 80.31656087
        0.4592273162 -0.1537981841 0.8749036465 -429.7652873
        -0.857629701 -0.3334045552 0.391551655 437.0023656
        -0 0 -0 1  ]};
end
intrMat = eye(3);
intrMat(1,1) = camParam(1);
intrMat(2,2) = camParam(2);
intrMat(1,3) = camParam(3);
intrMat(2,3) = camParam(4);
dist = camParam(5:end)';


for i =  1 : size(camParams,1)
    intrMat_ = eye(3);
    intrMat_(1,1) = camParams(i,1);
    intrMat_(2,2) = camParams(i,2);
    intrMat_(1,3) = camParams(i,3);
    intrMat_(2,3) = camParams(i,4);
    
    intrMat_undist = intrMat_;
    intrMat_undist(1,1) = intrMat_(1,1)*1.0;
    intrMat_undist(2,2) = intrMat_(2,2)*1.0;
    
    dist_ = camParams(i,5:end)';
    camIntr{i,1} = intrMat_;
    camIntr{i,2} = intrMat; %intrMat_undist;
    camDist{i,1} = dist_;
    DirInfo{i,1} = dir(fullfile(inputDir0,dirInfo0(1+i).name,'*.png'));
end

[xMat, yMat] = meshgrid(1:640, 1:480);
pix = [xMat(:) yMat(:)];


% pt3d = unprojectKB8( intrMat(1,1),  intrMat(2,2),  intrMat(1,3),  intrMat(2,3),  dist(1),  dist(2),  dist(3),  dist(4), pix);
% [~,norm_3d] = NormalizeVector(pt3d);
% [pt2d, inImageFlag] = projectKB8(pextend(pt3d')', intrMat, dist, eye(3), [0;0;0]);
% 
% [~,err] = NormalizeVector(pt2d - pix);
% validId = abs(norm_3d-1) < 0.0000001;
% mask = zeros(480, 640);
% mask(validId) = 1;
% figure(30),imshow(mask)
if useRS
    pixUndist_gt = Orig2Rect(pix, intrMat, intrMat, eye(3), dist); % 畸变图整数->无畸变图浮点
else
    pixUndist_gt = remapRect(pix', intrMat, intrMat, dist, eye(3));% 无畸变图整数->畸变图浮点
end
for i = 1 : length(dirInfo)
    
    img = imread(fullfile(inputDir0,dirInfo0(1).name, dirInfo(i).name));
    depth = double(imread(fullfile(inputDepthDir, dirDepthInfo(i).name)))./8;
    %     depth = RemoveDepthHoles(depth0==0, areaRatio);
    depth_mask = depth~=0;
    se = strel('disk',3);
    depth_mask = imerode(depth_mask,se);
    %         depth_mask = imdilate(depth_mask,se);
    depth = immultiply(depth, depth_mask);
        
    img_undist = interp2(xMat,yMat,double(img),reshape(pixUndist_gt(:,1),480,640),reshape(pixUndist_gt(:,2),480,640));
    depth_undist = interp2(xMat,yMat,double(depth),reshape(pixUndist_gt(:,1),480,640),reshape(pixUndist_gt(:,2),480,640));
%     depth_undist = 0.96.*depth_undist;
    depth_undist = 1.0 *depth_undist;
    if 1
        depth_mask = depth_undist~=0;
        se = strel('disk',3);
        depth_mask = imerode(depth_mask,se);
%         depth_mask = imdilate(depth_mask,se);
        depth_undist = immultiply(depth_undist, depth_mask);
        depth_undist = RemoveDepthHoles(depth_undist, areaRatio);
        depth_undist(isnan(depth_undist)) = 0;
    end
    
    pix_raw = [40 134; 155 156; 144 264; 20 250];
    ind = sub2ind(size(depth), pix_raw(:,2), pix_raw(:,2));
    depth_list = depth(ind);
    
    pix_select_undist = Orig2Rect(pix_raw, intrMat, intrMat, eye(3), dist);
    pix_select_undist_all = Orig2Rect(pix, intrMat, intrMat, eye(3), dist);
    
    
    [XYZ] = GetXYZFromDepth(intrMat, pix_select_undist,depth_list);
    [XYZ] = GetXYZFromDepth(intrMat, pix_raw,depth_list);
    
    [XYZ_all] = GetXYZFromDepth(intrMat, pix_select_undist_all,depth(:));
    [XYZ_all_warp] = GetXYZFromDepth(intrMat, pix,depth_undist(:));
    [XYZ_all_orig] = GetXYZFromDepth(intrMat, pix,depth(:));
    
    figure(50),clf;subplot(1,3,1);pcshow(XYZ_all_warp);title('warpped pix');subplot(1,3,2);pcshow(XYZ_all_orig);title('original pix');subplot(1,3,3);pcshowpair(pointCloud(XYZ_all_warp), pointCloud(XYZ_all_orig));
    figure(51),clf;subplot(1,3,1);pcshow(XYZ_all);title('undist pix');subplot(1,3,2);pcshow(XYZ_all_orig);title('original pix');subplot(1,3,3);pcshowpair(pointCloud(XYZ_all), pointCloud(XYZ_all_orig));
%     figure,subplot(1,3,1);pcshow(XYZ_all);title('undist pix');subplot(1,3,2);pcshow(XYZ_all_orig);title('original pix');subplot(1,3,3);pcshowpair(pointCloud(XYZ_all), pointCloud(XYZ_all_orig));
    %    figure,imshow(img_undist,[]);hold on;plot(pix_select_undist(:,1), pix_select_undist(:,2),'.r')
    
    %    norm(XYZ(1,:) - XYZ(2,:))
    %    norm(XYZ(2,:) - XYZ(3,:))
    %    norm(XYZ(3,:) - XYZ(4,:))
    %    norm(XYZ(4,:) - XYZ(1,:))
    
    figure(10),subplot(1,3,1);imshow([double(img) img_undist], []);subplot(1,3,2),imshowpair(depth, img);subplot(1,3,3),imshowpair(depth_undist, img_undist);
   
%     for j = 3:3%size(camParams,1)
    for j = 3:3%size(camParams,1)    
        %%%%% 对鱼眼图做去畸变，计算remap矩阵
        pixDist = remapRectKB8(pix', camIntr{j,2}, camIntr{j,1},camDist{j,1}(1:4), eye(3));
        
        img_fish = imread(fullfile(inputDir0,dirInfo0(1+j).name, DirInfo{j}(i).name));
        img_fish_undist = interp2(xMat,yMat,double(img_fish),reshape(pixDist(:,1),480,640),reshape(pixDist(:,2),480,640));
        pose_fish_depth = camExs{j};
        % 把rgb的内容渲染到depth坐标系
        % 把6目的内容渲染到0坐标系
        % imgCur->6目 / rgb
        % imgPrv->0目 / depth
        % optPose-> depth2rgb / 026 
        % [warped_image1, valid_mask11, warped_z1, warped_z_mask1] = warpImage(double(imgCur), intrMatL(1)*baseline./prv_disp, optPose, intrMatL, intrMatRGB, double(imgPrv));
        
        
        
        
%         [pose_rel, score, warped_image, valid_mask, weightMap, warped_z, validMask] = estimateVisualOdometry(img_curr, img_prev, dep_curr, dep_prev, K, num_levels, pose_init, isVerbose,initial_sigma,default_dof,para);
        [pose_rel, score, warped_image1, valid_mask, weightMap, warped_z, validMask] = estimateVisualOdometry(img_fish_undist, img_undist, depth_undist, depth_undist,intrMat, num_levels, pose_fish_depth, 0,initial_sigma,default_dof,para);
        [warped_image, valid_mask, warped_z, validMask] = warpImage(img_fish_undist, depth_undist, pose_rel, intrMat, camIntr{j,2}, img_undist);
        [warped_image2, valid_mask2, warped_z2, validMask2] = warpImage(img_fish_undist, depth_undist, pose_fish_depth, intrMat, camIntr{j,2}, img_undist);
        
        figure(60),clf;subplot(1,2,1);imshowpair(warped_image, img_undist);title(sprintf('opt, frame: %d, useRS: %d, num-levels: %d', i,useRS, num_levels));subplot(1,2,2);imshowpair(warped_image2, img_undist);title('cam Ex');
        figure(15),clf;imshow([img_fish img_fish_undist img_undist], [])
    end
    
end

end



