% function [pt2dL, pt2dR, pt3dL, intrMatL, intrMatR] = genCalibData()
% function [ptTrackL, ptTrackR, imgSize, intrMatL, intrMatR, Marker2CamL, Marker2CamR, L2R, markerSize, markerRow, rt2, rt3] = genCalibData(ind,ptTrackL, ptTrackR)
function [ptTrackL, ptTrackR, imgSize, intrMatL, intrMatR, Marker2CamL, Marker2CamR, L2R, markerSize, markerRow, rt2, rt3, kcL, kcR] = genCalibData()

% close all

use_inverse_BC = true;

sigma_pixel = 0; 2; 0; 2; 0; 2; 1; 0; 2; 1; 0;

coord_offset = 0; -600;



flip = 1; 0; 1; 0; 1;

imgWidth = 640;
imgHeight = 480;
imgSize = [imgHeight imgWidth];
if 1
    intrMatL = [501 0 321; 0 499 241; 0 0 1];
    intrMatR = [502 0 319; 0 498 239; 0 0 1];
    camParam = [475.399, 475.399 318.445, 245.336 [0.153221,0.101447,0.00371338,0.00427771,0.002554]];
    
    if use_inverse_BC
        intrMatL = [camParam(1,1) 0 camParam(1,3); 0 camParam(1,2) camParam(1,4); 0 0 1];
        intrMatR = intrMatL;
    end
    
else
    intrMatL = [501 0 319.5; 0 501 239.5; 0 0 1];
    intrMatR = intrMatL; %[502 0 319; 0 501 239; 0 0 1];
end
kcL = [-0.0762618615083300;0.0897966202834870;-0.00270054859097500;-0.00304068736751200;0];
kcR = [-0.0290863163581273;0.0164066955445654;-0.00293902113085841;-0.00354834432702059;0];

kcL = -[-0.5762618615083300;-0.2897966202834870;-0.00270054859097500;-0.00304068736751200;-0];
kcR = -[-0.5290863163581273;-0.2164066955445654;-0.00293902113085841;-0.00354834432702059;-0.];

if ~use_inverse_BC
    kcL = -[-0.1762618615083300;-0.0897966202834870;-0.00270054859097500;-0.00304068736751200;-0];
    kcR = -[-0.190863163581273;-0.064066955445654;-0.00293902113085841;-0.00354834432702059;-0.];
else
    kcL = camParam(5:end)';
    kcR = camParam(5:end)';
end

rotVecRef = [0.00343335720939389;0.00358147150176189;0.00395848519477192];
transVecRef = [-25.2881486549827;0.119516509409703;0.0881998428058557];

L2R = [rodrigues(rotVecRef) transVecRef; 0 0 0 1];
if 0
    Cam2MarkerRot = rodrigues([-0.6 0.8 0.2]);
    Cam2MarkerTrans = [-1000 -800 -900]';
else
    if 0
        Cam2MarkerRot = rodrigues([-0.6 -0.3 0.2]);
    else
        if 0
            Cam2MarkerRot = rodrigues([-0.6 0.8 0.2] + 1.0*(rand(1,3)-0.5));
        else
            Cam2MarkerRot = rodrigues([-0.6 0.8 0.2] + 0.5*(rand(1,3)-0.5));
        end
    end
    Cam2MarkerTrans = [-1000 -800 -900]';
    
end
Cam2Marker = [Cam2MarkerRot Cam2MarkerTrans; 0 0 0 1];
Marker2CamL = inv(Cam2Marker);
Marker2CamR = L2R*Marker2CamL;

markerSize = 150; 100; 200;
markerRow = 15; 4; 7; 15;20; 7;
markerCol = 15; 4; 7; 15;20; 7;

[xMat0, yMat0] = meshgrid(markerSize : markerSize : (markerCol - 0) * markerSize, markerSize : markerSize : (markerRow - 0) * markerSize);
xMat = xMat0 + coord_offset;
yMat = yMat0 + coord_offset;

marker1 = [xMat(:) yMat(:) zeros(markerRow*markerCol,1)];
if 1
    if 0
        marker2_ = (rotx(90)*marker1')';
        marker2x = reshape(marker2_(:,1), markerRow, markerRow)';
        marker2y = reshape(marker2_(:,2), markerRow, markerRow)';
        marker2z = reshape(marker2_(:,3), markerRow, markerRow)';
        marker2 = [marker2x(:) marker2y(:) marker2z(:)];
        [tform,~,rmse] = pcregrigid(pointCloud(marker1), pointCloud(marker2),'InlierRatio',1,'MaxIterations',1000,'Extrapolate',true); %,'Metric','pointToPlane');
    else
        if flip
            marker2 = flipud((rotx(-90)*rotz(-90)*marker1')');
            
        else
            marker2 = ((rotx(-90)*rotz(-90)*marker1')');
        end
        rt2 = AlignCoord(marker2,marker1);
    end
    
    
    if 0
        marker3_ = (roty(-90)*marker1')';
        marker3x = reshape(marker3_(:,1), markerRow, markerRow)';
        marker3y = reshape(marker3_(:,2), markerRow, markerRow)';
        marker3z = reshape(marker3_(:,3), markerRow, markerRow)';
        marker3 = [marker3x(:) marker3y(:) marker3z(:)];
    elseif 0
        marker3 = (rotz(90)*marker2')';
    else
        if flip
            marker3 = flipud((roty(90)*rotz(90)*marker1')');
            
        else
            marker3 = ((roty(90)*rotz(90)*marker1')');
        end
        rt3 = AlignCoord(marker3,marker1);
    end
else
    
    marker2 = (rotx(-90)*marker1')';
    if 1
        marker3 = (roty(90)*marker1')';
    else
        marker3 = (rotz(90)*marker2')';
    end
end
marker = [marker1; marker2; marker3];
[pt2dL00, pt3dL] = TransformAndProject(marker, intrMatL, Marker2CamL(1:3,1:3), Marker2CamL(1:3,4));
if ~use_inverse_BC
    pt2dL0 = remapRect(pt2dL00', intrMatL, intrMatL, kcL, eye(3)); % gennerate distorted pixel
else
    pt2dL0 = Orig2Rect(pt2dL00, intrMatL, intrMatL,  eye(3), kcL); % gennerate distorted pixel
end
[pt2dR00, pt3dR] = TransformAndProject(marker, intrMatR, Marker2CamR(1:3,1:3), Marker2CamR(1:3,4));
if ~use_inverse_BC
    pt2dR0 = remapRect(pt2dR00', intrMatR, intrMatR, kcR,eye(3)); % gennerate distorted pixel
else
    pt2dR0 = Orig2Rect(pt2dR00, intrMatR, intrMatR,eye(3), kcR); % gennerate distorted pixel
end
err2 = [rotx(-90)*rotz(-90) [0 0 0]';0 0 0 1]*inv(Marker2CamL)*pextend(pt3dL(1:markerRow*markerCol,:)') - pextend(flipud(marker2)');
err3 = [roty(90)*rotz(90) [0 0 0]';0 0 0 1]*inv(Marker2CamL)*pextend(pt3dL(1:markerRow*markerCol,:)') - pextend(flipud(marker3)');

err2_ = [rotx(-90)*rotz(-90) [0 0 0]';0 0 0 1]*inv(Marker2CamL)*pextend(pt3dL(1:markerRow*markerCol,:)') - pextend((marker2)');
err3_ = [roty(90)*rotz(90) [0 0 0]';0 0 0 1]*inv(Marker2CamL)*pextend(pt3dL(1:markerRow*markerCol,:)') - pextend((marker3)');


aa = (L2R*pextend(pt3dL'))';
err = aa(:,1:3) - pt3dR;

err_lx = normrnd(0, sigma_pixel, size(pt2dL0,1),1);
err_ly = normrnd(0, sigma_pixel, size(pt2dL0,1),1);
err_rx = normrnd(0, sigma_pixel, size(pt2dR0,1),1);
err_ry = normrnd(0, sigma_pixel, size(pt2dR0,1),1);


pt2dL = pt2dL0 + [err_lx err_ly];
pt2dR = pt2dR0 + [err_rx err_ry];

ind1 = [1:size(marker1,1)]';
ind2 = [size(marker1,1)+1:2*size(marker1,1)]';
ind3 = [2*size(marker1,1)+1:3*size(marker1,1)]';



ptTrackL{1,1} = [{ind1} {[pt2dL(ind1,:) marker(ind1,:)]}; {ind2} {[pt2dL(ind2,:) marker(ind2,:)]}; {ind3} {[pt2dL(ind3,:) marker(ind3,:)]}];
ptTrackR{1,1} = [{ind1} {[pt2dR(ind1,:) marker(ind1,:)]}; {ind2} {[pt2dR(ind2,:) marker(ind2,:)]}; {ind3} {[pt2dR(ind3,:) marker(ind3,:)]}];


% figure,subplot(2,2,1);pcshow(pix1, [1 0 0]);hold on;pcshow(pix2, [0 1 0]);pcshow(pix3, [0 0 1]);
% subplot(2,2,2);pcshow(pt3d(ind1,:), [1 0 0]);hold on;pcshow(pt3d(ind2,:), [0 1 0]);pcshow(pt3d(ind3,:), [0 0 1]);

if 0
    figure(1),clf;
    subplot(2,2,1);plotPath(marker(ind1,:), [1 0 0]);plotPath(marker(ind2,:), [0 1 0]);plotPath(marker(ind3,:), [0 0 1]);
    subplot(2,2,2);plotPath(pt3dL(ind1,:), [1 0 0]);plotPath(pt3dL(ind2,:), [0 1 0]);plotPath(pt3dL(ind3,:), [0 0 1]);
    
    subplot(2,2,3);
    imshow(zeros(imgHeight, imgWidth)); hold on;
    plot(pt2dL(ind1,1), pt2dL(ind1,2),'-xr');plot(pt2dL(ind1(1),1), pt2dL(ind1(1),2),'or','MarkerSize',5, 'LineWidth',5);
    plot(pt2dL(ind2,1), pt2dL(ind2,2),'-xg');plot(pt2dL(ind2(1),1), pt2dL(ind2(1),2),'og','MarkerSize',5, 'LineWidth',5);
    plot(pt2dL(ind3,1), pt2dL(ind3,2),'-xb');plot(pt2dL(ind3(1),1), pt2dL(ind3(1),2),'ob','MarkerSize',5, 'LineWidth',5);
    title('L');
    
    subplot(2,2,4);
    imshow(zeros(imgHeight, imgWidth)); hold on;
    plot(pt2dR(ind1,1), pt2dR(ind1,2),'-xr');plot(pt2dR(ind1(1),1), pt2dR(ind1(1),2),'or','MarkerSize',5, 'LineWidth',5);
    plot(pt2dR(ind2,1), pt2dR(ind2,2),'-xg');plot(pt2dR(ind2(1),1), pt2dR(ind2(1),2),'og','MarkerSize',5, 'LineWidth',5);
    plot(pt2dR(ind3,1), pt2dR(ind3,2),'-xb');plot(pt2dR(ind3(1),1), pt2dR(ind3(1),2),'ob','MarkerSize',5, 'LineWidth',5);
    title('R');
    drawnow;
end


end