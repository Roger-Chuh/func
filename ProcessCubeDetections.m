function ProcessCubeDetections()

% close all;
if 1
    img = imread('C:\vs\test\ConsoleApplication2\ConsoleApplication2\5\img_002.jpg');
    corner = load('C:\vs\test\ConsoleApplication2\ConsoleApplication2\5\corner_00001.txt');
elseif 0
    
    img = imread('C:\vs\test\ConsoleApplication2\ConsoleApplication2\11\marker0.png');
    corner = load('C:\vs\test\ConsoleApplication2\ConsoleApplication2\11\corner_00001.txt');
elseif 0
    img = imread('C:\vs\test\ConsoleApplication2\ConsoleApplication2\12\marker0.png');
    corner = load('C:\vs\test\ConsoleApplication2\ConsoleApplication2\12\corner_00001.txt');
    
    
elseif 0
    img = imread('C:\vs\test\ConsoleApplication2\ConsoleApplication2\7\img_001.jpg');
    corner = load('C:\vs\test\ConsoleApplication2\ConsoleApplication2\7\corner_00001.txt');
elseif 0
    img = imread('C:\vs\test\ConsoleApplication2\ConsoleApplication2\8\img_001.png');
    corner = load('C:\vs\test\ConsoleApplication2\ConsoleApplication2\8\corner_00001.txt');
elseif 0
    img = imresize(imread('C:\vs\test\ConsoleApplication2\ConsoleApplication2\9\imgL_0001.png'),3);
    corner = load('C:\vs\test\ConsoleApplication2\ConsoleApplication2\9\corner_00001.txt');
elseif 0
    img = imresize(imread('C:\vs\test\ConsoleApplication2\ConsoleApplication2\10\imgL_0001.png'),3);
    corner = load('C:\vs\test\ConsoleApplication2\ConsoleApplication2\10\corner_00001.txt');
elseif 0
    img = imresize(imread('C:\vs\test\ConsoleApplication4\data\10\imgL_0001.png'),3);
    corner = load('C:\vs\test\ConsoleApplication4\data\10\cornersL_00001.txt');
else
    img = imread('C:\vs\test\ConsoleApplication2\ConsoleApplication2\13\marker.jpg');
    corner = load('C:\vs\test\ConsoleApplication2\ConsoleApplication2\13\corner_00001.txt');
end


rotVecRef = [0.00343335720939389;0.00358147150176189;0.00395848519477192];
transVecRef = [-25.2881486549827;0.119516509409703;0.0881998428058557];

L2R = [rodrigues(rotVecRef) transVecRef; 0 0 0 1];
Cam2MarkerRot = rodrigues([-0.6 0.8 0.2] + 1.0*(rand(1,3)-0.5));
Cam2MarkerTrans = [-1000 -800 -900]';
Cam2Marker = [Cam2MarkerRot Cam2MarkerTrans; 0 0 0 1];
Marker2CamL = inv(Cam2Marker);
Marker2CamR = L2R*Marker2CamL;




% bigFace = 1;
bigFace = 0;
midFace = 1;0;

corner(:,2:end) = corner(:,2:end) + 1;
id1 = find(corner(:,1) < 1000);
id2 = find(corner(:,1) > 1000);

idOffset = 10000;
coordOffset = 1.234; 50;

wrongOrder = 1;
wrongRt = 0;
wringId = 0;


figure,imshow(img);hold on

for i = 1 : length(id1)
    text(corner(id1(i),2)+3,corner(id1(i),3)+3,num2str(corner(id1(i),1)), 'Color',[1 0 0],'FontSize',15);
    plot(corner(id1(i),2), corner(id1(i),3),'or')
    
end

for i = 1 : length(id2)
    text(corner(id2(i),2)+3,corner(id2(i),3)+3,num2str(corner(id2(i),1)-idOffset), 'Color',[0 0 1],'FontSize',15);
    plot(corner(id2(i),2), corner(id2(i),3),'ob');
    
end

if ~bigFace
    if ~midFace
        faceNum = 8;
        rowIdStep  =28;
        startBottomRowId = [0 10 20];
    else
        faceNum = 13;
        rowIdStep  = 43;
        startBottomRowId = [0 15 30];
    end
else
    faceNum = 19;
    rowIdStep  = 61;
    startBottomRowId = [0 21 42];
end
markerSize = 50;
markerCol = faceNum;
markerRow = faceNum;

redFaceMat = genFaceData(startBottomRowId(1), rowIdStep,faceNum);
greenFaceMat = genFaceData(startBottomRowId(2), rowIdStep,faceNum);
blueFaceMat = genFaceData(startBottomRowId(3), rowIdStep,faceNum);


[xMat0, yMat0] = meshgrid(markerSize : markerSize : (markerCol - 0) * markerSize, markerSize : markerSize : (markerRow - 0) * markerSize);
xMat = xMat0 + coordOffset;
yMat = yMat0 + coordOffset;
marker1 = [xMat(:) yMat(:) zeros(markerRow*markerCol,1)];



if ~wrongOrder
    marker2 = flipud((rotx(-90)*rotz(-90)*marker1')');
    marker3 = flipud((roty(90)*rotz(90)*marker1')');
else
    marker3 = flipud((rotx(-90)*rotz(-90)*marker1')');
    marker2 = flipud((roty(90)*rotz(90)*marker1')');
end
if ~wrongRt
    rtOut21 = AlignCoord(marker2,marker1);
    rtOut31 = AlignCoord(marker3,marker1);
else
    rtOut31 = AlignCoord(marker2,marker1);
    rtOut21 = AlignCoord(marker3,marker1);
end
if 0
    % bug solved
    markers1 = [marker1; marker2; marker3];
else
    markers1 = [marker3; marker2; marker1];
%     markers1 = [marker2; marker3; marker1];
end
markers11 = [marker1;marker1;marker1];


% 看成列先就是凸的

figure,pcshow(marker1,[1 0 0]);hold on;pcshow(marker2,[0 1 0]);pcshow(marker3,[0 0 1]);
plot3(marker1(:,1),marker1(:,2),marker1(:,3),'-or');plot3(marker1(1,1),marker1(1,2),marker1(1,3),'*r','MarkerSize',10,'LineWidth', 10);
plot3(marker2(:,1),marker2(:,2),marker2(:,3),'-og');plot3(marker2(1,1),marker2(1,2),marker2(1,3),'*g','MarkerSize',10,'LineWidth', 10);
plot3(marker3(:,1),marker3(:,2),marker3(:,3),'-ob');plot3(marker3(1,1),marker3(1,2),marker3(1,3),'*b','MarkerSize',10,'LineWidth', 10);
if ~wringId
    faceIdList = 10000 + [blueFaceMat(:); greenFaceMat(:); redFaceMat(:)];
else
    %      faceIdList = 10000 + [greenFaceMat(:); blueFaceMat(:);  redFaceMat(:)];
    faceIdList = 10000 + [blueFaceMat(:); redFaceMat(:); greenFaceMat(:)];
end


% 看成行先就是凹的

marker11 = [marker1(:,1) marker1(:,2) marker1(:,3)];
if ~wrongOrder
    marker22 = flipud((rotx(90)*rotz(-90)*marker11')');
    marker33 = flipud((roty(-90)*rotz(90)*marker11')');
else
    marker33 = flipud((rotx(90)*rotz(-90)*marker11')');
    marker22 = flipud((roty(-90)*rotz(90)*marker11')');
end
if ~wrongRt
    rtIn21 = AlignCoord(marker22,marker11);
    rtIn31 = AlignCoord(marker33,marker11);
else
    rtIn31 = AlignCoord(marker22,marker11);
    rtIn21 = AlignCoord(marker33,marker11);
end
figure,pcshow(marker11,[1 0 0]);hold on;pcshow(marker22,[0 1 0]);pcshow(marker33,[0 0 1]);
plot3(marker11(:,1),marker11(:,2),marker11(:,3),'-or');plot3(marker11(1,1),marker11(1,2),marker11(1,3),'*r','MarkerSize',10,'LineWidth', 10);
plot3(marker22(:,1),marker22(:,2),marker22(:,3),'-og');plot3(marker22(1,1),marker22(1,2),marker22(1,3),'*g','MarkerSize',10,'LineWidth', 10);
plot3(marker33(:,1),marker33(:,2),marker33(:,3),'-ob');plot3(marker33(1,1),marker33(1,2),marker33(1,3),'*b','MarkerSize',10,'LineWidth', 10);
if 0
    % bug solved
    markers2 = [marker11; marker22; marker33];
else
    if 0
        markers2 = [ marker33; marker22; marker11];
    else
        % don't know why i changed to this, it just solved the bug
        % 20211208
        markers2 = [ marker22; marker33; marker11];
    end
end
markers22 = [marker11; marker11; marker11];

a = load('D:\facTest3\testCube\06\data.txt');
img = imread('D:\facTest3\testCube\06\imgL_0001.jpg');


inputDir = pwd;

intrMatL = [400 0 320;0 400 240; 0 0 1];
intrMatR = intrMatL;
[pt2dL0, pt3dL] = TransformAndProject(markers2, intrMatL, Marker2CamL(1:3,1:3), Marker2CamL(1:3,4));
[pt2dR0, pt3dR] = TransformAndProject(markers2, intrMatR, Marker2CamR(1:3,1:3), Marker2CamR(1:3,4));



if ~bigFace
    if ~midFace
        fid1 = fopen(fullfile(inputDir,'gridSmall.yaml'),'w');
    else
        fid1 = fopen(fullfile(inputDir,'gridMid.yaml'),'w');
    end
else
    fid1 = fopen(fullfile(inputDir,'gridBig.yaml'),'w');
end

% -----------------------------
fprintf(fid1,'%%YAML:1.0\n');
fprintf(fid1,'\n');

% -----------------------------
fprintf(fid1,sprintf('faceIdList: !!opencv-matrix\n   rows: %d\n   cols: %d\n   dt: d\n', length(faceIdList), 1));
fprintf(fid1,sprintf('   data: [ %d, \n',faceIdList(1)));
for i = 2 : length(faceIdList)-1
    fprintf(fid1,sprintf('           %d, \n',faceIdList(i)));
end
fprintf(fid1,sprintf('           %d] \n',faceIdList(end)));

% ------------------------------

fprintf(fid1,sprintf('outsideGrid: !!opencv-matrix\n   rows: %d\n   cols: %d\n   dt: d\n',size(markers1,1),size(markers1,2)));
fprintf(fid1,sprintf('   data: [ %0.15f, %0.15f,  %0.15f, \n',markers1(1,1), markers1(1,2),markers1(1,3)));
for i = 2 : size(markers1,1)-1
    fprintf(fid1,sprintf('           %0.15f, %0.15f, %0.15f, \n',markers1(i,1), markers1(i,2), markers1(i,3)));
end
fprintf(fid1,sprintf('           %0.15f, %0.15f, %0.15f] \n',markers1(end,1), markers1(end,2), markers1(end,3)));


fprintf(fid1,sprintf('outsideGridSingle: !!opencv-matrix\n   rows: %d\n   cols: %d\n   dt: d\n',size(markers11,1),size(markers11,2)));
fprintf(fid1,sprintf('   data: [ %0.15f, %0.15f,  %0.15f, \n',markers11(1,1), markers11(1,2),markers11(1,3)));
for i = 2 : size(markers11,1)-1
    fprintf(fid1,sprintf('           %0.15f, %0.15f, %0.15f, \n',markers11(i,1), markers11(i,2), markers11(i,3)));
end
fprintf(fid1,sprintf('           %0.15f, %0.15f, %0.15f] \n',markers11(end,1), markers11(end,2), markers11(end,3)));


fprintf(fid1,sprintf('rtOut21: !!opencv-matrix\n   rows: %d\n   cols: %d\n   dt: d\n', 4, 4));
fprintf(fid1,sprintf('   data: [ %0.15f, %0.15f,  %0.15f, %0.15f, \n',rtOut21(1,1), rtOut21(1,2),rtOut21(1,3), rtOut21(1,4)));
for i = 2 : 3
    fprintf(fid1,sprintf('           %0.15f, %0.15f, %0.15f, %0.15f, \n',rtOut21(i,1), rtOut21(i,2), rtOut21(i,3),rtOut21(i,4)));
end
fprintf(fid1,sprintf('           %0.15f, %0.15f, %0.15f, %0.15f] \n',rtOut21(end,1), rtOut21(end,2), rtOut21(end,3), rtOut21(end,4)));

fprintf(fid1,sprintf('rtOut31: !!opencv-matrix\n   rows: %d\n   cols: %d\n   dt: d\n', 4, 4));
fprintf(fid1,sprintf('   data: [ %0.15f, %0.15f,  %0.15f, %0.15f, \n',rtOut31(1,1), rtOut31(1,2),rtOut31(1,3), rtOut31(1,4)));
for i = 2 : 3
    fprintf(fid1,sprintf('           %0.15f, %0.15f, %0.15f, %0.15f, \n',rtOut31(i,1), rtOut31(i,2), rtOut31(i,3),rtOut31(i,4)));
end
fprintf(fid1,sprintf('           %0.15f, %0.15f, %0.15f, %0.15f] \n',rtOut31(end,1), rtOut31(end,2), rtOut31(end,3), rtOut31(end,4)));

% ---------------------------------

fprintf(fid1,sprintf('insideGrid: !!opencv-matrix\n   rows: %d\n   cols: %d\n   dt: d\n',size(markers2,1),size(markers2,2)));
fprintf(fid1,sprintf('   data: [ %0.15f, %0.15f,  %0.15f, \n',markers2(1,1), markers2(1,2),markers2(1,3)));
for i = 2 : size(markers2,1)-1
    fprintf(fid1,sprintf('           %0.15f, %0.15f, %0.15f, \n',markers2(i,1), markers2(i,2), markers2(i,3)));
end
fprintf(fid1,sprintf('           %0.15f, %0.15f, %0.15f] \n',markers2(end,1), markers2(end,2), markers2(end,3)));


fprintf(fid1,sprintf('insideGridSingle: !!opencv-matrix\n   rows: %d\n   cols: %d\n   dt: d\n',size(markers22,1),size(markers22,2)));
fprintf(fid1,sprintf('   data: [ %0.15f, %0.15f,  %0.15f, \n',markers22(1,1), markers22(1,2),markers22(1,3)));
for i = 2 : size(markers22,1)-1
    fprintf(fid1,sprintf('           %0.15f, %0.15f, %0.15f, \n',markers22(i,1), markers22(i,2), markers22(i,3)));
end
fprintf(fid1,sprintf('           %0.15f, %0.15f, %0.15f] \n',markers22(end,1), markers22(end,2), markers22(end,3)));


fprintf(fid1,sprintf('rtIn21: !!opencv-matrix\n   rows: %d\n   cols: %d\n   dt: d\n', 4, 4));
fprintf(fid1,sprintf('   data: [ %0.15f, %0.15f,  %0.15f, %0.15f, \n',rtIn21(1,1), rtIn21(1,2),rtIn21(1,3), rtIn21(1,4)));
for i = 2 : 3
    fprintf(fid1,sprintf('           %0.15f, %0.15f, %0.15f, %0.15f, \n',rtIn21(i,1), rtIn21(i,2), rtIn21(i,3),rtIn21(i,4)));
end
fprintf(fid1,sprintf('           %0.15f, %0.15f, %0.15f, %0.15f] \n',rtIn21(end,1), rtIn21(end,2), rtIn21(end,3), rtIn21(end,4)));

fprintf(fid1,sprintf('rtIn31: !!opencv-matrix\n   rows: %d\n   cols: %d\n   dt: d\n', 4, 4));
fprintf(fid1,sprintf('   data: [ %0.15f, %0.15f,  %0.15f, %0.15f, \n',rtIn31(1,1), rtIn31(1,2),rtIn31(1,3), rtIn31(1,4)));
for i = 2 : 3
    fprintf(fid1,sprintf('           %0.15f, %0.15f, %0.15f, %0.15f, \n',rtIn31(i,1), rtIn31(i,2), rtIn31(i,3),rtIn31(i,4)));
end
fprintf(fid1,sprintf('           %0.15f, %0.15f, %0.15f, %0.15f] \n',rtIn31(end,1), rtIn31(end,2), rtIn31(end,3), rtIn31(end,4)));



fclose(fid1);









end

function FaceIdMat = genFaceData(startBottomRowId, rowIdStep,faceNum)

FaceId0 = startBottomRowId:rowIdStep:(100)*rowIdStep;
FaceId = FaceId0(1:faceNum);
idOfst = 0:faceNum-1;

FaceIds = FaceId + idOfst';

FaceIds = sort(FaceIds(:),'descend');
FaceIdMat = fliplr(reshape(FaceIds, faceNum, faceNum)');

end