function Mat2Yaml(inputDir, scale)
global cfg

img_width = cfg.img_width;
img_height = cfg.img_height;

load(fullfile(inputDir,'calib.mat'));


if scale >= 1

fid1 = fopen(fullfile(inputDir,'extrinsics.yml'),'w');

intrMatL = [stereoParam.focLeft(1) 0 stereoParam.cenLeft(1);...
    0 stereoParam.focLeft(2) stereoParam.cenLeft(2); 0 0 1];
intrMatR = [stereoParam.focRight(1) 0 stereoParam.cenRight(1);...
    0 stereoParam.focRight(2) stereoParam.cenRight(2); 0 0 1];
kcL = stereoParam.kcLeft;
kcR = stereoParam.kcRight;


rotMat = rodrigues(stereoParam.rotVecRef);
transVec = stereoParam.transVecRef;

fprintf(fid1,'%%YAML:1.0\n');
fprintf(fid1,'---\n');
fprintf(fid1,'R: !!opencv-matrix\n   rows: 3\n   cols: 3\n   dt: d\n');
fprintf(fid1,sprintf('   data: [ %0.15f, %0.15f, \n',rotMat(1,1), rotMat(1,2)));
fprintf(fid1,sprintf('       %0.15f, %0.15f, \n',rotMat(1,3), rotMat(2,1)));
fprintf(fid1,sprintf('       %0.15f, %0.15f, \n',rotMat(2,2), rotMat(2,3)));
fprintf(fid1,sprintf('       %0.15f, %0.15f, \n',rotMat(3,1), rotMat(3,2)));
fprintf(fid1,sprintf('       %0.15f ]\n',rotMat(3,3)));

fprintf(fid1,'T: !!opencv-matrix\n   rows: 3\n   cols: 1\n   dt: d\n');
fprintf(fid1,sprintf('   data: [ %0.15f, %0.15f, \n',transVec(1), transVec(2)));
fprintf(fid1,sprintf('       %0.15f ]\n',transVec(3)));


fclose(fid1);


fid1 = fopen(fullfile(inputDir,'intrinsics.yml'),'w');
fprintf(fid1,'%%YAML:1.0\n');
fprintf(fid1,'---\n');
fprintf(fid1,'M1: !!opencv-matrix\n   rows: 3\n   cols: 3\n   dt: d\n');
fprintf(fid1,sprintf('   data: [ %0.15f, %0.15f,  %0.15f, %0.15f, \n',intrMatL(1,1), intrMatL(1,2), intrMatL(1,3), intrMatL(2,1)));
fprintf(fid1,sprintf('       %0.15f, %0.15f, %0.15f, %0.15f, %0.15f ]\n',intrMatL(2,2), intrMatL(2,3),intrMatL(3,1), intrMatL(3,2), intrMatL(3,3)));
fprintf(fid1,'D1: !!opencv-matrix\n   rows: 1\n   cols: 5\n   dt: d\n');
fprintf(fid1,sprintf('   data: [ %0.15f, %0.15f, \n',kcL(1), kcL(2)));
fprintf(fid1,sprintf('       %0.15f, %0.15f, %0.15f ]\n',kcL(3), kcL(4), kcL(5)));

fprintf(fid1,'M2: !!opencv-matrix\n   rows: 3\n   cols: 3\n   dt: d\n');
fprintf(fid1,sprintf('   data: [ %0.15f, %0.15f,  %0.15f, %0.15f, \n',intrMatR(1,1), intrMatR(1,2), intrMatR(1,3), intrMatR(2,1)));
fprintf(fid1,sprintf('       %0.15f, %0.15f, %0.15f, %0.15f, %0.15f ]\n',intrMatR(2,2), intrMatR(2,3),intrMatR(3,1), intrMatR(3,2), intrMatR(3,3)));
fprintf(fid1,'D2: !!opencv-matrix\n   rows: 1\n   cols: 5\n   dt: d\n');
fprintf(fid1,sprintf('   data: [ %0.15f, %0.15f, \n',kcR(1), kcR(2)));
fprintf(fid1,sprintf('       %0.15f, %0.15f, %0.15f ]\n',kcR(3), kcR(4), kcR(5)));

fclose(fid1);

end


if scale == 2
    
    fid1 = fopen(fullfile(inputDir,'extrinsics2.yml'),'w');
    
    intrMatL2 = [stereoParam.focLeft(1)/2 0 stereoParam.cenLeft(1)/2;...
        0 stereoParam.focLeft(2)/2 stereoParam.cenLeft(2)/2; 0 0 1];
    intrMatR2 = [stereoParam.focRight(1)/2 0 stereoParam.cenRight(1)/2;...
        0 stereoParam.focRight(2)/2 stereoParam.cenRight(2)/2; 0 0 1];
    kcL = stereoParam.kcLeft;
    kcR = stereoParam.kcRight;
    
    
    rotMat = rodrigues(stereoParam.rotVecRef);
    transVec = stereoParam.transVecRef;
    
    fprintf(fid1,'%%YAML:1.0\n');
    fprintf(fid1,'---\n');
    fprintf(fid1,'R: !!opencv-matrix\n   rows: 3\n   cols: 3\n   dt: d\n');
    fprintf(fid1,sprintf('   data: [ %0.15f, %0.15f, \n',rotMat(1,1), rotMat(1,2)));
    fprintf(fid1,sprintf('       %0.15f, %0.15f, \n',rotMat(1,3), rotMat(2,1)));
    fprintf(fid1,sprintf('       %0.15f, %0.15f, \n',rotMat(2,2), rotMat(2,3)));
    fprintf(fid1,sprintf('       %0.15f, %0.15f, \n',rotMat(3,1), rotMat(3,2)));
    fprintf(fid1,sprintf('       %0.15f ]\n',rotMat(3,3)));
    
    fprintf(fid1,'T: !!opencv-matrix\n   rows: 3\n   cols: 1\n   dt: d\n');
    fprintf(fid1,sprintf('   data: [ %0.15f, %0.15f, \n',transVec(1), transVec(2)));
    fprintf(fid1,sprintf('       %0.15f ]\n',transVec(3)));
    
    
    fclose(fid1);
    
    
    fid1 = fopen(fullfile(inputDir,'intrinsics2.yml'),'w');
    fprintf(fid1,'%%YAML:1.0\n');
    fprintf(fid1,'---\n');
    fprintf(fid1,'M1: !!opencv-matrix\n   rows: 3\n   cols: 3\n   dt: d\n');
    fprintf(fid1,sprintf('   data: [ %0.15f, %0.15f,  %0.15f, %0.15f, \n',intrMatL2(1,1), intrMatL2(1,2), intrMatL2(1,3), intrMatL2(2,1)));
    fprintf(fid1,sprintf('       %0.15f, %0.15f, %0.15f, %0.15f, %0.15f ]\n',intrMatL2(2,2), intrMatL2(2,3),intrMatL2(3,1), intrMatL2(3,2), intrMatL2(3,3)));
    fprintf(fid1,'D1: !!opencv-matrix\n   rows: 1\n   cols: 5\n   dt: d\n');
    fprintf(fid1,sprintf('   data: [ %0.15f, %0.15f, \n',kcL(1), kcL(2)));
    fprintf(fid1,sprintf('       %0.15f, %0.15f, %0.15f ]\n',kcL(3), kcL(4), kcL(5)));
    
    fprintf(fid1,'M2: !!opencv-matrix\n   rows: 3\n   cols: 3\n   dt: d\n');
    fprintf(fid1,sprintf('   data: [ %0.15f, %0.15f,  %0.15f, %0.15f, \n',intrMatR2(1,1), intrMatR2(1,2), intrMatR2(1,3), intrMatR2(2,1)));
    fprintf(fid1,sprintf('       %0.15f, %0.15f, %0.15f, %0.15f, %0.15f ]\n',intrMatR2(2,2), intrMatR2(2,3),intrMatR2(3,1), intrMatR2(3,2), intrMatR2(3,3)));
    fprintf(fid1,'D2: !!opencv-matrix\n   rows: 1\n   cols: 5\n   dt: d\n');
    fprintf(fid1,sprintf('   data: [ %0.15f, %0.15f, \n',kcR(1), kcR(2)));
    fprintf(fid1,sprintf('       %0.15f, %0.15f, %0.15f ]\n',kcR(3), kcR(4), kcR(5)));
    
    fclose(fid1);
    
    
    
    
end

end