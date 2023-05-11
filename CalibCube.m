function CalibCube()

% close all


%���ڱ궨��ͼ��
imgNum = 100; 10; 1; 100;

errStack = [];
stereoErrOld = [];
stereoErrNew = [];

ptTrackL = {};
ptTrackR = {};
%���÷������100
for i = 1 : 500 %100 % 500
    % [pt2dL, pt2dR, pt3dL, intrMatL, intrMatR] = genCalibData();
    %���ɷ�������
    for j = 1 : imgNum
%         [ptTrackL{j,1}, ptTrackR{j,1}, imgSize, intrMatL, intrMatR, Marker2CamL{j,1}, Marker2CamR{j,1}, L2R, cbSize, markerRow, rt2, rt3] = genCalibData(j, ptTrackL, ptTrackR);
        [ptTrackL{j,1}, ptTrackR{j,1}, imgSize, intrMatL, intrMatR, Marker2CamL{j,1}, Marker2CamR{j,1}, L2R, cbSize, markerRow, rt2, rt3, kcL, kcR] = genCalibData();
        checkDiff = inv(L2R*Marker2CamL{j})*Marker2CamR{j};
    end
%      [ptTrackL, ptTrackR, imgSize, intrMatL, intrMatR, Marker2CamL, Marker2CamR, L2R, cbSize, markerRow, rt2, rt3] = genCalibData();
    config.dX = cbSize;
    config.dY = cbSize;
    config.gridSize = markerRow;
    config.rt2 = rt2;
    config.rt3 = rt3;
    config.estDistortion = [1;1;1;1;1];
    
    %�궨��Ŀ�ڲ�
    [camParamL, camParamL1, cbcXYL, cbGridL, cbGridUseL, cbcXYL_, cbGridL_, configL, errCalibL] = CalibCubeSingle(ptTrackL, imgSize, config, Marker2CamL, intrMatL, kcL);
    %�궨��Ŀ�ڲ�
    [camParamR, camParamR1, cbcXYR, cbGridR, cbGridUseR, cbcXYR_, cbGridR_, configR, errCalibR] = CalibCubeSingle(ptTrackR, imgSize, config, Marker2CamR, intrMatR, kcR);
%     [camParamR, camParamR1, cbcXYR, cbGridR, cbGridUseR, cbcXYR_, cbGridR_, configR, errCalibR] = CalibCubeSingle(ptTrackR, imgSize, config, L2R*Marker2CamL, intrMatR);
    
    %�õ����Ż�˫Ŀ�����
    stereoParam1 = CbCalibStereo(camParamL1, camParamR1, cbcXYL_, cbcXYR_, cbGridL_, configL, configR);
    %�������Ż�˫Ŀ�����
    stereoParam = CalibCubeStereo(camParamL, camParamR, cbcXYL, cbcXYR, cbGridUseL, configL, configR);
    
    %���Ƚ�
    stereoErr = [stereoParam.focLeft - [intrMatL(1,1); intrMatL(2,2)];...
                 stereoParam.cenLeft - [intrMatL(1,3); intrMatL(2,3)];...
                 rad2deg(norm(rodrigues(rodrigues(stereoParam.optPose(1:3,1))'*Marker2CamL{1}(1:3,1:3))));...
                 norm(stereoParam.optPose(4:6,1) - Marker2CamL{1}(1:3,4));...
            stereoParam.kcLeft - kcL;...
                 stereoParam.focRight - [intrMatR(1,1); intrMatR(2,2)];...
                 stereoParam.cenRight - [intrMatR(1,3); intrMatR(2,3)];...
                stereoParam.kcRight - kcR;...
                 rad2deg(norm(rodrigues(rodrigues(stereoParam.rotVecRef)'*L2R(1:3,1:3))));...
                 norm(stereoParam.transVecRef - L2R(1:3,4))];


    stereoErr1 = [stereoParam1.focLeft - [intrMatL(1,1); intrMatL(2,2)];...
                  stereoParam1.cenLeft - [intrMatL(1,3); intrMatL(2,3)];...
                  rad2deg(norm(rodrigues(rodrigues(stereoParam1.optPose(1:3,1))'*Marker2CamL{1}(1:3,1:3))));...
                  norm(stereoParam1.optPose(4:6,1) - Marker2CamL{1}(1:3,4));...
            stereoParam1.kcLeft - kcL;...
                 stereoParam1.focRight - [intrMatR(1,1); intrMatR(2,2)];...
                 stereoParam1.cenRight - [intrMatR(1,3); intrMatR(2,3)];...
                 
            stereoParam1.kcRight - kcR;...
                 rad2deg(norm(rodrigues(rodrigues(stereoParam1.rotVecRef)'*L2R(1:3,1:3))));...
                 norm(stereoParam1.transVecRef - L2R(1:3,4))];

    stereoErr12 = [stereoErr stereoErr1 [errCalibL; errCalibR]];
    stereoErrNew = [stereoErrNew stereoErr];
    stereoErrOld = [stereoErrOld stereoErr1];
    fprintf('=====================================================================');
% [camParamR, cbcXYR, cbGridR, configR] = CalibCubeSingle(ptTrackR, imgSize, config, Marker2CamR);
    errStack = [errStack; [errCalibL;  errCalibR]];
end


if 0
    figure,subplot(1,2,1);plot(errStack(:,[1 3]));title('rot diff (deg)');subplot(1,2,2);plot(errStack(:,[2 4]));title('trans diff (mm)');
end

%���ͼʾ����ɫ��ʾ���������ɫ��ʾ����궨��
figure,subplot(2,4,1);hist(errStack(1:11:end,:), 100);title('fx');legend('new','old');
subplot(2,4,5);hist(errStack(2:11:end,:), 100);title('fy');
subplot(2,4,2);hist(errStack(3:11:end,:), 100);title('cx');
subplot(2,4,6);hist(errStack(4:11:end,:), 100);title('cy');
subplot(2,4,3);hist(errStack(5:11:end,:), 100);title('ang');
subplot(2,4,7);hist(errStack(6:11:end,:), 100);title('trans');

figure,subplot(2,2,1);hist(errStack(7:11:end,:), 100);title('k1');legend('new','old');
subplot(2,2,2);hist(errStack(8:11:end,:), 100);title('k2');
subplot(2,2,3);hist(errStack(9:11:end,:), 100);title('p1');
subplot(2,2,4);hist(errStack(10:11:end,:), 100);title('p2');
% subplot(1,5,5);hist(errStack(11:11:end,:), 100);title('k3');



end