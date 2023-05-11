function [KL, KR, rectParamL, rectParamR, rotMatLeft, rotMatRight,  intrMatLeftNew, intrMatRightNew, meanError, KK_newL] = CheckStereoCalibReprojection(stereoParam, imgSize, imgListL, imgListR, camParamL, camParamR, cbcXYL, cbcXYR, cbGridL, cbGridR, markerIdList, xyzCellGlobal)
global cfg


goodId = [1:size(imgListL,1)]';

nr = imgSize(1);
nc = imgSize(2);



cfg.Corner = {};

KL = [stereoParam.focLeft(1) 0 stereoParam.cenLeft(1);0 stereoParam.focLeft(2) stereoParam.cenLeft(2);0 0 1];
KR = [stereoParam.focRight(1) 0 stereoParam.cenRight(1);0 stereoParam.focRight(2) stereoParam.cenRight(2);0 0 1];


cbcXYLL = cbcXYL;
cbcXYRR = cbcXYR;





% �궨������ϵ�µ�3d��
pt3d = cell2mat(xyzCellGlobal);



[rectParamL, rectParamR, rotMatLeft, rotMatRight,  intrMatLeftNew, intrMatRightNew] = GetRectifyParam2(stereoParam, [nr nc]);
KRectL = rectParamL.rectIntrMat;







KK_newL  = intrMatLeftNew;

princpPtL = intrMatLeftNew(1:2,3);
princpPtR = intrMatRightNew(1:2,3);
cfg.Corner = cell(length(goodId), 3);
cfg.LRIntrMat = intrMatLeftNew;
cfg.LRBaseline = norm(stereoParam.transVecRef);
try
    for i = 1 :  length(goodId) %length(cbcXYL)
        
        % �Խǵ���undistortion
        [xr] = normalize_pixel(cbcXYRR{(i)},stereoParam.focRight,stereoParam.cenRight,stereoParam.kcRight,stereoParam.alphaRight);
        [xl] = normalize_pixel(cbcXYLL{(i)},stereoParam.focLeft,stereoParam.cenLeft,stereoParam.kcLeft,stereoParam.alphaLeft);
        
        
        % ת����������
        xrr = pflat((KR)*pextend(xr));
        xll = pflat((KL)*pextend(xl));
        
        % �Խǵ���rectification
        pixRectR = Orig2Rect(cbcXYRR{(i)}', KR, intrMatRightNew, rotMatRight, stereoParam.kcRight);
        pixRectL = Orig2Rect(cbcXYLL{(i)}', KL, intrMatLeftNew, rotMatLeft, stereoParam.kcLeft);
        
        
        
        
        if ~cfg.isRotate
            dispTemp= abs(pixRectL(:,1) - pixRectR(:,1));
            % ���㼫�߷�����Խ�ӽ�0Խ��
            epilineErr = abs(pixRectL(:,2) - pixRectR(:,2));
            cfg.Corner{i,1} = goodId;
            cfg.Corner{i,2} = pixRectL;
            cfg.Corner{i,3} = imgListL{i};
            cfg.CornerCheck{i,1} = imgListL{goodId(i)};
            
            
        else
            dispTemp= abs(pixRectL(:,2) - pixRectR(:,2));
            % ���㼫�߷�����Խ�ӽ�0Խ��
            epilineErr = abs(pixRectL(:,1) - pixRectR(:,1));
            cfg.Corner{i,1} = goodId;
            cfg.Corner{i,2} = pixRectR;
            cfg.Corner{i,3} = imgListR{i};
            cfg.CornerCheck{i,1} = imgListR{goodId(i)};
            
        end
        
        % ���Ӳ�ָ���ȣ������ָ�3d������
        depthTemp = intrMatLeftNew(1,1) * norm(stereoParam.transVecRef) ./ (dispTemp + (princpPtR(1) - princpPtL(1)));
        [XYZ_all] = GetXYZFromDepth(intrMatLeftNew, pixRectL, depthTemp);
        
        % ��3d����ָ���ÿ����ά��ˮƽ����ʹ�ֱ�����ϵı߳���������������Ķ�ά��߳��Ƚϣ��õ�ˮƽ����ֱ�����ϵ����
        X_all = reshape(XYZ_all(:,1), [], 4);
        Y_all = reshape(XYZ_all(:,2), [], 4);
        Z_all = reshape(XYZ_all(:,3), [], 4);
        XYZ1 = [X_all(:,1) Y_all(:,1) Z_all(:,1)];
        XYZ2 = [X_all(:,2) Y_all(:,2) Z_all(:,2)];
        XYZ3 = [X_all(:,3) Y_all(:,3) Z_all(:,3)];
        XYZ4 = [X_all(:,4) Y_all(:,4) Z_all(:,4)];
        [~, distHori1] = NormalizeVector(XYZ1 - XYZ2);
        [~, distHori2] = NormalizeVector(XYZ3 - XYZ4);
        [~, distVerti1] = NormalizeVector(XYZ1 - XYZ4);
        [~, distVerti2] = NormalizeVector(XYZ2 - XYZ3);
        horiErr = abs(mean([distHori1; distHori2]) - cfg.aruco_size);
        vertiErr = abs(mean([distVerti1; distVerti2]) - cfg.aruco_size);
        
        
        xyErr(i,:) = [horiErr vertiErr];
        
        if 0
            [rectImgLtemp, rectImgRtemp] = RectifyImagePair(stereoParam, imread(imgListL{i}), imread(imgListR{i}));
            figure,subplot(1,2,1);imshow(rectImgRtemp);hold on;plot(pixRectR(:,1), pixRectR(:,2),'.r'); subplot(1,2,2),imshow(imread(imgListR{i}));hold on;plot(cbcXYRR{(i)}(1,:), cbcXYRR{(i)}(2,:),'.r')
            figure,subplot(1,2,1);imshow(rectImgLtemp);hold on;plot(pixRectL(:,1), pixRectL(:,2),'.r'); subplot(1,2,2),imshow(imread(imgListL{i}));hold on;plot(cbcXYLL{(i)}(1,:), cbcXYLL{(i)}(2,:),'.r')
        end
        
        
        
        rt = [rodrigues(stereoParam.optPose(1:3,(i))) stereoParam.optPose(4:6,(i));0 0 0 1];
        
        
        
        
        
        
        pt_3d_aruto = cbGridL{i}';
        pt_3d_aruto(:,3) = 0;
        cbPt = rt(1:3,1:3)*pt_3d_aruto' + repmat(rt(1:3,4),1,size(pt_3d_aruto,1));
        
        
        ptIcs = TransformAndProject(cbPt', KR, rodrigues(stereoParam.rotVecRef), stereoParam.transVecRef);
        if 0
            figure(10),imshow(zeros(nr,nc));hold on;plot(xrr(1,:),xrr(2,:),'or');plot(ptIcs(:,1),ptIcs(:,2),'.g');
            drawnow;
            pause(0.5);
        end
        
        % ��Ŀ����Ŀ����ͶӰ���
        errMat = xrr(1:2,:)-ptIcs';
        err1(i,1) = mean(abs(errMat(:)));
        err(i,1) = norm(mean(abs(xrr(1:2,:)'-ptIcs)));
        
    end
catch
    err = 1.11111111;
    askjlh = 1;
end

meanError1 = mean(err);

meanError = [meanError1 mean(xyErr,1)];

fprintf(sprintf('\n\n\n### average reprojection error: %0.4f pixel ###\n\n\n',mean(err)));



end