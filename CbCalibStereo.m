function stereoParam = CbCalibStereo(paramLeft, paramRight, cbcXYLeft, cbcXYRight, cbGrid, ...
    configLeft, configRight)
% CbCalibStereo calibrate a pair of cameras using results from
% CbCalibSingle
% paramLeft and paramRight are structures including intrinsic and N set of
% extrinsic parameters of left and right cameras separately, where N is the
% number of images used for calibrating one camera.
% cbcXYLeft and cbcXYRight are two cell arrays each of which N arrays of
% image coordinates (0 based) of chessboard coner points
% cbGrid is a cell array which is the rectanglur coordinates (0 based) of
% checkerboard corner points
% configLeft and configRight are configurations used for single calibration
% of left and right cameras
% By Ji Zhou

nImgPair = length(cbGrid);
if (length(cbcXYLeft) ~= nImgPair)
    error('Number of images for left image calibration should be the number of chessboard grids');
end
if (length(cbcXYRight) ~= nImgPair)
    error('Number of images for right image calibration should be the number of chessboard grids');
end

rotVecRef = zeros(3, nImgPair);
transVecRef = zeros(3, nImgPair);
for iImgPair = 1:nImgPair
    try
        rotVecLeft = paramLeft.rotVec(:, iImgPair);
    catch
        fhav = 1;
    end
    rotMatLeft = rodrigues(rotVecLeft);
    transVecLeft = paramLeft.tranVec(:, iImgPair);
    rotVecRight = paramRight.rotVec(:, iImgPair);
    rotMatRight = rodrigues(rotVecRight);
    transVecRight = paramRight.tranVec(:, iImgPair);
    
    rotMatRef = rotMatRight * rotMatLeft';
    transVecRef(:, iImgPair) = transVecRight - rotMatRef * transVecLeft;
    rotVecRef(:, iImgPair) = rodrigues(rotMatRef);
end
rotVecRef = median(rotVecRef, 2);
transVecRef = median(transVecRef, 2);

% % rotVecRef = [0;0;0];
% % transVecRef = [0;0;0];

focLeft = paramLeft.foc;
cenLeft = paramLeft.cen;
alphaLeft = paramLeft.alpha;
kcLeft = paramLeft.kc;
focRight = paramRight.foc;
cenRight = paramRight.cen;
alphaRight = paramRight.alpha;
kcRight = paramRight.kc;
splitId = 27;
if 1
    param = [focLeft;cenLeft;alphaLeft;kcLeft; ...
            focRight;cenRight;alphaRight;kcRight; ...
            rotVecRef;transVecRef; ...
            reshape([paramLeft.rotVec; paramLeft.tranVec],[],1)];
else
    param = [focLeft;cenLeft;alphaLeft;kcLeft; ...
            focRight;cenRight;alphaRight;kcRight; ...
            [0;0;0];[0;0;0]; ...
            reshape([paramLeft.rotVec; paramLeft.tranVec],[],1)];
end
if 1
    estAspectRatioLeft = configLeft.estAspectRatio;     
    estFocalLenLeft = configLeft.estFocalLen;
    centerOptimLeft = configLeft.centerOptim;
    estAlphaLeft = configLeft.estAlpha;
    estDistortionLeft = configLeft.estDistortion;
    estAspectRatioRight = configRight.estAspectRatio;     
    estFocalLenRight = configRight.estFocalLen;
    centerOptimRight = configRight.centerOptim;
    estAlphaRight = configRight.estAlpha;
    estDistortionRight = configRight.estDistortion;
else
    estAspectRatioLeft = 0;  %configLeft.estAspectRatio;     
    estFocalLenLeft = [0 0]';  %configLeft.estFocalLen;
    centerOptimLeft = 0;  %configLeft.centerOptim;
    estAlphaLeft = configLeft.estAlpha;
    estDistortionLeft = [0 0 0 0 0]';  %configLeft.estDistortion;
    estAspectRatioRight = 0; %configRight.estAspectRatio;     
    estFocalLenRight = [0 0]';  %configRight.estFocalLen;
    centerOptimRight = 0';  %configRight.centerOptim;
    estAlphaRight = configRight.estAlpha;
    estDistortionRight = [0 0 0 0 0]';  %configRight.estDistortion;
end

ind_Jac = find([estFocalLenLeft & [1;estAspectRatioLeft];centerOptimLeft*ones(2,1);estAlphaLeft;estDistortionLeft;...
                estFocalLenRight & [1;estAspectRatioRight];centerOptimRight*ones(2,1);estAlphaRight;estDistortionRight;...
                ones(6,1);ones(6*nImgPair,1)]);

maxIter = 100;
change = 1;
iter = 1;
threshold = 50;
while ((change > 5e-6) && (iter <= maxIter))
    focLeft = param(1:2);
    cenLeft = param(3:4);
    alphaLeft = param(5);
    kcLeft = param(6:10);
    focRight = param(11:12);
    cenRight = param(13:14);
    alphaRight = param(15);
    kcRight = param(16:20);
    
    rotVecRefOld = param(1+20:3+20);
    transVecRefOld = param(4+20:6+20);
    
    J = [];
    e = [];

    for iImgPair = 1:nImgPair
        nPt = size(cbGrid{iImgPair}, 2);
        Jkk = sparse(4*nPt,20+(1+nImgPair)*6);
        ekk = zeros(4*nPt,1);
        
        % Project the structure onto the left view:
        rotVecLeft = param(6*(iImgPair-1)+7+20:6*(iImgPair-1)+7+2+20);
        transVecLeft = param(6*(iImgPair-1)+7+3+20:6*(iImgPair-1)+7+5+20);
        gridHomo = [cbGrid{iImgPair}; zeros(1, nPt)];
        [xl,dxldrl,dxldtl,dxldfl,dxldcl,dxldkl,dxldalphal] = project_points2(gridHomo,rotVecLeft,transVecLeft,focLeft,cenLeft,kcLeft,alphaLeft);
        
        ekk(1:2*nPt) = cbcXYLeft{iImgPair}(:) - xl(:);
        Jkk(1:2*nPt,6*(iImgPair-1)+7+20:6*(iImgPair-1)+7+2+20) = sparse(dxldrl);
        Jkk(1:2*nPt,6*(iImgPair-1)+7+3+20:6*(iImgPair-1)+7+5+20) = sparse(dxldtl);

        Jkk(1:2*nPt,1:2) = sparse(dxldfl);
        Jkk(1:2*nPt,3:4) = sparse(dxldcl);
        Jkk(1:2*nPt,5) = sparse(dxldalphal);
        Jkk(1:2*nPt,6:10) = sparse(dxldkl);
        
        % Project the structure onto the right view:
        [rotVecRight,transVecRight,drrdrl,drrdtl,drrdrref,drrdtref,dtrdrl,dtrdtl,dtrdrref,dtrdtref] = compose_motion(rotVecLeft,transVecLeft,rotVecRef,transVecRef);
        [xr,dxrdrr,dxrdtr,dxrdfr,dxrdcr,dxrdkr,dxrdalphar] = project_points2(gridHomo,rotVecRight,transVecRight,focRight,cenRight,kcRight,alphaRight);
        
        ekk(2*nPt+1:end) = cbcXYRight{iImgPair}(:) - xr(:);
        emax = max(abs(ekk));
        if (emax > threshold)
            error('Inconsistent image pairs #%d', iImgPair);
        end
        
        dxrdrref = dxrdrr * drrdrref + dxrdtr * dtrdrref;
        dxrdtref = dxrdrr * drrdtref + dxrdtr * dtrdtref;

        dxrdrl = dxrdrr * drrdrl + dxrdtr * dtrdrl;
        dxrdtl = dxrdrr * drrdtl + dxrdtr * dtrdtl;
        
        Jkk(2*nPt+1:end,1+20:3+20) =  sparse(dxrdrref);
        Jkk(2*nPt+1:end,4+20:6+20) =  sparse(dxrdtref);


        Jkk(2*nPt+1:end,6*(iImgPair-1)+7+20:6*(iImgPair-1)+7+2+20) = sparse(dxrdrl);
        Jkk(2*nPt+1:end,6*(iImgPair-1)+7+3+20:6*(iImgPair-1)+7+5+20) = sparse(dxrdtl);

        Jkk(2*nPt+1:end,11:12) = sparse(dxrdfr);
        Jkk(2*nPt+1:end,13:14) = sparse(dxrdcr);
        Jkk(2*nPt+1:end,15) = sparse(dxrdalphar);
        Jkk(2*nPt+1:end,16:20) = sparse(dxrdkr);
        
        J = [J; Jkk];
        e = [e; ekk];
    end
    
    J = J(:,ind_Jac);
    J2 = J'*J;
    J2_inv = inv(J2);
    
    param_update = J2_inv*J'*e;
    param(ind_Jac) = param(ind_Jac) + param_update;
    
    rotVecRef = param(1+20:3+20);
    transVecRef = param(4+20:6+20);
    change = norm([transVecRef;rotVecRef] - [transVecRefOld;rotVecRefOld])/norm([transVecRef;rotVecRef]);
    
    iter = iter+1;
end

% Computation of the error of estimation:

sigma_x = std(e(:));
paramError = zeros(20 + (1+nImgPair)*6,1);
paramError(ind_Jac) =  3*sqrt(full(diag(J2_inv)))*sigma_x;

% Extract parameters
stereoParam.focLeft = param(1:2);
stereoParam.cenLeft = param(3:4);
stereoParam.alphaLeft = param(5);
stereoParam.kcLeft = param(6:10);
stereoParam.focRight = param(11:12);
stereoParam.cenRight = param(13:14);
stereoParam.alphaRight = param(15);
stereoParam.kcRight = param(16:20);
stereoParam.rotVecRef = param(21:23);
stereoParam.transVecRef = param(24:26);

stereoParam.focLeftError = paramError(1:2);
stereoParam.cenLeftError = paramError(3:4);
stereoParam.alphaLeftError = paramError(5);
stereoParam.kcLeftError = paramError(6:10);
stereoParam.focRightError = paramError(11:12);
stereoParam.cenRightError = paramError(13:14);
stereoParam.alphaRightError = paramError(15);
stereoParam.kcRightError = paramError(16:20);
stereoParam.rotVecRefError = paramError(21:23);
stereoParam.transVecRefError = paramError(24:26);
stereoParam.optPose = reshape(param(splitId : end),6,[]);
end

