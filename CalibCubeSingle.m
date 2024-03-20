function [camParam, camParam1, cbcXY, cbGrid, cbGridUse, cbcXY_, cbGrid_, config, errCalib] = CalibCubeSingle(ptTrack, imgSize, config, poseGT, intrGT, kcGT)
% camParam = 1;
config = SetDefaultConfig(config);
estAspectRatio = config.estAspectRatio;

%整理匹配数据
[cbcXY, cbGrid, nImgC, nImgR] = ComposeData(ptTrack, imgSize);

%获得参数初值
[initFoc, initCen, initK, initAlpha, cbcXY_, cbGrid_, cbGridUse] = InitIntrinsic(cbcXY, cbGrid, nImgC, nImgR, estAspectRatio, config.gridSize);
[rotVec, tranVec, ~] = InitExtrinsic(cbcXY_, cbGrid_, initFoc, initCen, initK, initAlpha);


%整理数据
rtVec1 = [rotVec; tranVec];
initParam1 = [initFoc; initCen; initAlpha; initK; zeros(5,1); rtVec1(:)];
%用单面优化参数
param1 = OptimizeIter_(cbcXY_, cbGrid_, initParam1, nImgC, nImgR, config);
paramErr1 = EstimateError_(param1, cbcXY_, cbGrid_, config);
camParam1 = ConstructParam(param1, paramErr1, length(cbcXY_));
if 0
    if 1
        T1 = [rodrigues(rotVec(:,1)) tranVec(:,1); 0 0 0 1];
        T22 = [rodrigues(rotVec(:,2)) tranVec(:,2); 0 0 0 1];
        T33 = [rodrigues(rotVec(:,3)) tranVec(:,3); 0 0 0 1];
    else
        T11 = [rodrigues(rotVec(:,1)) tranVec(:,1); 0 0 0 1];
        T22 = [rodrigues(rotVec(:,2)) tranVec(:,2); 0 0 0 1];
        T3 = [rodrigues(rotVec(:,3)) tranVec(:,3); 0 0 0 1];
    end
    if 0
        T2 = T1*[rotx(-90) [0 0 0]';0 0 0 1];
        T3 = T2*[rotz(-90) [0 0 0]';0 0 0 1];
    elseif 0
        T2 = T1*[(rotx(-90)*rotz(-90))' [0 0 0]';0 0 0 1];
        T3 = T1*[(roty(90)*rotz(90))' [0 0 0]';0 0 0 1];
    elseif 1
        
        if 0
            T2 = inv([roty(90)*rotz(90) [0 0 0]';0 0 0 1]*inv(T1));
            T3 = inv([rotx(-90)*rotz(-90) [0 0 0]';0 0 0 1]*inv(T1));
        else
            if 0
                T2 = inv(inv(config.rt2)*inv(T1));
                T3 = inv(inv(config.rt3)*inv(T1));
            elseif 1
                T2 = T1*config.rt2;
                T3 = T1*config.rt3;
            else
                
                
                T1 = inv(config.rt3)*T3;
                T2 = inv(config.rt2)*config.rt3*T3;
            end
        end
        
        if 0
            err2 = [rotx(-90)*rotz(-90) [0 0 0]';0 0 0 1]*inv(T1)*pextend(pt3dL(1:49,:)') - pextend(marker2');
            err3 = [roty(90)*rotz(90) [0 0 0]';0 0 0 1]*inv(T1)*pextend(pt3dL(1:49,:)') - pextend(marker3');
        end
    else
        T2 = inv([rotx(-90)*rotz(-90) [0 0 0]';0 0 0 1]*inv(T1));
        T3 = inv([roty(90)*rotz(90) [0 0 0]';0 0 0 1]*inv(T1));
        
    end
    if 1
        dif = [T2 - T22; T3 - T33]
        diff = [{T2*inv(T22)}; {T3*inv(T33)}];
        err = [rad2deg(norm(rodrigues(diff{1}(1:3,1:3)))) norm(diff{1}(1:3,4)) rad2deg(norm(rodrigues(diff{2}(1:3,1:3)))) norm(diff{2}(1:3,4))];
    else
        dif = [T1 - T11; T2 - T22]
        diff = [{T1*inv(T11)}; {T2*inv(T22)}];
        err = [rad2deg(norm(rodrigues(diff{1}(1:3,1:3)))) norm(diff{1}(1:3,4)) rad2deg(norm(rodrigues(diff{2}(1:3,1:3)))) norm(diff{2}(1:3,4))];
        
    end
else
    dif = [];
    for j = 1 : length(cbcXY)
        T1 = [rodrigues(rotVec(:,1+3*(j-1))) tranVec(:,1+3*(j-1)); 0 0 0 1];
        T22 = [rodrigues(rotVec(:,2+3*(j-1))) tranVec(:,2+3*(j-1)); 0 0 0 1];
        T33 = [rodrigues(rotVec(:,3+3*(j-1))) tranVec(:,3+3*(j-1)); 0 0 0 1];
        T2 = T1*config.rt2;
        T3 = T1*config.rt3;
        dif = [dif [T2 - T22; T3 - T33]];
%         diff = [{T2*inv(T22)}; {T3*inv(T33)}];
%         err = [rad2deg(norm(rodrigues(diff{1}(1:3,1:3)))) norm(diff{1}(1:3,4)) rad2deg(norm(rodrigues(diff{2}(1:3,1:3)))) norm(diff{2}(1:3,4))];
    end
end
if 0
    dif
end
%  [rotVecRight,transVecRight,drrdrl, drrdtl, drrdrref,drrdtref,dtrdrl, dtrdtl, dtrdrref,dtrdtref] = compose_motion(rotVecLeft,                    transVecLeft,     rotVecRef,   transVecRef );
%测试用，算法里没用到
[rotVec2,    transVec2,    dr2dr12,dr2dt12,dr2dr1,  dr2dt1,  dt2dr12,dt2dt12,dt2dr1,  dt2dt1]   = compose_motion(rodrigues(config.rt2(1:3,1:3)),config.rt2(1:3,4),rotVec(:,1), tranVec(:,1));
[rotVec3,    transVec3,    dr3dr13,dr3dt13,dr3dr1,  dr3dt1,  dt3dr13,dt3dt13,dt3dr1,  dt3dt1]   = compose_motion(rodrigues(config.rt3(1:3,1:3)),config.rt3(1:3,4),rotVec(:,1), tranVec(:,1));

rtVec = [];
for j = 1 : length(cbcXY)
    rtVec = [rtVec [rotVec(:,1+3*(j-1)); tranVec(:,1+3*(j-1))]];
    
end
initParam = [initFoc; initCen; initAlpha; initK; zeros(5,1); rtVec(:)];
%用三面优化参数
param = OptimizeIterCube(cbcXY, cbGridUse, initParam, nImgC, nImgR, config);
%整理误差
intrDiffCat = ([[param(1:4)  param1(1:4)]] - [intrGT(1,1); intrGT(2,2); intrGT(1,3); intrGT(2,3)]);
diffCalib = [{[rodrigues(param(16:18)) param(19:21); 0 0 0 1]*inv(poseGT{1})}; {[rodrigues(param1(16:18)) param1(19:21); 0 0 0 1]*inv(poseGT{1})}];
extrDiffCat = [[rad2deg(norm(rodrigues(diffCalib{1}(1:3,1:3))));norm(diffCalib{1}(1:3,4))] [rad2deg(norm(rodrigues(diffCalib{2}(1:3,1:3)))); norm(diffCalib{2}(1:3,4))]];
kcDiffCat = [[param(6:10) - kcGT] [param1(6:10) - kcGT]];
errCalib = [intrDiffCat; extrDiffCat; kcDiffCat];


paramErr = EstimateError(param, cbcXY, cbGridUse, config);
%整理优化后的参数并输出
camParam = ConstructParam(param, paramErr, length(cbcXY));

end

function [cbcXY, cbGrid, nImgC, nImgR] = ComposeData(ptTrack, imgSize)
cbcXY = cell(length(ptTrack),1);
cbGrid = cell(length(ptTrack),1);
nImgR = imgSize(1);
nImgC = imgSize(2);
for i = 1 : length(ptTrack)
    
    for j = 1 : 3
        cbcXY{i,1} = [cbcXY{i,1} [ptTrack{i}{1}{j,2}(:,1:2)';j*ones(1,length(ptTrack{i}{1}{j,1}))]];
        cbGrid{i,1} = [cbGrid{i,1} [ptTrack{i}{1}{j,2}(:,3:5)';j*ones(1,length(ptTrack{i}{1}{j,1}))]];
    end
end




end
function [initFoc, initCen, initK, initAlpha, cbcXY, cbGrid, cbGridUse] = InitIntrinsic(cbcXY0, cbGrid0, nImgC, nImgR, estAspectRatio,gridSize)

cbcXY = {}; % cell(3*length(cbcXY0),1);
cbGrid = {}; % cell(3*length(cbGrid0),1);
% cbGridUse = {};
ofset = [0 gridSize^2 2*gridSize^2];
for i = 1 : length(cbcXY0)
    cbGridUse{i,1} = {};
    for j = 1 : 3
        idx = find(cbcXY0{i,1}(3,:) == j);
        cbcXY = [cbcXY; {cbcXY0{i,1}(1:2,idx)}];
        if 0
            if j == 1
                cbGrid = [cbGrid; {cbGrid0{i,1}([1 2],idx)}];
            elseif j == 2
                cbGrid = [cbGrid; {cbGrid0{i,1}([1 3],idx)}];
            else
                cbGrid = [cbGrid; {cbGrid0{i,1}([2 3],idx)}];
            end
            cbGridUse{i,1} = [cbGridUse{i,1}; [{idx} {cbGrid0{i,1}([1 2],idx - 0)}]];
        else
            if 0 % j == 1
                idx1 = idx;
            end
            cbGrid = [cbGrid; {cbGrid0{i,1}([1 2],idx - ofset(j))}];
            cbGridUse{i,1} = [cbGridUse{i,1}; [{idx} {cbGrid0{i,1}([1 2],idx - ofset(j))}]];
        end
    end
    
end

nImg = length(cbcXY);
% initialize at the center of the image
initCen = [nImgC; nImgR]/2-0.5;
% initialize to zero (no distortion)
initK = [0;0;0;0;0];
subCenMat = [1,0,-initCen(1)
    0,1,-initCen(2)
    0,0,1];

A = zeros(2*nImg, 2);
b = zeros(2*nImg, 1);
for iImg = 1:nImg
    homo = compute_homography(cbcXY{iImg}, cbGrid{iImg});
    % transform center to origin
    homo = subCenMat * homo;
    % Extract vanishing points (direct and diagonals)
    vanMat = homo * [1,0,0.5,0.5; 0,1,0.5,-0.5; 0,0,0,0];
    normVec = sqrt(sum(vanMat.^2));
    vanMat = vanMat./normVec(ones(3,1), :);
    A((iImg-1)*2+1:iImg*2, :) = [vanMat(1,1)*vanMat(1,2), vanMat(2,1)*vanMat(2,2); vanMat(1,3)*vanMat(1,4), vanMat(2,3)*vanMat(2,4)];
    b((iImg-1)*2+1:iImg*2) = -[vanMat(3,1)*vanMat(3,2); vanMat(3,3)*vanMat(3,4)];
end

% use all the vanishing points to estimate focal length
% Different fx, fy (two) or the same (one)
if (b'*(sum(A, 2)) < 0)
    % Use a two focals estimate:
    initFoc = sqrt(abs(1./((A'*A)\(A'*b)))); % if using a two-focal model for initial guess
else
    initFoc = sqrt(b'*(sum(A, 2))/(b'*b)) * ones(2,1); % if single focal length model is used
end

if ~estAspectRatio
    initFoc = mean(initFoc)*one(1,2);
end

initAlpha = 0;

initIntrinsicMat = [ ...
    initFoc(1)    initAlpha*initFoc(1)   initCen(1)
    0             initFoc(2)             initCen(2)
    0             0                      1];

end

function [rotVec, tranVec, rotMat] = InitExtrinsic(cbcXY, cbGrid, foc, cen, k, alpha)

condTh = 1e6; % threshold of conditional number

nImg = length(cbcXY);
rotVec = zeros(3, nImg);
tranVec = zeros(3, nImg);
rotMat = zeros(3,3,nImg);
for iImg = 1:nImg
    gridHomo = [cbGrid{iImg}; zeros(1, size(cbGrid{iImg}, 2))];
    [initRotVec,initTranVec] = compute_extrinsic_init(cbcXY{iImg},gridHomo,foc,cen,k,alpha);
    [rotVec(:,iImg),tranVec(:,iImg),rotMat(:,:,iImg),jacMat] = ...
        compute_extrinsic_refine(initRotVec,initTranVec,cbcXY{iImg},gridHomo,foc,cen,k,alpha,20,condTh);
    if (cond(jacMat)> condTh)
        error('View #%d ill-conditioned', iImg);
    end;
end

end

function param = OptimizeIterCube(cbcXY, cbGrid, initParam, imgW, imgH, config)

condTh = 1e6; % threshold of conditional number
maxIter = 30; 6; 8; 10; 30;
change = 1;
iter = 0;
alpha_smooth = 0.1;  % 0.05
estFocalLen = config.estFocalLen;
centerOptim = config.centerOptim;
estAlpha = config.estAlpha;
estDistortion = config.estDistortion;
param = initParam;

Change = [];
nImg = length(cbcXY);
selected_variables = [estFocalLen; centerOptim*ones(2,1); estAlpha; estDistortion; zeros(5,1); ...
    ones(6*nImg, 1)];
ind_Jac = find(selected_variables)';
while (change > 1e-9 && iter < maxIter)
    % while (change > 1e-5 && iter < maxIter)
    foc = param(1:2);
    cen = param(3:4);
    alpha = param(5);
    k = param(6:10);
    
    JJ3 = sparse([],[],[], 15+6*nImg, 15+6*nImg);
    ex3 = zeros(15+6*nImg, 1);
    
    % The first step consists of updating the whole vector of knowns
    % (intrinsic + extrinsic of active images) through a one step steepest
    % gradient descent.
    for iImg = 1:nImg
        for j = 1 : 1
            gridHomo = [cell2mat(cbGrid{iImg}(:,2)'); zeros(1, size(cell2mat(cbGrid{iImg}(:,2)'), 2))];
            r = param(15+6*(iImg-1) + 1:15+6*(iImg-1) + 3);
            t = param(15+6*(iImg-1) + 4:15+6*(iImg-1) + 6);
            T1 = [rodrigues(r) t; 0 0 0 1];
            %% important function
            cubeIdx = cbGrid{iImg}(:,1);
            if 1
                %构造三面残差约束
                [xPrj,dxdr,dxdt,dxdf,dxdc,dxdk,dxdalpha] = project_points_cube_flow(iImg, cbGrid, r,t,foc,cen,k,alpha,cubeIdx,config);
                %                 [xPrj,dxdr,dxdt,dxdf,dxdc,dxdk,dxdalpha] = project_points_cube2(gridHomo,r,t,foc,cen,k,alpha,cubeIdx,config);
                xErr = cbcXY{iImg}(1:2,:) - xPrj;
            else
                gridHomoi = [cbGrid{iImg}{iImg,2}; zeros(1, length(cubeIdx{iImg}))];
                if j == 1
                    T = T1;
                elseif j == 2
                    T = T1*config.rt2;
                else
                    T = T1*config.rt3;
                end
                [xPrj,dxdr,dxdt,dxdf,dxdc,dxdk,dxdalpha] = project_points2(gridHomoi,rodrigues(T(1:3,1:3)),T(1:3,4),foc,cen,k,alpha);
                xErr = cbcXY{iImg}(1:2,cubeIdx{j}) - xPrj;
            end
            
            A = [dxdf dxdc dxdalpha dxdk]';
            B = [dxdr dxdt]';
            if (cond(B') > condTh)
                error('View #%d ill-conditioned', iImg);
            end
            JJ3(1:10,1:10) = JJ3(1:10,1:10) + sparse(A*A');
            JJ3(15+6*(iImg-1) + 1:15+6*(iImg-1) + 6,15+6*(iImg-1) + 1:15+6*(iImg-1) + 6) = sparse(B*B');
            
            AB = sparse(A*B');
            JJ3(1:10,15+6*(iImg-1) + 1:15+6*(iImg-1) + 6) = AB;
            JJ3(15+6*(iImg-1) + 1:15+6*(iImg-1) + 6,1:10) = (AB)';
            
            ex3(1:10) = ex3(1:10) + A*xErr(:);
            ex3(15+6*(iImg-1) + 1:15+6*(iImg-1) + 6) = B*xErr(:);
        end
    end
    
    JJ3 = JJ3(ind_Jac,ind_Jac);
    ex3 = ex3(ind_Jac);
    alpha_smooth2 = 1-(1-alpha_smooth)^(iter+1);
    param_innov = alpha_smooth2*(JJ3\ex3);
    param(ind_Jac) = param(ind_Jac) + param_innov;
    
    % Second step: (optional) - It makes convergence faster, and the
    % region of convergence LARGER!!!
    % Recompute the extrinsic parameters only using compute_extrinsic.m
    % (this may be useful sometimes)
    % The complete gradient descent method is useful to precisely update
    % the intrinsic parameters.
    focCurr = param(1:2);
    cenCurr = param(3:4);
    if 0% (centerOptim && (cenCurr(1) < 0 || cenCurr(1) > imgW || cenCurr(2) < 0 || cenCurr(2) > imgH))
        error('Principal point cannot be estimated.');
    end
    
    alphaCurr = param(5);
    kcCurr = param(6:10);
    change = norm([focCurr;cenCurr] - [foc;cen])/norm([focCurr;cenCurr]);
    Change = [Change; change];
    
    % Recompute extrinsic
    for iImg = 1:nImg
        cubeIdx = cbGrid{iImg}(:,1);
        gridHomo1 = [cbGrid{iImg}{1,2}; zeros(1, length(cubeIdx{1}))];
        gridHomo = [cell2mat(cbGrid{iImg}(:,2)'); zeros(1, size(cell2mat(cbGrid{iImg}(:,2)'), 2))];
        [r, t] = compute_extrinsic_init(cbcXY{iImg}(1:2,cubeIdx{1}),gridHomo1,focCurr,cenCurr,kcCurr,alphaCurr);
        rt0 = [r t];
        %         [r, t] = compute_extrinsic_cube_init(cbcXY{iImg}(1:2,cubeIdx{1}),gridHomo1,focCurr,cenCurr,kcCurr,alphaCurr);
        %         [r, t, ~, jacMat] = compute_extrinsic_cube_refine(r, t, cbcXY{iImg}(1:2,:), gridHomo, focCurr,cenCurr,kcCurr,alphaCurr,20,condTh,cubeIdx,config);
        [r, t, ~, jacMat] = compute_extrinsic_cube_refine(r, t, cbcXY{iImg}(1:2,:), gridHomo, focCurr,cenCurr,kcCurr,alphaCurr,20,condTh,cubeIdx,config, iImg, cbGrid);
        if (cond(jacMat)> condTh)
            error('View #%d ill-conditioned', iImg);
        end;
        param(15+6*(iImg-1) + 1:15+6*(iImg-1) + 3) = r;
        param(15+6*(iImg-1) + 4:15+6*(iImg-1) + 6) = t;
    end
    
    iter = iter + 1;
end

end

function param = OptimizeIter_(cbcXY, cbGrid, initParam, imgW, imgH, config)

condTh = 1e6; % threshold of conditional number
maxIter = 30;
change = 1;
iter = 0;
alpha_smooth = 0.1;
estFocalLen = config.estFocalLen;
centerOptim = config.centerOptim;
estAlpha = config.estAlpha;
estDistortion = config.estDistortion;
param = initParam;


nImg = length(cbcXY);
selected_variables = [estFocalLen; centerOptim*ones(2,1); estAlpha; estDistortion; zeros(5,1); ...
    ones(6*nImg, 1)];
ind_Jac = find(selected_variables)';
while (change > 1e-9 && iter < maxIter)
    foc = param(1:2);
    cen = param(3:4);
    alpha = param(5);
    k = param(6:10);
    
    JJ3 = sparse([],[],[], 15+6*nImg, 15+6*nImg);
    ex3 = zeros(15+6*nImg, 1);
    
    % The first step consists of updating the whole vector of knowns
    % (intrinsic + extrinsic of active images) through a one step steepest
    % gradient descent.
    for iImg = 1:nImg
        gridHomo = [cbGrid{iImg}; zeros(1, size(cbGrid{iImg}, 2))];
        r = param(15+6*(iImg-1) + 1:15+6*(iImg-1) + 3);
        t = param(15+6*(iImg-1) + 4:15+6*(iImg-1) + 6);
        %% important function
        [xPrj,dxdr,dxdt,dxdf,dxdc,dxdk,dxdalpha] = project_points2(gridHomo,r,t,foc,cen,k,alpha);
        xErr = cbcXY{iImg} - xPrj;
        A = [dxdf dxdc dxdalpha dxdk]';
        B = [dxdr dxdt]';
        if (cond(B') > condTh)
            error('View #%d ill-conditioned', iImg);
        end
        JJ3(1:10,1:10) = JJ3(1:10,1:10) + sparse(A*A');
        JJ3(15+6*(iImg-1) + 1:15+6*(iImg-1) + 6,15+6*(iImg-1) + 1:15+6*(iImg-1) + 6) = sparse(B*B');
        
        AB = sparse(A*B');
        JJ3(1:10,15+6*(iImg-1) + 1:15+6*(iImg-1) + 6) = AB;
        JJ3(15+6*(iImg-1) + 1:15+6*(iImg-1) + 6,1:10) = (AB)';
        
        ex3(1:10) = ex3(1:10) + A*xErr(:);
        ex3(15+6*(iImg-1) + 1:15+6*(iImg-1) + 6) = B*xErr(:);
    end
    
    JJ3 = JJ3(ind_Jac,ind_Jac);
    ex3 = ex3(ind_Jac);
    alpha_smooth2 = 1-(1-alpha_smooth)^(iter+1);
    param_innov = alpha_smooth2*(JJ3\ex3);
    param(ind_Jac) = param(ind_Jac) + param_innov;
    
    % Second step: (optional) - It makes convergence faster, and the
    % region of convergence LARGER!!!
    % Recompute the extrinsic parameters only using compute_extrinsic.m
    % (this may be useful sometimes)
    % The complete gradient descent method is useful to precisely update
    % the intrinsic parameters.
    focCurr = param(1:2);
    cenCurr = param(3:4);
    if 0% (centerOptim && (cenCurr(1) < 0 || cenCurr(1) > imgW || cenCurr(2) < 0 || cenCurr(2) > imgH))
        error('Principal point cannot be estimated.');
    end
    
    alphaCurr = param(5);
    kcCurr = param(6:10);
    change = norm([focCurr;cenCurr] - [foc;cen])/norm([focCurr;cenCurr]);
    
    % Recompute extrinsic
    for iImg = 1:nImg
        gridHomo = [cbGrid{iImg}; zeros(1, size(cbGrid{iImg}, 2))];
        [r, t] = compute_extrinsic_init(cbcXY{iImg},gridHomo,focCurr,cenCurr,kcCurr,alphaCurr);
        [r, t, ~, jacMat] = compute_extrinsic_refine(r, t, cbcXY{iImg}, gridHomo, focCurr,cenCurr,kcCurr,alphaCurr,20,condTh);
        if (cond(jacMat)> condTh)
            error('View #%d ill-conditioned', iImg);
        end;
        param(15+6*(iImg-1) + 1:15+6*(iImg-1) + 3) = r;
        param(15+6*(iImg-1) + 4:15+6*(iImg-1) + 6) = t;
    end
    
    iter = iter + 1;
end

end

function paramErr = EstimateError(param, cbcXY, cbGrid, config)

foc = param(1:2);
cen = param(3:4);
alpha = param(5);
kc = param(6:10);
estFocalLen = config.estFocalLen;
centerOptim = config.centerOptim;
estAlpha = config.estAlpha;
estDistortion = config.estDistortion;

nImg = length(cbcXY);
xErr = [];

JJ3 = sparse([],[],[], 15+6*nImg, 15+6*nImg);
selected_variables = [estFocalLen; centerOptim*ones(2,1); estAlpha; estDistortion; zeros(5,1); ...
    ones(6*nImg, 1)];
ind_Jac = find(selected_variables)';
for iImg = 1:nImg
%     gridHomo = [cbGrid{iImg}; zeros(1, size(cbGrid{iImg}, 2))];
    rotVec = param(15+6*(iImg-1) + 1:15+6*(iImg-1) + 3);
    tranVec = param(15+6*(iImg-1) + 4:15+6*(iImg-1) + 6);
    cubeIdx = cbGrid{iImg}(:,1);
    if 0
        [xPrj,dxdr,dxdt,dxdf,dxdc,dxdk,dxdalpha] = project_points2(gridHomo,rotVec,tranVec,foc,cen,kc,alpha);
        xErr = [xErr, cbcXY{iImg} - xPrj];
    else
        [xPrj,dxdr,dxdt,dxdf,dxdc,dxdk,dxdalpha] = project_points_cube_flow(iImg, cbGrid, rotVec,tranVec,foc,cen,kc,alpha,cubeIdx,config);
        xErr = [xErr, cbcXY{iImg}(1:2,:) - xPrj];
    end
    
    
    A = [dxdf dxdc dxdalpha dxdk]';
    B = [dxdr dxdt]';
    JJ3(1:10,1:10) = JJ3(1:10,1:10) + sparse(A*A');
    JJ3(15+6*(iImg-1) + 1:15+6*(iImg-1) + 6,15+6*(iImg-1) + 1:15+6*(iImg-1) + 6) = sparse(B*B');
    
    AB = sparse(A*B');
    JJ3(1:10,15+6*(iImg-1) + 1:15+6*(iImg-1) + 6) = AB;
    JJ3(15+6*(iImg-1) + 1:15+6*(iImg-1) + 6,1:10) = (AB)';
end

sigma_x = std(xErr(:));

JJ3 = JJ3(ind_Jac,ind_Jac);

JJ2_inv = inv(JJ3); % not bad for sparse matrices!!

paramErr = zeros(6*nImg+15,1);
paramErr(ind_Jac) =  3*sqrt(full(diag(JJ2_inv)))*sigma_x;

end

function paramErr = EstimateError_(param, cbcXY, cbGrid, config)

foc = param(1:2);
cen = param(3:4);
alpha = param(5);
kc = param(6:10);
estFocalLen = config.estFocalLen;
centerOptim = config.centerOptim;
estAlpha = config.estAlpha;
estDistortion = config.estDistortion;

nImg = length(cbcXY);
xErr = [];

JJ3 = sparse([],[],[], 15+6*nImg, 15+6*nImg);
selected_variables = [estFocalLen; centerOptim*ones(2,1); estAlpha; estDistortion; zeros(5,1); ...
    ones(6*nImg, 1)];
ind_Jac = find(selected_variables)';
for iImg = 1:nImg
    gridHomo = [cbGrid{iImg}; zeros(1, size(cbGrid{iImg}, 2))];
    rotVec = param(15+6*(iImg-1) + 1:15+6*(iImg-1) + 3);
    tranVec = param(15+6*(iImg-1) + 4:15+6*(iImg-1) + 6);
    
    [xPrj,dxdr,dxdt,dxdf,dxdc,dxdk,dxdalpha] = project_points2(gridHomo,rotVec,tranVec,foc,cen,kc,alpha);
    xErr = [xErr, cbcXY{iImg} - xPrj];
    
    A = [dxdf dxdc dxdalpha dxdk]';
    B = [dxdr dxdt]';
    JJ3(1:10,1:10) = JJ3(1:10,1:10) + sparse(A*A');
    JJ3(15+6*(iImg-1) + 1:15+6*(iImg-1) + 6,15+6*(iImg-1) + 1:15+6*(iImg-1) + 6) = sparse(B*B');
    
    AB = sparse(A*B');
    JJ3(1:10,15+6*(iImg-1) + 1:15+6*(iImg-1) + 6) = AB;
    JJ3(15+6*(iImg-1) + 1:15+6*(iImg-1) + 6,1:10) = (AB)';
end

sigma_x = std(xErr(:));

JJ3 = JJ3(ind_Jac,ind_Jac);

JJ2_inv = inv(JJ3); % not bad for sparse matrices!!

paramErr = zeros(6*nImg+15,1);
paramErr(ind_Jac) =  3*sqrt(full(diag(JJ2_inv)))*sigma_x;

end

function camParam = ConstructParam(param, paramErr, nImg)
    camParam.foc = param(1:2);
    camParam.focErr = paramErr(1:2);
    camParam.cen = param(3:4);
    camParam.cenErr = paramErr(3:4);
    camParam.alpha = param(5);
    camParam.alphaErr = paramErr(5);
    camParam.kc = param(6:10);
    camParam.kcErr = paramErr(6:10);
    % extrinc param and error
    camParam.rotVec = zeros(3,nImg);
    camParam.rotVecErr = zeros(3,nImg);
    camParam.tranVec = zeros(3,nImg);
    camParam.tranVecErr = zeros(3,nImg);
    for iImg = 1:nImg
        camParam.rotVec(:,iImg) = param(15+6*(iImg-1) + 1:15+6*(iImg-1) + 3);
        camParam.rotVecErr(:,iImg) = paramErr(15+6*(iImg-1) + 1:15+6*(iImg-1) + 3);
        camParam.tranVec(:,iImg) = param(15+6*(iImg-1) + 4:15+6*(iImg-1) + 6);
        camParam.tranVecErr(:,iImg) = paramErr(15+6*(iImg-1) + 4:15+6*(iImg-1) + 6);
    end
end

function config = SetDefaultConfig(config)

% dX and dY is the width and height in mm of a white or black rectangle
dxdy = isfield(config, {'dX', 'dY'});
if all(~dxdy)
    config.dX = 100;
    config.dY = 100;
elseif (dxdy(1) && ~dxdy(2))
    config.dY = config.dX;
elseif (dxdy(2) && ~dxdy(1))
    config.dX = config.dY;
end

% halfWinW and halfWinH are the width and height of the rectangular area
% for corner search
wh = isfield(config, {'halfWinW', 'halfWinH'});
if all(~wh)
    config.halfWinW = 5;
    config.halfWinH = 5;
elseif (wh(1) && ~wh(2))
    config.halfWinH = config.halfWinW;
elseif (wh(2) && ~wh(1))
    config.halfWinW = config.halfWinH;
end

if ~isfield(config, 'estAspectRatio')
    config.estAspectRatio = 1;
end

if ~isfield(config, 'estFocalLen')
    config.estFocalLen = [1;1];
end

if ~isfield(config, 'centerOptim')
    config.centerOptim = 1;
end

if ~isfield(config, 'estAlpha')
    config.estAlpha = 0;
end

if ~isfield(config, 'estDistortion')
    config.estDistortion = [1;1;1;1;0];
    %         config.estDistortion = [1;1;1;1;1];
end

end