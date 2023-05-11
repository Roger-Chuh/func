function [camParam, cbcXY, cbGrid, config] = CalibArucoSingle(imgList, ptTrack, imgSize, config)
config = SetDefaultConfig(config);
estAspectRatio = config.estAspectRatio;
nImg = length(ptTrack);
[cbcXY, cbGrid, nImgC, nImgR] = ComposeData(ptTrack, imgSize);
[initFoc, initCen, initK, initAlpha] = InitIntrinsic(cbcXY, cbGrid, nImgC, nImgR, estAspectRatio);
[rotVec, tranVec, ~] = InitExtrinsic(cbcXY, cbGrid, initFoc, initCen, initK, initAlpha);

rtVec = [rotVec; tranVec];
initParam = [initFoc; initCen; initAlpha; initK; zeros(5,1); rtVec(:)];
try
    param = OptimizeIter(cbcXY, cbGrid, initParam, nImgC, nImgR, config);
    
    paramErr = EstimateError(param, cbcXY, cbGrid, config);
    
    % intrinsic param and error
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
catch
    camParam = [];
end

end

function [cbcXY, cbGrid, nImgC, nImgR] = ComposeData(ptTrack, imgSize)
cbcXY = cell(length(ptTrack),1);
cbGrid = cell(length(ptTrack),1);
nImgR = imgSize(1);
nImgC = imgSize(2);
for i = 1 : length(ptTrack)
    
    for j = 1 : 4
    cbcXY{i,1} = [cbcXY{i,1} ptTrack{i}{2}{j}(:,1:2)'];
    cbGrid{i,1} = [cbGrid{i,1} ptTrack{i}{2}{j}(:,3:4)'];
    end
end




end
function [initFoc, initCen, initK, initAlpha] = InitIntrinsic(cbcXY, cbGrid, nImgC, nImgR, estAspectRatio)

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

function param = OptimizeIter(cbcXY, cbGrid, initParam, imgW, imgH, config)

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



noSparse = 0;1;


nImg = length(cbcXY);
selected_variables = [estFocalLen; centerOptim*ones(2,1); estAlpha; estDistortion; zeros(5,1); ...
    ones(6*nImg, 1)];
ind_Jac = find(selected_variables)';
while (change > 1e-9 && iter < maxIter)
    foc = param(1:2);
    cen = param(3:4);
    alpha = param(5);
    k = param(6:10);
    if ~noSparse
        JJ3 = sparse([],[],[], 15+6*nImg, 15+6*nImg);
    else
        JJ3 = zeros( 15+6*nImg, 15+6*nImg);
    end
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
        if  ~noSparse
            JJ3(1:10,1:10) = JJ3(1:10,1:10) + sparse(A*A');
            JJ3(15+6*(iImg-1) + 1:15+6*(iImg-1) + 6,15+6*(iImg-1) + 1:15+6*(iImg-1) + 6) = sparse(B*B');
        else
            JJ3(1:10,1:10) = JJ3(1:10,1:10) + (A*A');
            JJ3(15+6*(iImg-1) + 1:15+6*(iImg-1) + 6,15+6*(iImg-1) + 1:15+6*(iImg-1) + 6) = (B*B');
        end
        
        if ~noSparse
            AB = sparse(A*B');
        else
            AB = (A*B');
        end
        JJ3(1:10,15+6*(iImg-1) + 1:15+6*(iImg-1) + 6) = AB;
        JJ3(15+6*(iImg-1) + 1:15+6*(iImg-1) + 6,1:10) = (AB)';
        
        ex3(1:10) = ex3(1:10) + A*xErr(:);
        ex3(15+6*(iImg-1) + 1:15+6*(iImg-1) + 6) = B*xErr(:);
    end
    JJ33 = JJ3;
    ex33 = ex3;
    JJ3 = JJ3(ind_Jac,ind_Jac);
    ex3 = ex3(ind_Jac);
    alpha_smooth2 = 1-(1-alpha_smooth)^(iter+1);
    param_innov = alpha_smooth2*(JJ3\ex3);
%     param_innov0 = alpha_smooth2*(JJ33\ex33);
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