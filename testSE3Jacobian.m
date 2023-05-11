function testSE3Jacobian()


% n = 10;
% 
% X = 10*randn(3,n);
% om = randn(3,1);
% T = [10*randn(2,1);40];
% f = 1000*rand(2,1);
% c = 1000*randn(2,1);
% k = 0.5*randn(5,1);
% alpha = 0.01*randn(1,1);
% 
% [x,dxdom,dxdT,dxdf,dxdc,dxdk,dxdalpha] = project_points2(X,om,T,f,c,k,alpha);
% 
% 
% % Test on om: OK
% 
% dom = 0.000000001 * norm(om)*randn(3,1);
% om2 = om + dom;
% 
% [x2] = project_points2(X,om2,T,f,c,k,alpha);
% 
% x_pred = x + reshape(dxdom * dom,2,n);
% 
% 
% norm(x2-x)/norm(x2 - x_pred)

pt3d0 = 1000*(rand(10,3)-0.5);
pt3d0(:,3) = 2000+2000*rand(10,1);

rotVec0 = randn(3,1);
transVec0 = [10;20;30];
% transVec0 = [0;0;0]


[pt2d0,dedt0, dedr0, dedp0] = projJac(pt3d0,rotVec0, transVec0);
dom = 0.000001 *randn(3,1);
rotVec = rotVec0 + dom;
dT = 0.00 * norm(transVec0)*randn(3,1);
transVec = transVec0 + dT;

[pt2d] = projJac(pt3d0,rotVec, transVec);

 x_pred1 = pt2d0 + reshape(dedt0 * [dT],2,10)';
errT = norm(x_pred1 - pt2d);
 x_pred2 = pt2d0 + reshape(dedr0 * [dom],2,10)';
 errR = norm(x_pred2 - pt2d);
 
 x_pred3 = pt2d0 + reshape(dedp0 * [dT;dom],2,10)';
 errRT = norm(x_pred3 - pt2d);
 
 [errT errR errRT]
 
%  x_pred = pt2d0 + reshape(dedr0 * [dom],2,10)';
%  norm(pt2d-pt2d0)/norm(pt2d - x_pred);

norm(x_pred1 - pt2d);

% x_pred = pt2d0 + reshape(dedp0 * [dT;dom],2,10)';
% norm(pt2d-pt2d0)/norm(pt2d - x_pred);
% norm(x_pred - pt2d)
end
function [pt2d, dedp] = projJac2(pt3d,rotVec, transVec)


intrMat = [500 0 320;0 500 240;0 0 1];

rotMat = rodrigues(rotVec);

pt3d_ =( rotMat*pt3d' + transVec)';
pt2d = pflat(intrMat*pt3d_');
pt2d = pt2d(1:2,:)';

fx = intrMat(1,1);
fy = intrMat(2,2);
dedp = [];
for i = 1 : size(pt3d,1)
    xyz = pt3d(i,:);
    xyz_ = pt3d_(i,:);
    dedp_ = -[fx/xyz_(3) 0 -fx*xyz_(1)/xyz_(3)/xyz_(3) -fx*xyz_(1)*xyz_(2)/xyz_(3)/xyz_(3) fx+fx*xyz_(1)*xyz(1)/xyz_(3)/xyz_(3) -fx*xyz_(2)/xyz_(3);
        0 fy/xyz_(3) -fy*xyz_(2)/xyz_(3)/xyz_(3)  -fy-fy*xyz_(2)*xyz_(2)/xyz_(3)/xyz_(3) fy*xyz_(1)*xyz_(2)/xyz_(3)/xyz_(3) fy*xyz_(1)/xyz_(3)];
    dedp = [dedp; dedp_];
end
end
function [pt2d, dedt, dedr, dedp] = projJac(pt3d,rotVec, transVec)


intrMat = [500 0 320;0 500 240;0 0 1];

rotMat = rodrigues(rotVec);

pt3d_ =( rotMat*pt3d' + transVec)';
pt2d = pflat(intrMat*pt3d_');
pt2d = pt2d(1:2,:)';

fx = intrMat(1,1);
fy = intrMat(2,2);



dedp = [];
dedt = [];
dedr = [];
for i = 1 : size(pt3d,1)
    xyz = pt3d(i,:);
    xyz_ = pt3d_(i,:);
    
    dpixdp = [fx/xyz_(3) 0 -fx*xyz_(1)/xyz_(3)/xyz_(3);0 fy/xyz_(3) -fy*xyz_(2)/xyz_(3)/xyz_(3)];
    
    dedp_ = [fx/xyz_(3) 0 -fx*xyz_(1)/xyz_(3)/xyz_(3) -fx*xyz_(1)*xyz_(2)/xyz_(3)/xyz_(3) fx+fx*xyz_(1)*xyz_(1)/xyz_(3)/xyz_(3) -fx*xyz_(2)/xyz_(3);
        0 fy/xyz_(3) -fy*xyz_(2)/xyz_(3)/xyz_(3)  -fy-fy*xyz_(2)*xyz_(2)/xyz_(3)/xyz_(3) fy*xyz_(1)*xyz_(2)/xyz_(3)/xyz_(3) fy*xyz_(1)/xyz_(3)];
    dedp = [dedp; dedp_];
    dedt_ = dpixdp*eye(3); %rotMat;
    dedt = [dedt; dedt_];
    dedr_ =  dpixdp*(-SkewSymMat(rotMat*xyz'));
    dedr = [dedr; dedr_];
    
end
end