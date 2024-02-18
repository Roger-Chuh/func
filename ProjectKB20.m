function [pt2d, d_uv_d_pt3d, d_uv_d_param] = ProjectKB20(X,Y,Z, fx,fy,cx,cy,k1,k2,k3,k4,k5,k6,p1,p2,s1,s2,s3,s4,s5,s6,t1,t2)
d_uv_d_pt3d = []; d_uv_d_param = [];
param = [fx,fy,cx,cy,k1,k2,k3,k4,k5,k6,p1,p2,s1,s2,s3,s4,s5,s6,t1,t2];
% fx = param(1);
% fy = param(2);
% cx = param(3);
% cy = param(4);
% k1 = param(5);
% k2 = param(6);
% k3 = param(7);
% k4 = param(8);
% k5 = param(9);
% k6 = param(10);
% p1 = param(11);
% p2 = param(12);
% s1 = param(13);
% s2 = param(14);
% s3 = param(15);
% s4 = param(16);

if 1
    a = X./Z;
    b = Y./Z;
    r = sqrt(a.^2 + b.^2);
    th = atan(r);
else
    a = X;
    b = Y;
    r = sqrt(X.^2 + Y.^2);
    th = atan2(r, Z);   
end


th2 = th.^2;
th4 = th2.^2;
th6 = th2.*th4;
th8 = th4.*th4;
th10 = th6.*th4;
th12 = th6.*th6;

thd = th.*(1 + k1.*th2 + k2.*th4 + k3.*th6 + k4.*th8 + k5.*th10 + k6.*th12);
x_r = a./r.*thd;
y_r = b./r.*thd;
r_d = sqrt(x_r.^2 + y_r.^2);
uvDistorted = [x_r      +      p1.*(2.*x_r.^2 + r_d.^2) + 2.*x_r.*y_r.*p2     +     s1.*r_d.^2 + s2.*r_d.^4   + s5.*r_d.^6 ...
               y_r      +      p2.*(2.*y_r.^2 + r_d.^2) + 2.*x_r.*y_r.*p1     +     s3.*r_d.^2 + s4.*r_d.^4   + s6.*r_d.^6];

           
[matTilt, dMatTiltdTauX, dMatTiltdTauY, invMatTilt, dInvMatTiltdTauX, dInvMatTiltdTauY] = computeTiltProjectionMatrix(t1, t2);
uvTilted =    (matTilt * [uvDistorted ones(size(uvDistorted,1),1)]')';
uvTilted = uvTilted./repmat(uvTilted(:,3),1,3);
      

% pt2d = [fx.*uvDistorted(:,1) + cx ...
%     fy.*uvDistorted(:,2) + cy];
pt2d = [fx.*uvTilted(:,1) + cx ...
        fy.*uvTilted(:,2) + cy];
 



duvDistorted_dxryr = compute_duvDistorted_dxryr([x_r y_r], x_r.^2 + y_r.^2, param);


d_thd_d_th = 1 + 3*k1.*th2 + 5*k2.*th4 + 7*k3.*th6 + 9*k4.*th8 + 11*k5.*th10 + 13*k6.*th12;

d_x_r_d_a = b.^2./r.^3.*thd + a.^2./(r.^2 + r.^4)*d_thd_d_th;
d_x_r_d_b = -a.*b./r.^3.*thd + a.*b./(r.^2 + r.^4)*d_thd_d_th;
d_y_r_d_a = -a.*b./r.^3.*thd + a.*b./(r.^2 + r.^4)*d_thd_d_th;
d_y_r_d_b = a.^2./r.^3.*thd + b.^2./(r.^2 + r.^4)*d_thd_d_th;
d_xr_yr_d_ab = [d_x_r_d_a d_x_r_d_b; d_y_r_d_a d_y_r_d_b];

d_ab_d_xyz = [ 1./Z    0     -X./Z.^2;
                 0    1./Z   -Y./Z.^2];
         
             
             

stx = sin(t1);
ctx = cos(t1);
sty = sin(t2);
cty = cos(t2);

tt1 = ctx;
tt4 = -stx*sty;
tt5 = cty;
tt7 = sty;
tt8 = -stx*cty;
tt9 = ctx*cty;
             
u_d = uvDistorted(:,1);
v_d = uvDistorted(:,2);
 
d_ut_d_ud = tt1/(tt7*u_d + tt8*v_d + tt9) - tt1*u_d*tt7/(tt7*u_d + tt8*v_d + tt9)^2;

d_ut_d_vd = -tt1*u_d*tt8/(tt7*u_d + tt8*v_d + tt9)^2;

d_vt_d_ud = tt4/(tt7*u_d + tt8*v_d + tt9) - (tt4*u_d+tt5*v_d)*tt7/(tt7*u_d + tt8*v_d + tt9)^2;

d_vt_d_vd = tt5/(tt7*u_d + tt8*v_d + tt9) - (tt4*u_d+tt5*v_d)*tt8/(tt7*u_d + tt8*v_d + tt9)^2;
             
d_uvTilted_d_uvDistorted = [d_ut_d_ud d_ut_d_vd; d_vt_d_ud d_vt_d_vd];
d_uv_d_uvTilted = [fx 0;0 fy];

d_uv_d_pt3d = d_uv_d_uvTilted * d_uvTilted_d_uvDistorted * duvDistorted_dxryr * d_xr_yr_d_ab * d_ab_d_xyz;




%%%%
d_uv_d_fxfycxcy = [uvTilted(:,1)     0            1   0;
                         0         uvTilted(:,2)  0   1];

d_xr_yr_d_thd = [a./r; b./r];
                     
d_thd_d_k1k2k3k4k5k6 = [th.^3   th.^5   th.^7   th.^9   th.^11   th.^13];

d_uv_d_k1k2k3k4k5k6 = d_uv_d_uvTilted * d_uvTilted_d_uvDistorted * duvDistorted_dxryr * d_xr_yr_d_thd * d_thd_d_k1k2k3k4k5k6;

duvDistorted_d_p1p2 = [(2.*x_r.^2 + r_d.^2)        2.*x_r.*y_r;
                          2.*x_r.*y_r         (2.*y_r.^2 + r_d.^2)];

duvDistorted_d_s1s2s3s4s5s6 = [r_d.^2   r_d.^4    0        0     r_d.^6     0 ;
                             0        0     r_d.^2      r_d.^4     0      r_d.^6];
                         
d_uv_d_p1p2 = d_uv_d_uvTilted * d_uvTilted_d_uvDistorted  * duvDistorted_d_p1p2;

d_uv_d_s1s2s3s4s5s6 = d_uv_d_uvTilted * d_uvTilted_d_uvDistorted  * duvDistorted_d_s1s2s3s4s5s6;


temp = sty*u_d-stx*cty*v_d+ctx*cty;

d_ut_d_tx = -stx*u_d/temp - ctx*u_d*(-ctx*cty*v_d-stx*cty)/temp^2;

d_ut_d_ty = -ctx*u_d*(cty*u_d+stx*sty*v_d-ctx*sty)/temp^2;

d_vt_d_tx = (-ctx*sty*u_d+cty*v_d)/temp - (cty*v_d-stx*sty*u_d)*(-ctx*cty*v_d-stx*cty)/temp^2;

d_vt_d_ty = (-stx*cty*u_d-sty*v_d)/temp-(cty*v_d-stx*sty*u_d)*(cty*u_d+stx*sty*v_d-ctx*sty)/temp^2;

d_uvTilted_d_t1t2 = [d_ut_d_tx d_ut_d_ty; d_vt_d_tx d_vt_d_ty];

d_uv_d_t1t2 = d_uv_d_uvTilted * d_uvTilted_d_t1t2;

d_uv_d_param = [d_uv_d_fxfycxcy d_uv_d_k1k2k3k4k5k6 d_uv_d_p1p2 d_uv_d_s1s2s3s4s5s6 d_uv_d_t1t2];





end
function [matTilt, dMatTiltdTauX, dMatTiltdTauY, invMatTilt, dInvMatTiltdTauX, dInvMatTiltdTauY] =  computeTiltProjectionMatrix(tauX, tauY)

cTauX = cos(tauX);
sTauX = sin(tauX);
cTauY = cos(tauY);
sTauY = sin(tauY);
matRotX = [1,0,0;0,cTauX,sTauX;0,-sTauX,cTauX];
matRotY = [cTauY,0,-sTauY;0,1,0;sTauY,0,cTauY];
matRotXY = matRotY * matRotX;
matProjZ = [matRotXY(3,3),0,-matRotXY(1,3);0,matRotXY(3,3),-matRotXY(2,3);0,0,1];

% Matrix for trapezoidal distortion of tilted image sensor
matTilt = matProjZ * matRotXY;

% Derivative with respect to tauX
dMatRotXYdTauX = matRotY * [0,0,0;0,-sTauX,cTauX;0,-cTauX,-sTauX];
dMatProjZdTauX = [dMatRotXYdTauX(3,3),0,-dMatRotXYdTauX(1,3);0,dMatRotXYdTauX(3,3),-dMatRotXYdTauX(2,3);0,0,0];
dMatTiltdTauX = (matProjZ * dMatRotXYdTauX) + (dMatProjZdTauX * matRotXY);
dMatTiltdTauX_check = [-sin(tauX) 0 0; -cos(tauX)*sin(tauY) 0  0; 0 -cos(tauX)*cos(tauY) -sin(tauX)*cos(tauY)];
err1 = dMatTiltdTauX_check - dMatTiltdTauX;

% Derivative with respect to tauY
dMatRotXYdTauY = [-sTauY,0,-cTauY;0,0,0;cTauY,0,-sTauY] * matRotX;
dMatProjZdTauY = [dMatRotXYdTauY(3,3),0,-dMatRotXYdTauY(1,3);0,dMatRotXYdTauY(3,3),-dMatRotXYdTauY(2,3);0,0,0];
dMatTiltdTauY = (matProjZ * dMatRotXYdTauY) + (dMatProjZdTauY * matRotXY);
dMatTiltdTauY_check = [0 0 0; -sin(tauX)*cos(tauY) -sin(tauY) 0; cos(tauY) sin(tauX)*sin(tauY) -cos(tauX)*sin(tauY)];
err2 = dMatTiltdTauY_check - dMatTiltdTauY;

invZ = 1./matRotXY(3,3);
invMatProjZ = [invZ,0,invZ*matRotXY(1,3);0,invZ,invZ*matRotXY(2,3);0,0,1];
invMatTilt = matRotXY'*invMatProjZ;
invMatTilt_check = [1/cos(tauX) 0 0; tan(tauX)*tan(tauY) 1/cos(tauY) 0; -tan(tauY) tan(tauX)/cos(tauY) 1/(cos(tauX)*cos(tauY))];
err3 = invMatTilt_check - invMatTilt;

dInvMatTiltdTauX = [sin(tauX)/(cos(tauX) * cos(tauX))               0                                     0;
                   tan(tauY)/(cos(tauX) * cos(tauX))                0                                     0;
                                   0                   1/cos(tauX)/cos(tauX)/cos(tauY)  sin(tauX)/cos(tauX)/cos(tauX)/cos(tauY)];

dInvMatTiltdTauY = [            0                                   0                                        0;
                  tan(tauX)/cos(tauY)/cos(tauY)       sin(tauY)/cos(tauY)/cos(tauY)                          0;
                     -1/cos(tauY)/cos(tauY)     tan(tauX)*sin(tauY)/cos(tauY)/cos(tauY)    sin(tauY)/cos(tauX)/cos(tauY)/cos(tauY)];
end
function duvDistorted_dxryr = compute_duvDistorted_dxryr(xr_yr, xr_yr_squaredNorm, param)
fx = param(1);
fy = param(2);
cx = param(3);
cy = param(4);
k1 = param(5);
k2 = param(6);
k3 = param(7);
k4 = param(8);
k5 = param(9);
k6 = param(10);
p1 = param(11);
p2 = param(12);
s1 = param(13);
s2 = param(14);
s3 = param(15);
s4 = param(16);
s5 = param(17);
s6 = param(18);
t1 = param(19);
t2 = param(20);
duvDistorted_dxryr(1,1,:) = 1 + 6 .* xr_yr(:,1).*p1 + 2.*xr_yr(:,2).*p2   +   2.*(s1 + 2.*s2.*xr_yr_squaredNorm + 3.*s5.*xr_yr_squaredNorm.^2 ).*xr_yr(:,1);
duvDistorted_dxryr(1,2,:) = 2.*xr_yr(:,2).*p1 + 2.*xr_yr(:,1).*p2         +   2.*(s1 + 2.*s2.*xr_yr_squaredNorm + 3.*s5.*xr_yr_squaredNorm.^2 ).*xr_yr(:,2);
duvDistorted_dxryr(2,1,:) = 2.*xr_yr(:,2).*p1 + 2.*xr_yr(:,1).*p2         +   2.*(s3 + 2.*s4.*xr_yr_squaredNorm + 3.*s6.*xr_yr_squaredNorm.^2 ).*xr_yr(:,1);
duvDistorted_dxryr(2,2,:) = 1 + 6 .* xr_yr(:,2).*p2 + 2.*xr_yr(:,1).*p1   +   2.*(s3 + 2.*s4.*xr_yr_squaredNorm + 3.*s6.*xr_yr_squaredNorm.^2 ).*xr_yr(:,2);



end
function [th, dthD_dth] = getThetaFromNorm_xr_yr(th_radialDesired, param)
test_Jac = 0;


% fx = Param(1);
% fy = Param(2);
% cx = Param(3);
% cy = Param(4);
% k1 = Param(5);
% k2 = Param(6);
% k3 = Param(7);
% k4 = Param(8);
% k5 = Param(9);
% k6 = Param(10);
% p1 = Param(11);
% p2 = Param(12);
% s1 = Param(13);
% s2 = Param(14);
% s3 = Param(15);
% s4 = Param(16);
% t1 = param(17);
% t2 = param(18);

% param = Param;


%  initial guess
th = th_radialDesired;



startK = 5;
max_iter = 20;
for i = 1 : max_iter
    
    if 1
        thetaSq = th * th;
        th_radial = 1;
        dthD_dth = 1;
        theta2is = thetaSq;
        for j = 0:5
            th_radial = th_radial + theta2is * param(startK + j);
            dthD_dth = dthD_dth + (2 * j + 3) * param(startK + j) * theta2is;
            theta2is = theta2is * thetaSq;
        end
    else
        th2 = th.^2;
        th4 = th2.^2;
        th6 = th2.*th4;
        th8 = th4.*th4;
        th10 = th6.*th4;
        th12 = th6.*th6;
        dthD_dth = 1 + 3*k1.*th2 + 5*k2.*th4 + 7*k3.*th6 + 9*k4.*th8 + 11*k5.*th10 + 13*k6.*th12;
    end
    
    th_radial = th_radial * th;
%     if(~test_Jac)
        if (abs(dthD_dth) > 1e-20)
            step = (th_radialDesired - th_radial) / dthD_dth;
        else
            
            if (th_radialDesired - th_radial) * dthD_dth > 0.0
                step =  1e-19;
            else
                step = -1e-19;
            end
        end
        
        th = th + step;
        if(norm(step) < 1e-20)
            break;
        end
        
        if (abs(th) >=3.1415926 / 2.0)
            
            th = (0.999) * 3.1415926 / 2.0;
        end
%     else
%          step = (th_radialDesired - th_radial) / dthD_dth;
%           th = th + step;
%     end
end

end
function [xr_yr, duvDistorted_dxryr] = compute_xr_yr_from_uvDistorted(uvDistorted, param)
test_Jac = 0;
fx = param(1);
fy = param(2);
cx = param(3);
cy = param(4);
k1 = param(5);
k2 = param(6);
k3 = param(7);
k4 = param(8);
k5 = param(9);
k6 = param(10);
p1 = param(11);
p2 = param(12);
s1 = param(13);
s2 = param(14);
s3 = param(15);
s4 = param(16);
s5 = param(17);
s6 = param(18);
t1 = param(19);
t2 = param(20);
%initial guess:
xr_yr = uvDistorted;
max_iter = 20;
for i = 1 : max_iter
    uvDistorted_est = xr_yr;
    xr_yr_squaredNorm = xr_yr(:,1).^2 + xr_yr(:,2).^2;
    x_r = xr_yr(:,1);
    y_r = xr_yr(:,2);
    r_d = sqrt(x_r.^2 + y_r.^2);
    uvDistorted_est = uvDistorted_est + [ p1.*(2.*x_r.^2 + r_d.^2) + 2.*x_r.*y_r.*p2     +     s1.*r_d.^2 + s2.*r_d.^4 + s5.*r_d.^6 ...
                                          p2.*(2.*y_r.^2 + r_d.^2) + 2.*x_r.*y_r.*p1     +     s3.*r_d.^2 + s4.*r_d.^4 + s6.*r_d.^6];
    duvDistorted_dxryr = compute_duvDistorted_dxryr(xr_yr, xr_yr_squaredNorm, param);
    correction = (inv(duvDistorted_dxryr) * (uvDistorted' - uvDistorted_est'))';
    xr_yr  = xr_yr + correction;
    if(~test_Jac)
        if(norm(correction) < 1e-20)
            break;
        end
    end
end



end