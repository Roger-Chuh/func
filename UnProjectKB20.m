function [pt3d, d_pt3d_d_uv, d_pt3d_d_param] = UnProjectKB20(u,v, fx,fy,cx,cy,k1,k2,k3,k4,k5,k6,p1,p2,s1,s2,s3,s4,s5,s6,t1,t2)

d_pt3d_d_uv = []; d_pt3d_d_param = [];
% u = pt2d(:,1);
% v = pt2d(:,2);
% Z = pt3d(:,3);

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

uvTilted = [(u - cx)./fx (v - cy)./fy];

[matTilt, dMatTiltdTauX, dMatTiltdTauY, invMatTilt, dInvMatTiltdTauX, dInvMatTiltdTauY] = computeTiltProjectionMatrix(t1, t2);

uvDistorted =    (invMatTilt * [uvTilted ones(size(uvTilted,1),1)]')';
uvDistorted = uvDistorted./repmat(uvDistorted(:,3),1,3);

uvDistorted = uvDistorted(:,1:2);

if 1 %isempty(data.xr_yr)
    [xr_yr, duvDistorted_dxryr] = compute_xr_yr_from_uvDistorted(uvDistorted, param);
    data.xr_yr = xr_yr;
    data.duvDistorted_dxryr = duvDistorted_dxryr;
else
    xr_yr = data.xr_yr;
    duvDistorted_dxryr = data.duvDistorted_dxryr;
end
xr_yrNorm = norm(xr_yr);
if (xr_yrNorm == 0.0)
    pt3d = [0 0 1];
else
    [theta, dthD_dth] = getThetaFromNorm_xr_yr(xr_yrNorm, param);
end

if 0
    pt3d = [tan(theta) / xr_yrNorm * xr_yr 1];
else
    thd = xr_yrNorm;
    scaling = sin(theta)./thd;
    x = xr_yr(:,1)*scaling;
    y = xr_yr(:,2)*scaling;
    z = cos(theta);
    pt3d = [x y z];
    
    sin_theta = sin(theta);
    cos_theta = cos(theta);
    mx = xr_yr(:,1);
    my = xr_yr(:,2);
    
    x_r = mx;
    y_r = my;
    r_d = sqrt(x_r.^2 + y_r.^2);
    
    d_thetad_d_mx = mx / thd;
    d_thetad_d_my = my / thd;
    
end


theta2 = theta * theta;
d_scaling_d_thetad = (thd * cos_theta / dthD_dth - sin_theta) / (thd * thd);
d_cos_d_thetad = -sin_theta / dthD_dth;
d_scaling_d_k1 = -cos_theta * theta * theta2 / (dthD_dth * thd);
d_cos_d_k1 =  -d_cos_d_thetad * theta * theta2;


d_pt3d_d_k1 = [mx * d_scaling_d_k1;
               my * d_scaling_d_k1;
               d_cos_d_k1];
d_pt3d_d_k2 = d_pt3d_d_k1 .* theta2;
d_pt3d_d_k3 = d_pt3d_d_k2 .* theta2;
d_pt3d_d_k4 = d_pt3d_d_k3 .* theta2;
d_pt3d_d_k5 = d_pt3d_d_k4 .* theta2;
d_pt3d_d_k6 = d_pt3d_d_k5 .* theta2;

d_X_d_mx = scaling + mx * d_scaling_d_thetad * d_thetad_d_mx;
d_X_d_my = mx * d_scaling_d_thetad * d_thetad_d_my;

d_Y_d_mx = my * d_scaling_d_thetad * d_thetad_d_mx;
d_Y_d_my = scaling + my * d_scaling_d_thetad * d_thetad_d_my;

d_Z_d_mx = d_cos_d_thetad * d_thetad_d_mx;
d_Z_d_my = d_cos_d_thetad * d_thetad_d_my;


d_pt3d_d_mxmy = [d_X_d_mx d_X_d_my;
                d_Y_d_mx d_Y_d_my;
                d_Z_d_mx d_Z_d_my];

d_pt3d_d_xryr = d_pt3d_d_mxmy;

d_xryr_duvDistorted = inv(duvDistorted_dxryr);
if 0
    d_tilt_d_fxfycxcy = [-uvDistorted(:,1)./fx           0             -1./fx     0;
                                  0          -uvDistorted(:,2)./fy    0      -1./fy];
else
    if 0
        d_tilt_d_fxfycxcy = [-(uvDistorted(:,1)-cx)./fx./fx           0             -1./fx     0;
            0          -(uvDistorted(:,2)-cy)./fy./fy    0      -1./fy];
        
    else
        d_tilt_d_fxfycxcy = [-(u-cx)./fx./fx           0             -1./fx     0;
                                     0          -(v-cy)./fy./fy    0      -1./fy];
    end
end
d_tilt_d_uv = [1./fx      0;
               0        1./fy];

           
           
           
           
stx = sin(t1);
ctx = cos(t1);
sty = sin(t2);
cty = cos(t2);
ttx = tan(t1);
tty = tan(t2);

tt1 = 1/ctx;
tt4 = ttx*tty;
tt5 = 1/cty;
tt7 = -tty;
tt8 = ttx/cty;
tt9 = 1/ctx/cty;
             
u_t = uvTilted(:,1);
v_t = uvTilted(:,2);

d_ud_d_ut = tt1/(tt7*u_t + tt8*v_t + tt9) - tt1*u_t*tt7/(tt7*u_t + tt8*v_t + tt9)^2;

d_ud_d_vt = -tt1*u_t*tt8/(tt7*u_t + tt8*v_t + tt9)^2;

d_vd_d_ut = tt4/(tt7*u_t + tt8*v_t + tt9) - (tt4*u_t+tt5*v_t)*tt7/(tt7*u_t + tt8*v_t + tt9)^2;

d_vd_d_vt = tt5/(tt7*u_t + tt8*v_t + tt9) - (tt4*u_t+tt5*v_t)*tt8/(tt7*u_t + tt8*v_t + tt9)^2;
             
d_uvDistorted_d_uvTilted = [d_ud_d_ut d_ud_d_vt; d_vd_d_ut d_vd_d_vt];
           
           
           
           
           
d_udvd_d_tilt = d_uvDistorted_d_uvTilted;




temp = -tty*u_t + ttx/cty*v_t+1/ctx/cty;           


d_ud_d_tx = (stx/ctx^2*u_t/temp)-(1/ctx*u_t*(v_t/ctx^2/cty+cty*stx/(ctx*cty)^2)/temp^2);

d_ud_d_ty = -u_t/ctx*(-u_t/cty^2+ttx*sty*v_t/cty^2+ctx*sty/(ctx*cty)^2)/temp^2;

d_vd_d_tx = (1/ctx^2*tty*u_t/temp) - ((ttx*tty*u_t+1/cty*v_t)*(1/cty/ctx^2*v_t+cty*stx/(ctx*cty)^2)/temp^2);

d_vd_d_ty = ((ttx*u_t/cty^2 + sty*v_t/cty^2)/temp) - ((ttx*tty*u_t+1/cty*v_t)*(-u_t/cty^2+ttx*sty*v_t/cty^2+ctx*sty/(ctx*cty)^2)/temp^2);

d_udvd_d_t1t2 = [d_ud_d_tx d_ud_d_ty; d_vd_d_tx d_vd_d_ty];
%%
duvDistorted_d_p1p2 = -[(2.*x_r.^2 + r_d.^2)        2.*x_r.*y_r;
                          2.*x_r.*y_r         (2.*y_r.^2 + r_d.^2)]; 
%%
duvDistorted_d_s1s2s3s4s5s6 = -[r_d.^2   r_d.^4    0        0     r_d.^6    0;
                                0        0     r_d.^2   r_d.^4      0     r_d.^6];           
                         
                         
d_xryr_d_uv = d_xryr_duvDistorted * d_udvd_d_tilt * d_tilt_d_uv;
d_xryr_d_fxfycxcy = d_xryr_duvDistorted * d_udvd_d_tilt* d_tilt_d_fxfycxcy;
% d_xryr_d_t1t2 = d_xryr_duvDistorted * d_udvd_d_tilt* d_tilt_d_t1t2;
d_xryr_d_t1t2 = d_xryr_duvDistorted * d_udvd_d_t1t2;
d_xryr_d_p1p2 = d_xryr_duvDistorted * duvDistorted_d_p1p2;
d_xryr_d_s1s2s3s4s5s6 = d_xryr_duvDistorted * duvDistorted_d_s1s2s3s4s5s6;


d_pt3d_d_param = [d_pt3d_d_mxmy * d_xryr_d_fxfycxcy ...
                   [d_pt3d_d_k1 d_pt3d_d_k2 d_pt3d_d_k3 d_pt3d_d_k4 d_pt3d_d_k5 d_pt3d_d_k6] ...
                 d_pt3d_d_mxmy * [d_xryr_d_p1p2 d_xryr_d_s1s2s3s4s5s6 d_xryr_d_t1t2]];
             
d_pt3d_d_uv = d_pt3d_d_mxmy * d_xryr_d_uv;

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
% 
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