function undistirtImgs()


inputDir = 'G:\matlab\data\vib\rog463\Camera0\images';
outputDir = 'G:\matlab\data\vib\rog463\Camera0\undistort';
Param0_bak = [2.3368958079766023e+02 2.3380527733852011e+02 3.2031702730452713e+02 2.4308165803544867e+02 2.2758418025869792e-01 -1.9867565230558476e-01 1.0540751138748485e-01 -4.1615767245009774e-02 ...
    1.0647574394068167e-02 -1.2570315612344181e-03 2.4661233870034425e-03 4.7780930589446866e-03 -2.5535605305932623e-03 -2.5761049011873316e-04 -6.7601975948881647e-03 -2.6929525666274702e-05 ]';
Param0 = [233.5640589, 233.6785296 320.3490101 243.1172132  0.2324045045 -0.2160760084 0.1320611134 -0.06176198181 0.01798935924 -0.002277538976 ...
      0.002621474503 0.004857002128 -0.002787615464 -0.0002740460884 -0.007038871399 0.0001033324432]';
ParamKB8 = [2.3430601550262170e+02 2.3463508194410093e+02 3.1996092464657170e+02 2.4071656171198191e+02 2.2065917471151195e-01 -1.8019308429323991e-01 ...
            7.2749867808954344e-02 -1.3554955559408860e-02 0 0 0 0 0 0 0 0]';
% Param00 = Param0;
% Param00(3) = Param00(3)-1;
% Param00(4) = Param00(4)-1;

fx0 = Param0(1);
fy0 = Param0(2);
cx0 = Param0(3);
cy0 = Param0(4);
k10 = Param0(5);
k20 = Param0(6);
k30 = Param0(7);
k40 = Param0(8);
k50 = Param0(9);
k60 = Param0(10);
p10 = Param0(11);
p20 = Param0(12);
s10 = Param0(13);
s20 = Param0(14);
s30 = Param0(15);
s40 = Param0(16);


X = [-10.1234 -20.5678 3];
X_ = X;
X = [X(1)/X(3) X(2)/X(3) 1];
X = X./norm(X);
[pt2d, d_uv_d_pt3d, d_uv_d_param] = project(X_(1),X_(2), X_(3),fx0,fy0,cx0,cy0,k10,k20,k30,k40,k50,k60,p10,p20,s10,s20,s30,s40);

[pt3d, d_pt3d_d_uv, d_pt3d_d_param] = unproject(pt2d(1),pt2d(2), fx0,fy0,cx0,cy0,k10,k20,k30,k40,k50,k60,p10,p20,s10,s20,s30,s40);

err_ = X - pt3d;

[xMat, yMat] = meshgrid(1:640, 1:480);
pix = [xMat(:) yMat(:)];

[xMat_, yMat_] = meshgrid(0:639, 0:479);
pix_ = [xMat_(:) yMat_(:)];


K_undist = [fx0 0 cx0; 0 fy0 cy0; 0 0 1];

pixDist_ = remapRectFisheyeRadTanThinPrism(pix_', K_undist, Param0, eye(3));
pixDist = pixDist_+1;

pixDist_KB8 = remapRectFisheyeRadTanThinPrism(pix_', K_undist, ParamKB8, eye(3));
pixDistKB8 = pixDist_KB8+1;

% pixDist0 = remapRectFisheyeRadTanThinPrism(pix', K_undist, Param00, eye(3));
dirInfo = dir(fullfile(inputDir, '*.bmp'));


for i = 1 : length(dirInfo)
    
   img_fish = imread(fullfile(inputDir, dirInfo(i).name)); 
   img_fish_undist = uint8(interp2(xMat,yMat,double(img_fish),reshape(pixDist(:,1),480,640),reshape(pixDist(:,2),480,640)));
   if 0
       img_fish_undist_kb8 = uint8(interp2(xMat,yMat,double(img_fish),reshape(pixDistKB8(:,1),480,640),reshape(pixDistKB8(:,2),480,640)));
       figure,imshowpair(img_fish_undist_kb8, img_fish_undist)
   end
%    img_fish_undist0 = uint8(interp2(xMat,yMat,double(img_fish),reshape(pixDist0(:,1),480,640),reshape(pixDist0(:,2),480,640)));
   
%    figure(1),imshow([img_fish img_fish_undist], []);drawnow;
   imwrite(img_fish_undist,fullfile(outputDir, dirInfo(i).name));
end




end
function [pt2d, d_uv_d_pt3d, d_uv_d_param] = project(X,Y,Z, fx,fy,cx,cy,k1,k2,k3,k4,k5,k6,p1,p2,s1,s2,s3,s4)

d_uv_d_pt3d = []; d_uv_d_param = [];

% X = pt3d(:,1);
% Y = pt3d(:,2);
% Z = pt3d(:,3);

param = [fx,fy,cx,cy,k1,k2,k3,k4,k5,k6,p1,p2,s1,s2,s3,s4];
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
uvDistorted = [x_r      +      p1.*(2.*x_r.^2 + r_d.^2) + 2.*x_r.*y_r.*p2     +     s1.*r_d.^2 + s2.*r_d.^4 ...
               y_r      +      p2.*(2.*y_r.^2 + r_d.^2) + 2.*x_r.*y_r.*p1     +     s3.*r_d.^2 + s4.*r_d.^4];


pt2d = [fx.*uvDistorted(:,1) + cx ...
    fy.*uvDistorted(:,2) + cy];




duvDistorted_dxryr = compute_duvDistorted_dxryr([x_r y_r], x_r.^2 + y_r.^2, param);


d_thd_d_th = 1 + 3*k1.*th2 + 5*k2.*th4 + 7*k3.*th6 + 9*k4.*th8 + 11*k5.*th10 + 13*k6.*th12;

d_x_r_d_a = b.^2./r.^3.*thd + a.^2./(r.^2 + r.^4)*d_thd_d_th;
d_x_r_d_b = -a.*b./r.^3.*thd + a.*b./(r.^2 + r.^4)*d_thd_d_th;
d_y_r_d_a = -a.*b./r.^3.*thd + a.*b./(r.^2 + r.^4)*d_thd_d_th;
d_y_r_d_b = a.^2./r.^3.*thd + b.^2./(r.^2 + r.^4)*d_thd_d_th;
d_xr_yr_d_ab = [d_x_r_d_a d_x_r_d_b; d_y_r_d_a d_y_r_d_b];

d_ab_d_xyz = [ 1./Z    0     -X./Z.^2;
                 0    1./Z   -Y./Z.^2];
             
d_uv_d_uvDistorted = [fx 0;0 fy];

d_uv_d_pt3d = d_uv_d_uvDistorted * duvDistorted_dxryr * d_xr_yr_d_ab * d_ab_d_xyz;




%%%%
d_uv_d_fxfycxcy = [uvDistorted(:,1)     0            1   0;
                         0         uvDistorted(:,2)  0   1];

d_xr_yr_d_thd = [a./r; b./r];
                     
d_thd_d_k1k2k3k4k5k6 = [th.^3   th.^5   th.^7   th.^9   th.^11   th.^13];

d_uv_d_k1k2k3k4k5k6 = d_uv_d_uvDistorted * duvDistorted_dxryr * d_xr_yr_d_thd * d_thd_d_k1k2k3k4k5k6;

duvDistorted_d_p1p2 = [(2.*x_r.^2 + r_d.^2)        2.*x_r.*y_r;
                          2.*x_r.*y_r         (2.*y_r.^2 + r_d.^2)];

duvDistorted_d_s1s2s3s4 = [r_d.^2   r_d.^4    0        0;
                             0        0     r_d.^2   r_d.^4];
                         
d_uv_d_p1p2 = d_uv_d_uvDistorted * duvDistorted_d_p1p2;

d_uv_d_s1s2s3s4 = d_uv_d_uvDistorted * duvDistorted_d_s1s2s3s4;

d_uv_d_param = [d_uv_d_fxfycxcy d_uv_d_k1k2k3k4k5k6 d_uv_d_p1p2 d_uv_d_s1s2s3s4];



end
function [pt3d, d_pt3d_d_uv, d_pt3d_d_param] = unproject(u,v, fx,fy,cx,cy,k1,k2,k3,k4,k5,k6,p1,p2,s1,s2,s3,s4)


d_pt3d_d_uv = []; d_pt3d_d_param = [];
% u = pt2d(:,1);
% v = pt2d(:,2);
% Z = pt3d(:,3);

param = [fx,fy,cx,cy,k1,k2,k3,k4,k5,k6,p1,p2,s1,s2,s3,s4];
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

uvDistorted = [(u - cx)./fx (v - cy)./fy];

if 1 %isempty(data.xr_yr)
    [xr_yr, duvDistorted_dxryr] = compute_xr_yr_from_uvDistorted(uvDistorted, param);

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
d_udvd_d_fxfycxcy = [-uvDistorted(:,1)./fx           0             -1./fx     0;
                                 0          -uvDistorted(:,2)./fy    0      -1./fy];
d_udvd_d_uv = [1./fx      0;
               0        1./fy];

%%
duvDistorted_d_p1p2 = -[(2.*x_r.^2 + r_d.^2)        2.*x_r.*y_r;
                          2.*x_r.*y_r         (2.*y_r.^2 + r_d.^2)]; 
%%
duvDistorted_d_s1s2s3s4 = -[r_d.^2   r_d.^4    0        0;
                             0        0     r_d.^2   r_d.^4];           
                         
                         
d_xryr_d_uv = d_xryr_duvDistorted * d_udvd_d_uv;
d_xryr_d_fxfycxcy = d_xryr_duvDistorted * d_udvd_d_fxfycxcy;
d_xryr_d_p1p2 = d_xryr_duvDistorted * duvDistorted_d_p1p2;
d_xryr_d_s1s2s3s4 = d_xryr_duvDistorted * duvDistorted_d_s1s2s3s4;


d_pt3d_d_param = [d_pt3d_d_mxmy * d_xryr_d_fxfycxcy ...
                   [d_pt3d_d_k1 d_pt3d_d_k2 d_pt3d_d_k3 d_pt3d_d_k4 d_pt3d_d_k5 d_pt3d_d_k6] ...
                 d_pt3d_d_mxmy * [d_xryr_d_p1p2 d_xryr_d_s1s2s3s4]];
             
d_pt3d_d_uv = d_pt3d_d_mxmy * d_xryr_d_uv;

end
function [xr_yr, duvDistorted_dxryr] = compute_xr_yr_from_uvDistorted(uvDistorted, param)

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

%initial guess:
xr_yr = uvDistorted;
max_iter = 20;
for i = 1 : max_iter
    uvDistorted_est = xr_yr;
    xr_yr_squaredNorm = xr_yr(:,1).^2 + xr_yr(:,2).^2;
    x_r = xr_yr(:,1);
    y_r = xr_yr(:,2);
    r_d = sqrt(x_r.^2 + y_r.^2);
    uvDistorted_est = uvDistorted_est + [ p1.*(2.*x_r.^2 + r_d.^2) + 2.*x_r.*y_r.*p2     +     s1.*r_d.^2 + s2.*r_d.^4 ...
                                          p2.*(2.*y_r.^2 + r_d.^2) + 2.*x_r.*y_r.*p1     +     s3.*r_d.^2 + s4.*r_d.^4];
    duvDistorted_dxryr = compute_duvDistorted_dxryr(xr_yr, xr_yr_squaredNorm, param);
    correction = (inv(duvDistorted_dxryr) * (uvDistorted' - uvDistorted_est'))';
    xr_yr  = xr_yr + correction;
    
    if(norm(correction) < 1e-20)
        break;
    end
    
end



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

duvDistorted_dxryr(1,1,:) = 1 + 6 .* xr_yr(:,1).*p1 + 2.*xr_yr(:,2).*p2   +   2.*(s1 + 2.*s2.*xr_yr_squaredNorm).*xr_yr(:,1);
duvDistorted_dxryr(1,2,:) = 2.*xr_yr(:,2).*p1 + 2.*xr_yr(:,1).*p2         +   2.*(s1 + 2.*s2.*xr_yr_squaredNorm).*xr_yr(:,2);
duvDistorted_dxryr(2,1,:) = 2.*xr_yr(:,2).*p1 + 2.*xr_yr(:,1).*p2         +   2.*(s3 + 2.*s4.*xr_yr_squaredNorm).*xr_yr(:,1);
duvDistorted_dxryr(2,2,:) = 1 + 6 .* xr_yr(:,2).*p2 + 2.*xr_yr(:,1).*p1   +   2.*(s3 + 2.*s4.*xr_yr_squaredNorm).*xr_yr(:,2);



end
function [th, dthD_dth] = getThetaFromNorm_xr_yr(th_radialDesired, param)


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
