function [pt2d, d_uv_d_pt3d, d_uv_d_param] = projectFisheyeRadTanThinPrism(XYZ, param,  R, t)

XYZ(:,4) = 1;
XYZ = ([R t; 0 0 0 1] * XYZ')';
XYZ = XYZ(:,1:3);


d_uv_d_pt3d = []; d_uv_d_param = [];

X = XYZ(:,1);
Y = XYZ(:,2);
Z = XYZ(:,3);

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


if (size(XYZ,1) == 1)
    
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

end