function testFisheyeRadTanThinPrismOpenCV()
global Param test_Jac data
data.xr_yr = [];
Param = [1214.025993758243,
    1212.925993758243,
    1469.535869576514,
    1438.595858396595,
    0.4049855059199108,
    -0.4960655170834019,
    0.2366889743579386,
    0.9948191490966187,
    -1.573948083968366,
    0.6086728193625544,
    -0.0003672190538065101,
    -0.0007413905358394097,
    -0.0005839984838196491,
    -0.0002193412417918993,
    0.0001936268547567258,
    -0.0002042473404798113];

Param = [242.0503165415735,
          241.8503165415735,
          319.8053191211108,
          240.9679064324467,
          -0.02689165420917155,
          0.1019919163316408,
          -0.07180862612716554,
          0.01245898233008245,
          0.001242086400994195,
          -0.0004149914187376608,
          0.0005499786009926172,
          0.0005819661606400815,
          -0.001232007804086031,
          -0.0003868334860997226,
          -0.0005963449586486949,
          0.00001071922498768078];

% Param(9:10) = 0;
% Param(11:end) = 0;

Param = [235.5949267 235.9431433 317.3055017 237.6924467  0.05494971188 -0.02426120721 -0.06438938922 0.0666085306 -0.0250625348 0.003568605662 0.0002650098114 -0.005080557258 -0.000474610901 0.0009838698557 0.005597876514 0.0002624416604]';
Param = [233.6897149 233.8054466 320.318364 243.084033 0.2275629214 -0.1985754394 0.1052303283 -0.04146862157 0.01058943721 -0.001248198038 0.00246700086 0.004779722788 -0.002555542365 -0.0002580159852 -0.006763030277 -2.772998947e-05  ]';

Param = [234.7771905 234.7827464 330.6800546 233.8421571 0.2175176991 -0.1822643617 0.07829997996 -0.01315725164 -0.004524769258 0.001810066375 0.006967525661 -0.001222621116 -0.007796129172 -0.0001589226627 -0.000399970679 7.753565425e-05];


test_Jac = 0;


Param0 = Param;

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


fx = Param(1);
fy = Param(2);
cx = Param(3);
cy = Param(4);
k1 = Param(5);
k2 = Param(6);
k3 = Param(7);
k4 = Param(8);
k5 = Param(9);
k6 = Param(10);
p1 = Param(11);
p2 = Param(12);
s1 = Param(13);
s2 = Param(14);
s3 = Param(15);
s4 = Param(16);


X = [-10.1234 -20.5678 30];
X_ = X;
X = [X(1)/X(3) X(2)/X(3) 1];
X = X./norm(X);

[pt2d, d_uv_d_pt3d, d_uv_d_param] = project(X_(1),X_(2), X_(3),fx0,fy0,cx0,cy0,k10,k20,k30,k40,k50,k60,p10,p20,s10,s20,s30,s40);

[pt3d, d_pt3d_d_uv, d_pt3d_d_param] = unproject(pt2d(1),pt2d(2), fx0,fy0,cx0,cy0,k10,k20,k30,k40,k50,k60,p10,p20,s10,s20,s30,s40);

[pt3d11, d_pt3d_d_uv11, d_pt3d_d_param11] = unproject(100,200, fx0,fy0,cx0,cy0,k10,k20,k30,k40,k50,k60,p10,p20,s10,s20,s30,s40);
 


err_ = X - pt3d

disturbance = [1.*(rand(4,1)-0.5);
    0.0001.*(rand(6,1)-0.5);
    0.000001.*(rand(6,1)-0.5)];
% disturbance(5:end) = 0;

disturbance_pt3d = 0.01.*(rand(3,1)-0.5);
disturbance_pt2d = 1000.*(rand(2,1)-0.5);

% disturbance(11:end) = 0;

Param = Param + disturbance;
fx = Param(1);
fy = Param(2);
cx = Param(3);
cy = Param(4);
k1 = Param(5);
k2 = Param(6);
k3 = Param(7);
k4 = Param(8);
k5 = Param(9);
k6 = Param(10);
p1 = Param(11);
p2 = Param(12);
s1 = Param(13);
s2 = Param(14);
s3 = Param(15);
s4 = Param(16);
[pt2d_err, d_uv_d_pt3d_err, d_uv_d_param_err] = project(X_(1),X_(2), X_(3),fx,fy,cx,cy,k1,k2,k3,k4,k5,k6,p1,p2,s1,s2,s3,s4);
[pt2d_err1, d_uv_d_pt3d_err1, d_uv_d_param_err1] = project(X_(1)+disturbance_pt3d(1),X_(2)+disturbance_pt3d(2), X_(3)+disturbance_pt3d(3), fx0,fy0,cx0,cy0,k10,k20,k30,k40,k50,k60,p10,p20,s10,s20,s30,s40);

err1 = pt2d_err' - pt2d';
err_comp1 = d_uv_d_param *disturbance;
check_Jac1 = [err1 - err_comp1 err1  err_comp1];

err11 = pt2d_err1' - pt2d';
err_comp11 = d_uv_d_pt3d *disturbance_pt3d;
check_Jac11 = [err11 - err_comp11 err11  err_comp11];


[pt3d_err, d_pt3d_d_uv_err, d_pt3d_d_param_err] = unproject(pt2d(1),pt2d(2), fx,fy,cx,cy,k1,k2,k3,k4,k5,k6,p1,p2,s1,s2,s3,s4);
[pt3d_err2, d_pt3d_d_uv_err2, d_pt3d_d_param_err2] = unproject(pt2d(1)+disturbance_pt2d(1),pt2d(2)+disturbance_pt2d(2), fx0,fy0,cx0,cy0,k10,k20,k30,k40,k50,k60,p10,p20,s10,s20,s30,s40);


err2 = pt3d_err' - pt3d';
err_comp2 = d_pt3d_d_param*disturbance;
check_Jac2 = [err2 - err_comp2 err2  err_comp2];


err22 = pt3d_err2' - pt3d';
err_comp22 = d_pt3d_d_uv*disturbance_pt2d;
check_Jac22 = [err22 - err_comp22 err22  err_comp22];

err_cat = [check_Jac1; check_Jac11];
% err_cat(3,:) = 0;
err_cat = [err_cat; check_Jac2; check_Jac22];


[pt2d_all_err1, d_uv_d_pt3d_all_err1, d_uv_d_param_all_err1] = project(X_(1)+disturbance_pt3d(1),X_(2)+disturbance_pt3d(2), X_(3)+disturbance_pt3d(3), fx,fy,cx,cy,k1,k2,k3,k4,k5,k6,p1,p2,s1,s2,s3,s4);
err3 = pt2d_all_err1' - pt2d';
err_comp3 = [d_uv_d_pt3d, d_uv_d_param] * [disturbance_pt3d; disturbance];
check_Jac3 = [err3 - err_comp3 err3  err_comp3];

[pt3d_all_err2, d_pt3d_d_uv_all_err2, d_pt3d_d_param_all_err2] = unproject(pt2d(1)+disturbance_pt2d(1),pt2d(2)+disturbance_pt2d(2), fx,fy,cx,cy,k1,k2,k3,k4,k5,k6,p1,p2,s1,s2,s3,s4);
err33 = pt3d_all_err2' - pt3d';
err_comp33 = [ d_pt3d_d_uv d_pt3d_d_param] * [disturbance_pt2d; disturbance];
check_Jac33 = [err33 - err_comp33 err33  err_comp33];


err_cat = [err_cat; check_Jac3; check_Jac33];






[xMat, yMat] = meshgrid(1:4: 640, 1:4: 480);
pix = [xMat(:) yMat(:)];
cov_2d = zeros(2,2,size(pix,1));
cov_3d = zeros(3,3,size(pix,1));
for i = 1 : size(pix,1)
    [pt3d_i, d_pt3d_d_uv_i, d_pt3d_d_param_i] = unproject(pix(i,1),pix(i,2), fx0,fy0,cx0,cy0,k10,k20,k30,k40,k50,k60,p10,p20,s10,s20,s30,s40);
    [pt2d, d_uv_d_pt3d, d_uv_d_param] = project(pt3d_i(1),pt3d_i(2), pt3d_i(3),fx0,fy0,cx0,cy0,k10,k20,k30,k40,k50,k60,p10,p20,s10,s20,s30,s40);
    cov_2d(:,:,i) = d_pt3d_d_uv_i'*d_pt3d_d_uv_i;
    cov_3d(:,:,i) = d_uv_d_pt3d'*d_uv_d_pt3d;
end

cov_2d_11 = reshape(cov_2d(1,1,:), size(xMat));
cov_2d_12 = reshape(cov_2d(1,2,:), size(xMat));
cov_2d_21 = reshape(cov_2d(2,1,:), size(xMat));
cov_2d_22 = reshape(cov_2d(2,2,:), size(xMat));

asg = 1;

if 0
    
    syms Fx Fy Cx Cy K1 K2 K3 K4 K5 K6 P1 P2 S1 S2 S3 S4
    [pt2d, d_uv_d_pt3d, d_uv_d_param1] = project(X_(1),X_(2), X_(3),fx,fy,cx,cy,k1,K2,k3,k4,k5,k6,p1,p2,s1,s2,s3,s4);
    f(K2) = pt2d;
    d = (diff(f(K2), K2));
    diary('G:\matlab\data\vio\ci\outputs_temp\k1.txt');
    d
    diary off
    K2 = k2;
    dd1 = [ - (48275021572435495174898032358408269513219969837322525381732012136048209129685548970780841343577975855498889641513383742886785*K1)/1322111937580497197903830616065542079656809365928562438569297590548811582472622691650378420879430569695182424050046716608512 - (5400921095745737414864740567421*(((7538505412310750881198669132351350142507979705*K1)/5708990770823839524233143877797980545530986496 + 108616074649820597627439802580852225777621454955/22835963083295358096932575511191922182123945984)^2 + ((7658023570111082183785847262302949358552067749*K1)/2854495385411919762116571938898990272765493248 + 110338112698430371295717186394507058396835377199/11417981541647679048466287755595961091061972992)^2)*((291410363852946023596658517056327944386602229867452446331213491986081972889503335934635595029*K1)/16296287810675888690147565507275025288411747149327490005089123594835050398106693649467179008 + 4198690337516231856552518384262660734600025545759371363542075449882910274323304770522983108479/65185151242703554760590262029100101153646988597309960020356494379340201592426774597868716032))/10141204801825835211973625643008 - 9173334200905149215201177199746110203060254848275400734911070613524627812110541602474809405419712340332366084219678939233532355/5288447750321988791615322464262168318627237463714249754277190362195246329890490766601513683517722278780729696200186866434048, - (14570102498104968102687249988128790534735505934308818962012474138076873364262459461514613087559427073127488557231239696949549*K1)/330527984395124299475957654016385519914202341482140609642324397637202895618155672912594605219857642423795606012511679152128 - (1256175362649241410674650484283*(((7538505412310750881198669132351350142507979705*K1)/5708990770823839524233143877797980545530986496 + 108616074649820597627439802580852225777621454955/22835963083295358096932575511191922182123945984)^2 + ((7658023570111082183785847262302949358552067749*K1)/2854495385411919762116571938898990272765493248 + 110338112698430371295717186394507058396835377199/11417981541647679048466287755595961091061972992)^2)*((291410363852946023596658517056327944386602229867452446331213491986081972889503335934635595029*K1)/16296287810675888690147565507275025288411747149327490005089123594835050398106693649467179008 + 4198690337516231856552518384262660734600025545759371363542075449882910274323304770522983108479/65185151242703554760590262029100101153646988597309960020356494379340201592426774597868716032))/2535301200456458802993406410752 - 4512121468592493264397478556438384715937038910902384670412856992193050228834794986370962975451826433051050202303796805675114951/1322111937580497197903830616065542079656809365928562438569297590548811582472622691650378420879430569695182424050046716608512];
    dd2 = [ - (207960099903111898174347030742480337410707985713547487735713601464917560573278037239047498497134384296330328836845642613706785*K2)/1322111937580497197903830616065542079656809365928562438569297590548811582472622691650378420879430569695182424050046716608512 - (5400921095745737414864740567421*((1255343372322049086338669750089511338995115487397891181604146513306105782797421711189802303029*K2)/16296287810675888690147565507275025288411747149327490005089123594835050398106693649467179008 + 48740894432609441524824872929981085575585992223008904015911655044437931993929073374306174369595/260740604970814219042361048116400404614587954389239840081425977517360806369707098391474864128)*(((15646385323405317538376334653118246984114721545*K2)/5708990770823839524233143877797980545530986496 + 607498181066896240122845838517084812926392999975/91343852333181432387730302044767688728495783936)^2 + ((15894448705708353092422390404592575971634264501*K2)/2854495385411919762116571938898990272765493248 + 617129674247175268608511369635475545525414554555/45671926166590716193865151022383844364247891968)^2))/10141204801825835211973625643008 - 78457915915197855817226432569563050027678887314427276783590477663726073059890380439515460825593319044646792094011541810266282895/21153791001287955166461289857048673274508949854856999017108761448780985319561963066406054734070889115122918784800747465736192, - (62765377878869472407993739840518676831425525484341125755627919493843509009947178108998892804165690067214883676240388674697549*K2)/330527984395124299475957654016385519914202341482140609642324397637202895618155672912594605219857642423795606012511679152128 - (1256175362649241410674650484283*((1255343372322049086338669750089511338995115487397891181604146513306105782797421711189802303029*K2)/16296287810675888690147565507275025288411747149327490005089123594835050398106693649467179008 + 48740894432609441524824872929981085575585992223008904015911655044437931993929073374306174369595/260740604970814219042361048116400404614587954389239840081425977517360806369707098391474864128)*(((15646385323405317538376334653118246984114721545*K2)/5708990770823839524233143877797980545530986496 + 607498181066896240122845838517084812926392999975/91343852333181432387730302044767688728495783936)^2 + ((15894448705708353092422390404592575971634264501*K2)/2854495385411919762116571938898990272765493248 + 617129674247175268608511369635475545525414554555/45671926166590716193865151022383844364247891968)^2))/2535301200456458802993406410752 - 38154277130673311662549397021029024933420887439692192875566770858344905525858931019544992932713671553350078102082658069226036387/5288447750321988791615322464262168318627237463714249754277190362195246329890490766601513683517722278780729696200186866434048];
    d_uv_d_k2 = d_uv_d_param(:,6);
    err1 = d_uv_d_k2-dd2';
    
    sdfkb = 1;
    
elseif 0
      syms Fx Fy Cx Cy K1 K2 K3 K4 K5 K6 P1 P2 S1 S2 S3 S4
    [pt3d, d_pt3d_d_uv1, d_pt3d_d_param1] = unproject(pt2d(1),pt2d(2), fx,fy,cx,cy,K1,k2,k3,k4,k5,k6,p1,p2,s1,s2,s3,s4);
     f(K1) = pt3d;
     d = (diff(f(K1), K1));
    diary('G:\matlab\data\vio\ci\outputs_temp\k1.txt');
    d
    diary off
    K1 = k1;
end

end

function [pt2d, d_uv_d_pt3d, d_uv_d_param] = project(X,Y,Z, fx,fy,cx,cy,k1,k2,k3,k4,k5,k6,p1,p2,s1,s2,s3,s4)
global Param
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
global Param data

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
if 1
    d_udvd_d_fxfycxcy = [-uvDistorted(:,1)./fx           0             -1./fx     0;
                                 0          -uvDistorted(:,2)./fy    0      -1./fy];
else
    d_udvd_d_fxfycxcy = [-(uvDistorted(:,1)-cx)./fx./fx           0             -1./fx     0;
                                    0          -(uvDistorted(:,2)-cy)./fy./fy    0      -1./fy];
end
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
global test_Jac
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
    if(~test_Jac)
        if(norm(correction) < 1e-20)
            break;
        end
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
global test_Jac Param


fx = Param(1);
fy = Param(2);
cx = Param(3);
cy = Param(4);
k1 = Param(5);
k2 = Param(6);
k3 = Param(7);
k4 = Param(8);
k5 = Param(9);
k6 = Param(10);
p1 = Param(11);
p2 = Param(12);
s1 = Param(13);
s2 = Param(14);
s3 = Param(15);
s4 = Param(16);

param = Param;


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