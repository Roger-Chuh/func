function testFisheyeRadTanThinPrismTilt()
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
          0.00001071922498768078,
          -0.0003963449586486949,
          -0.0002963449586486949,
          0.00901,
          -0.0009];
 
      
%       Param = [236.8820281-0;
% 236.7437341-0;
% 325.0065396;
% 241.9557923;
% 0.2124214875;
% -0.1719834148;
% 0.06211290885;
% -0.009899091649;
% 0;
% 0;
% 0;
% 0;
% 0;
% 0;
% 0;
% 0];
      
% Param(9:10) = 0;
% Param(11:end) = 0;
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
s50 = Param0(17);
s60 = Param0(18);
t10 = Param0(19);
t20 = Param0(20);


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
s5 = Param(17);
s6 = Param(18);
t1 = Param(19);
t2 = Param(20);

X = [-10.1234 -20.5678 15];
X_ = X;
X = [X(1)/X(3) X(2)/X(3) 1];
X = X./norm(X);

[pt2d, d_uv_d_pt3d, d_uv_d_param] = project(X_(1),X_(2), X_(3),fx0,fy0,cx0,cy0,k10,k20,k30,k40,k50,k60,p10,p20,s10,s20,s30,s40,s50,s60,t10,t20);

[pt3d, d_pt3d_d_uv, d_pt3d_d_param] = unproject(pt2d(1),pt2d(2), fx0,fy0,cx0,cy0,k10,k20,k30,k40,k50,k60,p10,p20,s10,s20,s30,s40,s50,s60,t10,t20);

 
err_ = X - pt3d

disturbance = [5.*(rand(4,1)-0.5);
    0.1.*(rand(6,1)-0.5);
    0.1.*(rand(10,1)-0.5)];
% disturbance(5:end) = 0;


disturbance_pt3d = 0.01.*(rand(3,1)-0.5);
disturbance_pt2d = 1000.*(rand(2,1)-0.5);
disturbance_pt2d = 10.*(rand(2,1)-0.5);

disturbance = [-2.48512711138509;-2.12291452309554;2.28265094210089;2.853457473397897;0.0166720512777738;-0.0339520395422974;0.0264749552794402;0.0158897114300060;0.0390732504151609;0.0396184685167673;-0.0192979630493664;0.0405478614969866;0.0146714295416348;0.0356074008212624;0.0196533331327743;-0.00643988179087626;0.0231564768498888;0.0176666421228288;0.0207302735927681;0.0338140440971284];
disturbance = disturbance./100;
disturbance_pt3d = 20.0.*[-0.00426187183200918;-0.000800120347325337;-0.00132981489423306];
disturbance_pt2d = [-3.77737164288323;3.75650983611374]./10;
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
s5 = Param(17);
s6 = Param(18);
t1 = Param(19);
t2 = Param(20);
[pt2d_err, d_uv_d_pt3d_err, d_uv_d_param_err] = project(X_(1),X_(2), X_(3),fx,fy,cx,cy,k1,k2,k3,k4,k5,k6,p1,p2,s1,s2,s3,s4,s5,s6,t1,t2);
[pt2d_err1, d_uv_d_pt3d_err1, d_uv_d_param_err1] = project(X_(1)+disturbance_pt3d(1),X_(2)+disturbance_pt3d(2), X_(3)+disturbance_pt3d(3), fx0,fy0,cx0,cy0,k10,k20,k30,k40,k50,k60,p10,p20,s10,s20,s30,s40,s50,s60,t10,t20);

err1 = pt2d_err' - pt2d';
err_comp1 = d_uv_d_param *disturbance;
check_Jac1 = [err1 - err_comp1 err1  err_comp1];

err11 = pt2d_err1' - pt2d';
err_comp11 = d_uv_d_pt3d *disturbance_pt3d;
check_Jac11 = [err11 - err_comp11 err11  err_comp11];


[pt3d_err, d_pt3d_d_uv_err, d_pt3d_d_param_err] = unproject(pt2d(1),pt2d(2), fx,fy,cx,cy,k1,k2,k3,k4,k5,k6,p1,p2,s1,s2,s3,s4,s5,s6,t1,t2);
[pt3d_err2, d_pt3d_d_uv_err2, d_pt3d_d_param_err2] = unproject(pt2d(1)+disturbance_pt2d(1),pt2d(2)+disturbance_pt2d(2), fx0,fy0,cx0,cy0,k10,k20,k30,k40,k50,k60,p10,p20,s10,s20,s30,s40,s50,s60,t10,t20);


err2 = pt3d_err' - pt3d';
err_comp2 = d_pt3d_d_param*disturbance;
check_Jac2 = [err2 - err_comp2 err2  err_comp2];


err22 = pt3d_err2' - pt3d';
err_comp22 = d_pt3d_d_uv*disturbance_pt2d;
check_Jac22 = [err22 - err_comp22 err22  err_comp22];

err_cat = [check_Jac1; check_Jac11];
% err_cat(3,:) = 0;
err_cat = [err_cat; check_Jac2; check_Jac22];


[pt2d_all_err1, d_uv_d_pt3d_all_err1, d_uv_d_param_all_err1] = project(X_(1)+disturbance_pt3d(1),X_(2)+disturbance_pt3d(2), X_(3)+disturbance_pt3d(3), fx,fy,cx,cy,k1,k2,k3,k4,k5,k6,p1,p2,s1,s2,s3,s4,s5,s6,t1,t2);
err3 = pt2d_all_err1' - pt2d';
err_comp3 = [d_uv_d_pt3d, d_uv_d_param] * [disturbance_pt3d; disturbance];
check_Jac3 = [err3 - err_comp3 err3  err_comp3];

[pt3d_all_err2, d_pt3d_d_uv_all_err2, d_pt3d_d_param_all_err2] = unproject(pt2d(1)+disturbance_pt2d(1),pt2d(2)+disturbance_pt2d(2), fx,fy,cx,cy,k1,k2,k3,k4,k5,k6,p1,p2,s1,s2,s3,s4,s5,s6,t1,t2);
err33 = pt3d_all_err2' - pt3d';
err_comp33 = [ d_pt3d_d_uv d_pt3d_d_param] * [disturbance_pt2d; disturbance];
check_Jac33 = [err33 - err_comp33 err33  err_comp33];


err_cat = [err_cat; check_Jac3; check_Jac33];






[xMat, yMat] = meshgrid(1:4: 640, 1:4: 480);
pix = [xMat(:) yMat(:)];
cov_2d = zeros(2,2,size(pix,1));
cov_3d = zeros(3,3,size(pix,1));
for i = 1 : size(pix,1)
    [pt3d_i, d_pt3d_d_uv_i, d_pt3d_d_param_i] = unproject(pix(i,1),pix(i,2), fx0,fy0,cx0,cy0,k10,k20,k30,k40,k50,k60,p10,p20,s10,s20,s30,s40,s50,s60,t10,t20);
    [pt2d, d_uv_d_pt3d, d_uv_d_param] = project(pt3d_i(1),pt3d_i(2), pt3d_i(3),fx0,fy0,cx0,cy0,k10,k20,k30,k40,k50,k60,p10,p20,s10,s20,s30,s40,s50,s60,t10,t20);
    cov_2d(:,:,i) = d_pt3d_d_uv_i'*d_pt3d_d_uv_i;
    cov_3d(:,:,i) = d_uv_d_pt3d'*d_uv_d_pt3d;
end

cov_2d_11 = reshape(cov_2d(1,1,:), size(xMat));
cov_2d_12 = reshape(cov_2d(1,2,:), size(xMat));
cov_2d_21 = reshape(cov_2d(2,1,:), size(xMat));
cov_2d_22 = reshape(cov_2d(2,2,:), size(xMat));

asg = 1;

figure,surface(abs(cov_2d_11));

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

function [pt2d, d_uv_d_pt3d, d_uv_d_param] = project(X,Y,Z, fx,fy,cx,cy,k1,k2,k3,k4,k5,k6,p1,p2,s1,s2,s3,s4,s5,s6,t1,t2)
global Param
d_uv_d_pt3d = []; d_uv_d_param = [];

% X = pt3d(:,1);
% Y = pt3d(:,2);
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
function [pt3d, d_pt3d_d_uv, d_pt3d_d_param] = unproject(u,v, fx,fy,cx,cy,k1,k2,k3,k4,k5,k6,p1,p2,s1,s2,s3,s4,s5,s6,t1,t2)
global Param data

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
t1 = param(17);
t2 = param(18);

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