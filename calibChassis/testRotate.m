function testRotate()

w = 640;
h = 480;
r = rodrigues([0.01; 0.002; -0.003]);
t = [-300; 20; -3]; [0;0;0];

K = [500 0 319; 0 500 251; 0 0 1];
K_ = [K(1,1) 0 1+h-K(2,3); 0 K(1,1) K(1,3);0 0 1];

KK = [520 0 327; 0 520 259; 0 0 1];
KK_ = [KK(1,1) 0 1+h-KK(2,3); 0 KK(1,1) KK(1,3);0 0 1];

z1 = 1010;
pt1 = [300 200];

pt1_ = [1+h-pt1(2) pt1(1)];

metric1 = inv(K)*[pt1 1]';
xyz1 = z1.*metric1;

xyz2 = r*xyz1 + t;
pt2 = pflat(KK*xyz2);

xyz1_0 = rotz(90)*xyz1;

t_ = rotz(90)*t;
% t_ = [-r(:,2) r(:,1) r(:,3)]*t;
r_ = rotz(90)*r;

r__ = [-r_(:,2) r_(:,1) r_(:,3)];

metric1_ = inv(K_)*[pt1_ 1]';
xyz1_ = z1.*metric1_;

xyz2__ = r_*xyz1 + t_;
xyz2_ = r__*xyz1_ + t_;
pt2_ = pflat(KK_*xyz2_);

end