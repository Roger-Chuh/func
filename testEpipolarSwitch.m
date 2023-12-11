function testEpipolarSwitch()

num = 5;


K = [400 0 320; 0 400 240; 0 0 1];

T = [rodrigues([0.1 -0.2 0.3]) [100 200 -300]'; 0 0 0 1];

xyz = [200 300 1200];
xyz = xyz./xyz(3);
depth = 1200;


host_pix = pflat(K*[xyz]');
host_pix_norm = inv(K)*host_pix;
host_pix = host_pix(1:2)';
[target_pix, tgtPt3d] = TransformAndProject(depth.*xyz, K, T(1:3, 1:3), T(1:3, 4));
[target_pix_err, tgtPt3d] = TransformAndProject((depth+100).*xyz, K, T(1:3, 1:3), T(1:3, 4));
target_pix = target_pix_err;
if 0
    target_pix = target_pix + [1 1];
end
target_pix_norm = inv(K)*[target_pix 1]';


F = inv(K')*SkewSymMat(T(1:3,4)) * T(1:3,1:3)*inv(K);
E = SkewSymMat(T(1:3,4)) * T(1:3,1:3);
epiline_host_in_target = F*[host_pix 1]'; %这是target帧上极线的方向，但是极限的方向是用host的pixel坐标计算得到
epiline_host_in_target = epiline_host_in_target./norm(epiline_host_in_target(1:2));
target_pix_err = dot([target_pix 1], epiline_host_in_target');

epiline_host_in_target_norm = E*host_pix_norm;
epiline_host_in_target_norm = epiline_host_in_target_norm./norm(epiline_host_in_target_norm(1:2));
target_pix_err_norm =  dot(target_pix_norm, epiline_host_in_target_norm);




epiline_target_in_host = F'*[target_pix 1]';
epiline_target_in_host = epiline_target_in_host./norm(epiline_target_in_host(1:2));
host_pix_err = dot([host_pix 1], epiline_target_in_host');

epiline_target_in_host_norm = E'*target_pix_norm;
epiline_target_in_host_norm = epiline_target_in_host_norm./norm(epiline_target_in_host_norm(1:2));
host_pix_err_norm =  dot(host_pix_norm, epiline_target_in_host_norm);



for i = 1 : num
Norm(:,i) = rand(3,1);
Norm(:,i) = Norm(:,i)./norm(Norm(:,i));
end



rotMat = rodrigues(rand(3,1));
rotAxis = [0.1 0.2 0.3]';
rotAxis = rotAxis./norm(rotAxis);

for i = 1:num
Norm(:,i) = rodrigues(rodrigues(deg2rad(80*i)*rotAxis)*rotMat);
Norm(:,i) = Norm(:,i)./norm(Norm(:,i));
end
NN = Norm';
NN(11,:) = [0 0 0];
N = Norm';
N(:,4) = 1;
[A1, B1, C1] = svd(N,0);
[A, B, CC] = svd(Norm',0);
[x,y,z] = svd(Norm*Norm');
[xx,yy,zz] = svd(NN);
NN(:,4) = 1;
[xxx,yyy,zzz] = svd(NN);
figure,plotQuiver(Norm')
rank(Norm*Norm');
y;
[y(2,2) sum(sum(y)) y(2,2)/sum(sum(y))]
end