function CalibB2C()

close all;


R1 = rodrigues([0.1 0.2 0.3]);
R2 = rodrigues([0.1 0.2 -0.2]);
R3= R2*R1;
axis3 = rodrigues(R3)./norm(rodrigues(R3));
axis1 = rodrigues(R1)./norm(rodrigues(R1));

% axis2 = rodrigues(R2)./norm(rodrigues(R2));


b2c = rodrigues([-0.1 0.2 0.3]);

R3 = b2c*R1/b2c;

axis3 = rodrigues(R3)./norm(rodrigues(R3));
axis1 = rodrigues(R1)./norm(rodrigues(R1));

axis3 - b2c*axis1;



% R1_list = {rodrigues([0.1 0.2 0.3]); rodrigues([-0.3 0.2 0.3]); rodrigues([0.15 -0.2 0.3]); rodrigues([0.1 0.27 -0.33])};

R1_list = {rodrigues([0.1 0.2 0.3]);
    rodrigues([0.11 0.22 0.33]);
    rodrigues([0.15 0.19 0.25]);
    rodrigues([0.13 0.17 0.29]);
    rodrigues([0.13 0.21 0.28])};

b2c_list = {rodrigues([-0.3 0.22 -0.13]);
    rodrigues([0.11 -0.202 0.33]);
    rodrigues([-0.15 0.19 0.205]);
    rodrigues([0.13 -0.17 0.29]);
    rodrigues([0.103 0.21 -0.28])};
%     rodrigues([0.123 0.21 0.28]);
%     rodrigues([0.113 0.21 0.28]);
%     rodrigues([0.123 0.21 0.28])};


rot_axis = [-0.55; 0.6; -0.2];
rot_axis = rot_axis./norm(rot_axis);
axis_cross_mat2 =[];
for i = 1 : length(b2c_list)
    rotAxis(:,i) =  b2c_list{i}*rot_axis;
    cros_m = cross(rotAxis(:,i), rot_axis);
    cros_M(:,i) = cros_m./norm(cros_m);
    [axis_cross_mat2(:,i), rotAxis2(:,i)] = genData2(rotAxis(:,i), b2c);
end

[~,~,V3] = svd(pextend(cros_M)');
plane3 = [V3(:,4)];
plane3 = plane3./norm(plane3(1:3));

plane3 = sign(rot_axis(1))*sign(plane3(1)).*plane3;
error = plane3(1:3) - rot_axis
cros_M2 = cros_M;
cros_M2(4,:) = 0;
vert_van_dir_err3 = ((dot(repmat(plane3,1,size(cros_M2,2)), cros_M2)));

figure,plotQuiver(cros_M','r')







saddkj = 1;



axis_cross_mat22 = [axis_cross_mat2];  [[0;0;0]];
[~,~,V] = svd(pextend(axis_cross_mat22)');
plane = [V(:,4)];
plane = plane./norm(plane(1:3));

axis_cross_mat222 = axis_cross_mat22;
axis_cross_mat222(4,:) = 0;

vert_van_dir_err = ((dot(repmat(plane,1,size(axis_cross_mat222,2)), axis_cross_mat222)));

% figure,plot3(axis_cross_mat2(1,:), axis_cross_mat2(2,:), axis_cross_mat2(3,:),'or');axis equal;

figure,plotQuiver(axis_cross_mat2','r');
figure,plotQuiver(rotAxis','r');

figure,plot3(rotAxis(1,:), rotAxis(2,:), rotAxis(3,:),'or');hold on;
plot3(rotAxis2(1,:), rotAxis2(2,:), rotAxis2(3,:),'ob');
axis equal;





baseR_list = {};

% baseR_list{1,1} = eye(3); R1_list{1,1};
baseR_list{1,1} = R1_list{1}';
baseR_list{2,1} = (R1_list{2}*R1_list{1})';
baseR_list{3,1} = (R1_list{3}*R1_list{2}*R1_list{1})';
baseR_list{4,1} = (R1_list{4}*R1_list{3}*R1_list{2}*R1_list{1})';
baseR_list{5,1} = (R1_list{5}*R1_list{4}*R1_list{3}*R1_list{2}*R1_list{1})';
% for i = 2 : length(R1_list)
%     baseR_list{i,1} = baseR_list{i,1}*R1_list{i-1}*R1_list{i};   % (R1_list{i,1})*R1_list{1,1};
% end

R1_list_check = {};
R1_list_check{1,1} = baseR_list{1,1}' * eye(3);
for i = 2 : length(baseR_list)
    R1_list_check{i,1} = baseR_list{i,1}'*baseR_list{i-1,1};
    
end


axis1_mat = [];
R3_list = {};
axis3_mat = [];
axis_cross_mat = [];
err = [];

for i = 1 : length(R1_list)
    if i == 1
        [R3_list{i,1}, axis3_mat(:,i), axis1_mat(:,i),axis_cross_mat(:,i), err(:,i)] = genData(R1_list{i}, b2c, eye(3));
    else
        [R3_list{i,1}, axis3_mat(:,i), axis1_mat(:,i),axis_cross_mat(:,i), err(:,i)] = genData(R1_list{i}, b2c, baseR_list{i-1});
    end
    %         [R3_list{i,1}, axis3_mat(:,i), axis1_mat(:,i),axis_cross_mat(:,i), err(:,i)] = genData(baseR_list{i}, b2c, baseR_list{i});
end

axis_cross_mat = [axis_cross_mat];


[~,~,V2] = svd(pextend(axis_cross_mat)');
plane2 = [V2(:,4)];
plane2 = plane2./norm(plane2(1:3));


vert_van_dir_err = ((dot(repmat(plane2,1,size(axis_cross_mat,2)), pextend(axis_cross_mat))));


figure,plotQuiver(axis_cross_mat','r');
figure,plot3(axis_cross_mat(1,:), axis_cross_mat(2,:), axis_cross_mat(3,:),'or');axis equal;
figure,plot3(axis1_mat(1,:), axis1_mat(2,:), axis1_mat(3,:),'or');hold on;
plot3(axis3_mat(1,:), axis3_mat(2,:), axis3_mat(3,:),'ob');
axis equal;


end
function [R3, axis3, axis1,axis_cross, err] = genData(R1, b2c,baseR)

cur2base = baseR*R1;
R3 = b2c*R1/b2c;
if 0
    axis3 = cur2base'*rodrigues(R3)./norm(rodrigues(R3));
    axis1 = cur2base'*rodrigues(R1)./norm(rodrigues(R1));
else
    axis3 = rodrigues(R3)./norm(rodrigues(R3));
    axis1 = rodrigues(R1)./norm(rodrigues(R1));
end
err = axis3 - b2c*axis1;
axis_cross = cross(axis3,axis1);
axis_cross = axis_cross./norm(axis_cross);
end
function [axis_cross, rotAxis2] = genData2(rotAxis, b2c)

rotAxis2 = b2c*rotAxis;
axis_cross = cross(rotAxis,rotAxis2);
axis_cross = axis_cross./norm(axis_cross);
end
function plotQuiver(vec,color)

xyz = zeros(size(vec));

quiver3(xyz(:,1),xyz(:,2),xyz(:,3),vec(:,1),vec(:,2),vec(:,3),color);
end