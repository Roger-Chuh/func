function testTrifocalTensor()
intrMat = [500 0 320; 0 500 240;0 0 1];

xyz = [10 20 1000; 10 22 1100; 15 20 1020; -10 30 900; 17 17 800; 20 -12 950];

Tc1w = [rodrigues([0.01 0.02 -0.03]) [10 20 -30]';0 0 0 1];
Tc2w = [rodrigues([0.01 -0.02 0.03]) [-10 20 -30]';0 0 0 1];
Tc3w = [rodrigues([-0.01 0.02 -0.03]) [10 -20 30]';0 0 0 1];


pix1 = TransformAndProject(xyz, intrMat,Tc1w(1:3,1:3), Tc1w(1:3,4));
pix2 = TransformAndProject(xyz, intrMat,Tc2w(1:3,1:3), Tc2w(1:3,4));
pix3 = TransformAndProject(xyz, intrMat,Tc3w(1:3,1:3), Tc3w(1:3,4));


trace = {};
for i = 1 : 6
    trace{i,1}.a1 = inv(intrMat)*[pix1(i,:) 1]';
    trace{i,1}.a2 = inv(intrMat)*[pix2(i,:) 1]';
    trace{i,1}.a3 = inv(intrMat)*[pix3(i,:) 1]';
end

Tri = getTrifocal(trace{1}, trace{2}, trace{3}, trace{4}, trace{5}, trace{6});
% Tri = Tri./norm(Tri(:));

p3 = pointTransfer(Tri,trace{2,1}.a1, trace{2,1}.a2);
p3 = (p3/p3(3));

T_2_1 = inv(Tc2w * inv(Tc1w));
% T_2_1(1:3,4) = T_2_1(1:3,4)./norm(T_2_1(1:3,4));
T_3_1 = inv(Tc3w * inv(Tc1w));
% T_3_1(1:3,4) = T_3_1(1:3,4)./norm(T_3_1(1:3,4));
T=TFT_from_P(intrMat * eye(3,4),intrMat * T_2_1(1:3,:),intrMat * T_3_1(1:3,:));


% 'add' calibration from the tensor
T_pix= transform_TFT(T,inv(intrMat),inv(intrMat),inv(intrMat), 0);


[ point3 ] = pointTransfer( T, inv(intrMat)*[pix1 1]', inv(intrMat)*[pix2 1]' );
point3 = point3./point3(3);
pflat(intrMat * point3)
end
function T=TFT_from_P(P1,P2,P3)
%TFT_FROM_P Trifocal tensor from the projection matrices.
%
%  General formula to compute the trifocal tensor from any three projection
%  matrices.
%

T=zeros(3,3,3);
for i=1:3
    for j=1:3
        for k=1:3
            T(j,k,i)=(-1)^(i+1)*det([P1([1:(i-1) (i+1):3],:);P2(j,:);P3(k,:)]);
        end
    end
end
T=T/norm(T(:));
end
function T_new = transform_TFT(T_old,M1,M2,M3,inverse)
% Tranformed TFT
%
% short function to transform the TFT when an algebraic transformation
% has been aplied to the image points.
%
% if inverse==0 :
%   from a TFT T_old assossiated to P1_old, P2_old, P3_old, find the new TFT
%   T_new associated to P1_new=M1*P1_old, P2_new=M2*P2_old, P3_new=M3*P3_old.
% if inverse==1 :
%   from a TFT T_old assossiated to P1_old, P2_old, P3 _old, find the new TFT
%   T_new associated to M1*P1_new=P1_old, M2*P2_new=P2_old, M3*P3_new=P3_old.
%

% Copyright (c) 2017 Laura F. Julia <laura.fernandez-julia@enpc.fr>
% All rights reserved.
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


if nargin<5
    inverse=0;
end

if inverse==0
    M1i=inv(M1); T_new=zeros(3,3,3);
    T_new(:,:,1)=M2*(M1i(1,1)*T_old(:,:,1) + M1i(2,1)*T_old(:,:,2) + M1i(3,1)*T_old(:,:,3) )*M3.';
    T_new(:,:,2)=M2*(M1i(1,2)*T_old(:,:,1) + M1i(2,2)*T_old(:,:,2) + M1i(3,2)*T_old(:,:,3) )*M3.';
    T_new(:,:,3)=M2*(M1i(1,3)*T_old(:,:,1) + M1i(2,3)*T_old(:,:,2) + M1i(3,3)*T_old(:,:,3) )*M3.';

elseif inverse==1
    M2i=inv(M2); M3i=inv(M3); T_new=zeros(3,3,3);
    T_new(:,:,1)=M2i*(M1(1,1)*T_old(:,:,1) + M1(2,1)*T_old(:,:,2) + M1(3,1)*T_old(:,:,3) )*M3i.';
    T_new(:,:,2)=M2i*(M1(1,2)*T_old(:,:,1) + M1(2,2)*T_old(:,:,2) + M1(3,2)*T_old(:,:,3) )*M3i.';
    T_new(:,:,3)=M2i*(M1(1,3)*T_old(:,:,1) + M1(2,3)*T_old(:,:,2) + M1(3,3)*T_old(:,:,3) )*M3i.';
end

T_new=T_new/norm(T_new(:));

end
function [ point3 ] = pointTransfer( Tri, point1, point2 )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    point3 = zeros(3,1);
    i = 1;
    j = 2;
    
    for l=1:3,
        v1 = 0;
        v2 = 0;
        for k=1:3,
            v1 = v1+point1(k)*Tri(k,j,l);
            v2 = v2+point1(k)*Tri(k,i,l);
        end;
        
        point3(l) = point2(i)*v1 - point2(j)*v2;
    end;
end
function [Tri ] = getTrifocal( m1, m2, m3, m4, m5, m6 )
%m1, m2, m3, m4, m5, m6 are the 6 matches(structs) input
%   mi.a1 mi.a2 mi.a3 are the 3 image points
    
    lambda1 = inv([m1.a1, m2.a1, m3.a1])*m4.a1;
    B1 = inv([lambda1(1)*m1.a1, lambda1(2)*m2.a1, lambda1(3)*m3.a1]);
    %{
    fprintf('\n***** 1 ****\n\n');
    lambda1
    m1.a1
    m2.a1
    m3.a1
    m4.a1
    B1
    %}
    lambda2 = inv([m1.a2, m2.a2, m3.a2])*m4.a2;
    B2 = inv([lambda2(1)*m1.a2, lambda2(2)*m2.a2, lambda2(3)*m3.a2]);
    %{
    fprintf('\n***** 2 ****\n\n');
    lambda2
    m1.a2
    m2.a2
    m3.a2
    m4.a2
    B2
    %}
    lambda3 = inv([m1.a3, m2.a3, m3.a3])*m4.a3;
    B3 = inv([lambda3(1)*m1.a3, lambda3(2)*m2.a3, lambda3(3)*m3.a3]);
    %{
    fprintf('\n***** 3 ****\n\n');
    lambda3
    m1.a3
    m2.a3
    m3.a3
    m4.a2
    B3
    %}
    X5 = struct;
    X5.a1 = B1*m5.a1;
    X5.a2 = B2*m5.a2;
    X5.a3 = B3*m5.a3;
    
    X6 = struct;
    X6.a1 = B1*m6.a1;
    X6.a2 = B2*m6.a2;
    X6.a3 = B3*m6.a3;
    
    Matrix = zeros(3, 5); 
    x5 = X5.a1(1);
    y5 = X5.a1(2);
    w5 = X5.a1(3);
    
    x6 = X6.a1(1);
    y6 = X6.a1(2);
    w6 = X6.a1(3);
    
    Matrix(1,:) = [-x5*y6 + x5*w6, x6*y5 - y5*w6, -x6*w5 + y6*w5, ...
                    -x5*w6 + y5*w6, x5*y6 - y6*w5 ];
                
    x5 = X5.a2(1);
    y5 = X5.a2(2);
    w5 = X5.a2(3);
    
    x6 = X6.a2(1);
    y6 = X6.a2(2);
    w6 = X6.a2(3);
    
    Matrix(2,:) = [-x5*y6 + x5*w6, x6*y5 - y5*w6, -x6*w5 + y6*w5, ...
                    -x5*w6 + y5*w6, x5*y6 - y6*w5 ];
    x5 = X5.a3(1);
    y5 = X5.a3(2);
    w5 = X5.a3(3);
    
    x6 = X6.a3(1);
    y6 = X6.a3(2);
    w6 = X6.a3(3);
    
    Matrix(3,:) = [-x5*y6 + x5*w6, x6*y5 - y5*w6, -x6*w5 + y6*w5, ...
                    -x5*w6 + y5*w6, x5*y6 - y6*w5 ];
                    
    [~, ~, V] = svd(Matrix);
    
    T1 = V(:, end-1);
    T2 = V(:, end);
    
    a1 = T1(1);
    a2 = T1(2);
    a3 = T1(3);
    a4 = T1(4);
    a5 = T1(5);
    
    b1 = T2(1);
    b2 = T2(2);
    b3 = T2(3);
    b4 = T2(4);
    b5 = T2(5);
    
    polynomial = Fun(a1,b1,a2,b2,a5,b5)-Fun(a2,b2,a3,b3,a5,b5) ...
                -Fun(a2,b2,a4,b4,a5,b5)-Fun(a1,b1,a3,b3,a4,b4) ...
                +Fun(a2,b2,a3,b3,a4,b4)+Fun(a3,b3,a4,b4,a5,b5);
    
    alphas = roots(polynomial);
    alphas = alphas(imag(alphas) == 0);
    
    T = T1 + alphas(1)*T2;
   
    XbyW = (T(4)-T(5))/(T(2)-T(3));
    YbyW = T(4)/(T(1)-T(3));
    ZbyW = T(5)/(T(1)-T(2));
    
    Matrix1 = zeros(4,4);
    
    x5 = X5.a1(1);
    y5 = X5.a1(2);
    w5 = X5.a1(3);
    
    x6 = X6.a1(1);
    y6 = X6.a1(2);
    w6 = X6.a1(3);
    
    Matrix1(1,:) = [w5, 0, -x5, w5-x5];
    Matrix1(2,:) = [0, w5, -y5, w5-y5];
    Matrix1(3,:) = [w6*XbyW, 0, -x6*ZbyW, w6-x6];
    Matrix1(4,:) = [0, w6*YbyW, -y6*ZbyW, w6-y6];
    
    %{
    fprintf('\n***** Matrix1 *****\n\n');
    Matrix1
    %}
    Matrix2 = zeros(4,4);
    
    x5 = X5.a2(1);
    y5 = X5.a2(2);
    w5 = X5.a2(3);
   
    x6 = X6.a2(1);
    y6 = X6.a2(2);
    w6 = X6.a2(3);
    
    Matrix2(1,:) = [w5, 0, -x5, w5-x5];
    Matrix2(2,:) = [0, w5, -y5, w5-y5];
    Matrix2(3,:) = [w6*XbyW, 0, -x6*ZbyW, w6-x6];
    Matrix2(4,:) = [0, w6*YbyW, -y6*ZbyW, w6-y6];
    
    %{
    fprintf('\n***** Matrix2 *****\n\n');
    Matrix2
    %}
    Matrix3 = zeros(4,4);
    
    x5 = X5.a3(1);
    y5 = X5.a3(2);
    w5 = X5.a3(3);
    
    x6 = X6.a3(1);
    y6 = X6.a3(2);
    w6 = X6.a3(3);
    
    Matrix3(1,:) = [w5, 0, -x5, w5-x5];
    Matrix3(2,:) = [0, w5, -y5, w5-y5];
    Matrix3(3,:) = [w6*XbyW, 0, -x6*ZbyW, w6-x6];
    Matrix3(4,:) = [0, w6*YbyW, -y6*ZbyW, w6-y6];

    %{
    fprintf('\n***** Matrix3 *****\n\n');
    Matrix3
    %}
    
    [~, ~, V] = svd(Matrix1);
    abgd1 = V(:,end);
    
    [~, ~, V] = svd(Matrix2);
    abgd2 = V(:,end);
    
    [~, ~, V] = svd(Matrix3);
    abgd3 = V(:,end);
    
    P1 = zeros(3, 4);
    
    P1(1,:) = [abgd1(1), 0, 0, abgd1(4)];
    P1(2,:) = [0, abgd1(2), 0, abgd1(4)];
    P1(3,:) = [0, 0, abgd1(3), abgd1(4)];
    
    P2 = zeros(3, 4);
    
    P2(1,:) = [abgd2(1), 0, 0, abgd2(4)];
    P2(2,:) = [0, abgd2(2), 0, abgd2(4)];
    P2(3,:) = [0, 0, abgd2(3), abgd2(4)];
    
    P3 = zeros(3, 4);
    
    P3(1,:) = [abgd3(1), 0, 0, abgd3(4)];
    P3(2,:) = [0, abgd3(2), 0, abgd3(4)];
    P3(3,:) = [0, 0, abgd3(3), abgd3(4)];
    
    P1 = B1\P1;
    P2 = B2\P2;
    P3 = B3\P3;
    
    H = zeros(4, 4);
    H(1:3,:) = P1;
    H(4,:) = cross4(H(1,:)', H(2,:)', H(3,:)');
    
    H = inv(H);
    
    P2 = P2*H;
    P3 = P3*H;
    P1 = P1*H;
    
    Tri = zeros(3, 3, 3);
    for i=1:3,
        for j=1:3,
            for k=1:3,
                Tri(i,j,k) = P2(j,i)*P3(k,4) - P2(j,4)*P3(k,i);
            end;
        end;
    end;
    
end