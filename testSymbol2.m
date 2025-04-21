function testSymbol2()

% syms wx wy wz awx awy awz
% 
% w = [wx; wy; wz];
% aw = [awx; awy; awz];
% 
% 
% y = 2 * SkewSymMat(w) * SkewSymMat(aw) +  SkewSymMat(aw) * SkewSymMat(w);
% 
% 
% w1 = rand(3,1);
% aw1 = rand(3,1);
% y1 = 2 * SkewSymMat(w1) * SkewSymMat(aw1) +  SkewSymMat(aw1) * SkewSymMat(w1);





%% 20250417/1514
% syms w1 w2 w3 a1 a2 a3 R11 R12 R13 R21 R22 R23 R31 R32 R33 dt
% w = [w1; w2; w3];
% a = [a1; a2; a3];
% R_T = [R11 R21 R31; R12 R22 R32; R13 R23 R33];
% % R_T = R';
% % A = zeros(18, 18);
% % R = rodrigues(rand(3, 1));
% % w = rand(3,1);
% % dt = 0.01;
% % a = rand(3,1);
% 
% A(1:3,4:6) = R_T;
% A(4:6,7:9) = eye(3);
% A(10:12, 10:12) = SkewSymMat(-R_T * w);
% A(10:12, 13:15) = eye(3);
% A(13:15, 13:15) = SkewSymMat(-R_T * w);
% A(13:15, 16:18) = eye(3);
% A(16:18, 4:6) = -R_T * SkewSymMat(a);
% 
% % F = expm(A * dt);
% 
% 
% AA = zeros(6, 6);
% A(1,2) = a
% 
% 
% F = eye(18) + A * dt + 0.5 * A * A * dt * dt;



syms AA BB CC DD dt


AAA(1,2) = AA;
AAA(2,3) = BB;
AAA(4,4) = CC;
AAA(4, 5) = BB;
AAA(5,5) = CC;
AAA(5,6) = BB;
AAA(6,2) = DD;

 F = eye(6) + AAA * dt + 0.5 * AAA * AAA * dt * dt;

F_simple = simplify(F);

end