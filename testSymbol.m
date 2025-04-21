function testSymbol()
% syms x y z a1 a2 a3 p1 p2 p3
% expression = inv((SkewSymMat([x y z]) * SkewSymMat([x y z]) + SkewSymMat([a1 a2 a3]) * SkewSymMat([a1 a2 a3]))) * (SkewSymMat([x y z])* SkewSymMat([x y z])) * [p1;p2;p3]
% simplifiedExpr = simplify(expression)




syms Rbc1_00 Rbc1_01 Rbc1_02 Rbc1_10 Rbc1_11 Rbc1_12 Rbc1_20 Rbc1_21 Rbc1_22 ...
     Rbc2_00 Rbc2_01 Rbc2_02 Rbc2_10 Rbc2_11 Rbc2_12 Rbc2_20 Rbc2_21 Rbc2_22 ...
     Rbc3_00 Rbc3_01 Rbc3_02 Rbc3_10 Rbc3_11 Rbc3_12 Rbc3_20 Rbc3_21 Rbc3_22 ...
     tbc1_0 tbc1_1 tbc1_2 ...
     tbc2_0 tbc2_1 tbc2_2 ...
     tbc3_0 tbc3_1 tbc3_2 ...
     x1 y1 z1 ...
     x2 y2 z2 ...
     x3 y3 z3

Rbc1 = [Rbc1_00 Rbc1_01 Rbc1_02; Rbc1_10 Rbc1_11 Rbc1_12; Rbc1_20 Rbc1_21 Rbc1_22];
Rbc2 = [Rbc2_00 Rbc2_01 Rbc2_02; Rbc2_10 Rbc2_11 Rbc2_12; Rbc2_20 Rbc2_21 Rbc2_22];
Rbc3 = [Rbc3_00 Rbc3_01 Rbc3_02; Rbc3_10 Rbc3_11 Rbc3_12; Rbc3_20 Rbc3_21 Rbc3_22];

Rbc1_transpose = [Rbc1_00 Rbc1_10 Rbc1_20; Rbc1_01 Rbc1_11 Rbc1_21; Rbc1_02 Rbc1_12 Rbc1_22];
Rbc2_transpose = [Rbc2_00 Rbc2_10 Rbc2_20; Rbc2_01 Rbc2_11 Rbc2_21; Rbc2_02 Rbc2_12 Rbc2_22];
Rbc3_transpose = [Rbc3_00 Rbc3_10 Rbc3_20; Rbc3_01 Rbc3_11 Rbc3_21; Rbc3_02 Rbc3_12 Rbc3_22];

tbc1 = [tbc1_0; tbc1_1; tbc1_2];
tbc2 = [tbc2_0; tbc2_1; tbc2_2];
tbc3 = [tbc3_0; tbc3_1; tbc3_2];

X1 = [x1; y1; z1];
X1_transpose = [x1 y1 z1];
X2 = [x2; y2; z2];
X2_transpose = [x2 y2 z2];
X3 = [x3; y3; z3];

R21 = Rbc2_transpose * Rbc1;
R21_transpose = Rbc1_transpose * Rbc2;
t21 = Rbc2_transpose * (tbc1 - tbc2);
t21_transpose = [t21(1) t21(2) t21(3)];
R31 = Rbc3_transpose * Rbc1;
t31 = Rbc3_transpose * (tbc1 - tbc3);
t31_transpose = [t31(1) t31(2) t31(3)];

predict = (t31 * X1_transpose * R21_transpose - R31 * X1 * t21_transpose) * SkewSymMat(X2) * SkewSymMat(t21) * R21 * X1;
simplifiedExpr = simplify(predict);
end