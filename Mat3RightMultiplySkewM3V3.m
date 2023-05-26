function result = Mat3RightMultiplySkewM3V3(left_mat3, right_skew)
result = zeros(3,3);
ofst = 1;

result(0 + ofst, 0 + ofst) = left_mat3(0 + ofst, 1 + ofst) * right_skew(2 + ofst) - left_mat3(0 + ofst, 2 + ofst) * right_skew(1 + ofst);
result(0 + ofst, 1 + ofst) = left_mat3(0 + ofst, 2 + ofst) * right_skew(0 + ofst) - left_mat3(0 + ofst, 0 + ofst) * right_skew(2 + ofst);
result(0 + ofst, 2 + ofst) = left_mat3(0 + ofst, 0 + ofst) * right_skew(1 + ofst) - left_mat3(0 + ofst, 1 + ofst) * right_skew(0 + ofst);
result(1 + ofst, 0 + ofst) = left_mat3(1 + ofst, 1 + ofst) * right_skew(2 + ofst) - left_mat3(1 + ofst, 2 + ofst) * right_skew(1 + ofst);
result(1 + ofst, 1 + ofst) = left_mat3(1 + ofst, 2 + ofst) * right_skew(0 + ofst) - left_mat3(1 + ofst, 0 + ofst) * right_skew(2 + ofst);
result(1 + ofst, 2 + ofst) = left_mat3(1 + ofst, 0 + ofst) * right_skew(1 + ofst) - left_mat3(1 + ofst, 1 + ofst) * right_skew(0 + ofst);
result(2 + ofst, 0 + ofst) = left_mat3(2 + ofst, 1 + ofst) * right_skew(2 + ofst) - left_mat3(2 + ofst, 2 + ofst) * right_skew(1 + ofst);
result(2 + ofst, 1 + ofst) = left_mat3(2 + ofst, 2 + ofst) * right_skew(0 + ofst) - left_mat3(2 + ofst, 0 + ofst) * right_skew(2 + ofst);
result(2 + ofst, 2 + ofst) = left_mat3(2 + ofst, 0 + ofst) * right_skew(1 + ofst) - left_mat3(2 + ofst, 1 + ofst) * right_skew(0 + ofst);


end