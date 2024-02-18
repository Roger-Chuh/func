function testWeightCombination()
gradWeight = rand(10,1);
r_vec = rand(10,1);
r2_vec_ = r_vec.^2;


weight_vec_ = gradWeight .* (5 ./ (4 + r2_vec_ .* gradWeight));

weight_vec_./gradWeight


r_vec' * diag(ones(10,1)) * r_vec
r_vec' * r_vec
end