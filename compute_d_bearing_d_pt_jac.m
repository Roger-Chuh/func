function d_bearing_d_pt = compute_d_bearing_d_pt_jac(pt)
d_bearing_d_pt = (norm(pt).*eye(3) - pt * pt'./norm(pt))./(norm(pt)^2);
end