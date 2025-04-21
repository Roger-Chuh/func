function ret = IntegrateJlTau(w, t)
omega = norm(w);
eps = 1e-8;

if (omega < eps)
    ret = 0.5 * t * t * eye(3);
    return ;
end

omega_sq = omega * omega;
omega_t = omega * t;
sin_omega_t = sin(omega_t);
cos_omega_t = cos(omega_t);

skew_w = SkewSymMat(w);
skew_w_sq = skew_w * skew_w;

inv_omega_sq = 1.0 / omega_sq;
inv_omega_quad = inv_omega_sq * inv_omega_sq;
t_sq_half = 0.5 * t * t;

term3_part1 = t_sq_half * inv_omega_sq;
term3_part2 = (cos_omega_t - 1) * inv_omega_quad;
cos_term = term3_part1 + term3_part2;

sin_term = (t - sin_omega_t / omega) * inv_omega_sq;

ret = t_sq_half * eye(3) + sin_term * skew_w + cos_term * skew_w_sq;
end