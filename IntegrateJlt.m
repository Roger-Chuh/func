function intg_Jl_t = IntegrateJlt(w, t)
w_norm = norm(w);

intg_Jl_t2 = 0.5 * t^2 * eye(3) + SkewSymMat(w)/w_norm^2 * (t - sin(w_norm * t)/w_norm) + (t^2/(2 * w_norm^2) + (cos(w_norm * t)-1)/w_norm^4) * SkewSymMat(w) * SkewSymMat(w);

intg_Jl_t = w * w' * t^2/(2 * w_norm * w_norm) + SkewSymMat(w) * t / (w_norm * w_norm) + (-cos(w_norm * t)/(w_norm * w_norm) + 1/(w_norm * w_norm)) * eye(3) - (-w * w' * cos(w_norm * t) / (w_norm^4) + w * w' / (w_norm^4))...
    - (SkewSymMat(w) * sin(w_norm * t) / (w_norm^3));

intg_Jl_t - intg_Jl_t2
end