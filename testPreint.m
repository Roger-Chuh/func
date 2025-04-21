function testPreint()

if 0
    rot_vec = rand(3,1);
    rodrigues(rot_vec) - expm(SkewSymMat(rot_vec))
    
    
    
    w = [1 2 3]';
    a = [7 4 9]';
    
    
    dt = 0.01;
    
    
    num = 15 + 3 + 6 + 6;
    
    A = zeros( num, num );
    
    A(1:3,1:3) = - skew(w);
    A(1:3,10:12) = -eye(3);
    A(1:3,16:18) = - skew(w);
    A(1:3,19:24) = func_K(w);
    
    A(4:6, 1:3) = - skew(a);
    A(4:6, 4:6) = - skew(w);
    A(4:6, 13:15) = - eye(3);
    A(4:6, 25:30) = - func_K(a);
    
    A(7:9, 4:6) = eye(3);
    A(7:9,7:9) = - skew(w);
    
    
    expm(A*dt)
end


if 0
    w = rand(3, 1);
    a = rand(3, 1);
    num = 15;
    A = zeros( num, num );
    A(1:3,1:3) = - skew(w);
    A(1:3,10:12) = -eye(3);
    % A(1:3,16:18) = - skew(w);
    
    A(4:6, 1:3) = - skew(a);
    A(4:6, 4:6) = - skew(w);
    A(4:6, 13:15) = - eye(3);
    
    A(7:9, 4:6) = eye(3);
    A(7:9,7:9) = - skew(w);
    dt = 0.01;
    AA = expm(A*dt);
    
    
    rodrigues(w * dt)' - AA(1:3, 1:3)
    
    
    IntegrateJlt(w, dt);
end


% syms wx wy wz ax ay az w a dt;
% w = [wx; wy ;wz];
% a = [ax; ay ;az];

w = rand(3,1);
a = rand(3,1);
num = 15;
A = zeros( num, num );
A(1:3,1:3) = - skew(w);
A(1:3,10:12) = -eye(3);
% A(1:3,16:18) = - skew(w);

A(4:6, 1:3) = - skew(a);
A(4:6, 4:6) = - skew(w);
A(4:6, 13:15) = - eye(3);

A(7:9, 4:6) = eye(3);
A(7:9,7:9) = - skew(w);


dt = 0.01
AA = expm(A * dt);
B = zeros(15, 12);
B(1:3,1:3) = -eye(3);
B(4:6,4:6) = -eye(3);
B(10:12,7:9) = eye(3);
B(13:15,10:12) = eye(3);

BB = AA * B * dt;


end




function intg_Jl_t = IntegrateJlt(w, t)
w_norm = norm(w);

intg_Jl_t = 0.5 * t^2 * eye(3) + SkewSymMat(w)/w_norm^2 * (t - sin(w_norm * t)/w_norm) + (t^2/(2 * w_norm^2) + (cos(w_norm * t)-1)/w_norm^4) * SkewSymMat(w) * SkewSymMat(w);

intg_Jl_t2 = w * w' * t^2/(2 * w_norm * w_norm) + SkewSymMat(w) * t / (w_norm * w_norm) + (-cos(w_norm * t)/(w_norm * w_norm) + 1/(w_norm * w_norm)) * eye(3) - (-w * w' * cos(w_norm * t) / (w_norm^4) + w * w' / (w_norm^4))...
    - (SkewSymMat(w) * sin(w_norm * t) / (w_norm^3));

intg_Jl_t - intg_Jl_t2
end