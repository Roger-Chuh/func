function test2dofGravity()
rot_vec = rand(3, 1);
rot_vec = rot_vec./norm(rot_vec);
% rot_vec(1) = 0;
R = rodrigues(rot_vec);





normal = cross(R(:,3), [1;0;0]);
normal = normal./norm(normal);

dx = rand(2,1);
dx2 = rand(3,1);
dx2(3) = 0;
rot_vec2 = rot_vec + dx2;
Rnew = rodrigues(rot_vec2);
R_new = rodrigues(normal * dx(1)) * R * rodrigues([0;0;1] * dx(2));
% R_new = rodrigues([0;0;1] * dx(2)) * R * rodrigues(normal * dx(1));

normal_new = cross(R_new(:,3), [1;0;0]);
normal_new = normal_new./norm(normal_new);

normal - normal_new


dx3 = rand(3,1);
dx3(1) = 0;




R_align = GetRotationFromV1ToV2(R(:,3), [0; 0; 1]);
R_aligned = R_align * R;
rot_aligned = rodrigues(R_aligned);
rot_aligned_inc = rot_aligned + dx2;
R_aligned_inc = rodrigues(dx2) * R_aligned;

end