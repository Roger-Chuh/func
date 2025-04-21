function Rw0 = GetRotationFromV1ToV2(v1, v2)
a = v1 ./norm(v1);
b = v2 ./norm(v2);
n = cross(a, b);
n = n ./norm(n);
cos_theta = a' * b;
theta = 0;

if (cos_theta >= 1)
    theta = 0;
elseif (cos_theta <= -1)
    theta = pi;
else
    theta = acos(cos_theta);
end

Rw0 = rodrigues(theta * n);
end