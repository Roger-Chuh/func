function testImuIntr()

clc

clear


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

if 0
expm(A*dt)
end



a = rand(6, 1);

A = eye(3);
A(1,1) = A(1,1) + a(1);
A(1,2) = A(1,2) + a(2);
A(1,3) = A(1,3) + a(3);
A(2,2) = A(2,2) + a(4);
A(2,3) = A(2,3) + a(5);
A(3,3) = A(3,3) + a(6);


x = rand(3,1);


test1 = A * x;
test2 = x + func_K(x) * a;

test1 - test2


end