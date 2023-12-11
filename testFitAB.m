function testFitAB()

a = 2;
b = 20;

x = [1:10];
y = exp(a).*x - b;
Mat = [x' y'];
Mat(:,3) = 1;
[u,v,s] = svd(Mat' * Mat);
% [u,v,s] = svd(Mat);
para = s(:,3);
para = para./norm(para(1:2));


a_ = -para(1) / para(2);
b_ = para(3) / para(2);





y = [52      37 36.0773 37.5831      37 37.4317 38.1365 69.9194];
x = [40.6281 42.4221 42.6258 42.3077 37.0681 38.0537 37.2001 38.1114];

aa = 0;1.9;
bb = 0;50;
for iter = 1 : 20
    
    H = zeros(2,2);
    b = zeros(2,1);
    Err = 0;
    if 0
        for i = 1 : length(x)
            [err, d_err_d_a, d_err_d_b] = evaluate(x(i), y(i), aa, bb);
            H(1,1) = H(1,1) + d_err_d_a * d_err_d_a;
            H(2,2) = H(2,2) + d_err_d_b * d_err_d_b;
            H(1,2) = H(1,2) + d_err_d_a * d_err_d_b;
            H(2,1) = H(2,1) + d_err_d_b * d_err_d_a;
            b(1,1) = b(1,1) - d_err_d_a*err;
            b(2,1) = b(2,1) - d_err_d_b*err;
            Err = Err + norm(err);
        end
    else
        err = (exp(aa).*x-bb-y)';
        d_err_d_a = exp(aa)*x';
        d_err_d_b = -1 * ones(length(x),1);
        H(1,1) = H(1,1) + d_err_d_a' * d_err_d_a;
        H(2,2) = H(2,2) + d_err_d_b' * d_err_d_b;
        H(1,2) = H(1,2) + d_err_d_a' * d_err_d_b;
        H(2,1) = H(2,1) + d_err_d_b' * d_err_d_a;
        b(1,1) = b(1,1) - d_err_d_a'*err;
        b(2,1) = b(2,1) - d_err_d_b'*err;
        H(1,1) = H(1,1) + 1000;
        H(2,2) = H(2,2) + 1000;
        Err = Err + norm(err);
    end
    dx = inv(H) * b;
    aa = aa + dx(1);
    bb = bb + dx(2);
end


% figure,plot(host', target','or');axis equal
% figure,plot(host', target','or');axis equal;hold on;plot(a*target-b, target','ob');

host = [52 100.163 97.9899 108.528 102.717 14.7207      14      14];
target = [67.85 95.6043  92.807 103.069 94.1911 30.4497 13.0267      11];
Mat = [target' host'];
Mat(:,3) = 1;
[u,v,s] = svd(Mat' * Mat);
para = s(:,3);
para = para./norm(para(1:2));
a_ = -para(1) / para(2);
b_ = para(3) / para(2);
a = exp(0.469126);
norm(a_.*target-b_-host)
norm(target-host)
a1 = exp(0.090167);
b1 = 6.47648;
norm(a1.*target-b1-host)

figure,plot(host', target','or');axis equal;hold on;plot(a1*target-b1, target','ob');
end
function [err, d_err_d_a, d_err_d_b] = evaluate(x, y, a, b)

err = exp(a)*x-b-y;
d_err_d_a = x * exp(a);
d_err_d_b = -1;


end