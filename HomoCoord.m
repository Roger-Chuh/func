function y = HomoCoord(x,flag)
[m, n] = size(x);
assert(m~=0&&n~=0,'the input matrix should not be [empty]');
if flag  == 1
    one_ = ones(1,n);
    y = [x;one_];
elseif flag == 2
        one_ = ones(m,1);
    y = [x,one_];
else
    assert('the flag should be 1 or 2');
end
end