function [A] = ProduceOtherOthogonalBasis(n)
N = n;
if (N(0+1) < 0)
    N(0+1) = -N(0+1);
end
if (N(1+1) < 0)
    N(1+1) = -N(1+1);
end
if (N(2+1) < 0)
    N(2+1) = -N(2+1);
end
minIdx = 0+1;
if (N(0+1) <= N(1+1))
    if (N(0+1) <= N(2+1))
        minIdx = 0+1;
    else
        minIdx = 2+1;
    end
else
    if (N(1+1) <= N(2+1))
        minIdx = 1+1;
    else
        minIdx = 2+1;
    end
end
A = zeros(3,2);
switch (minIdx)
    case 1
        A(:,1) = [0, -n(2+1), n(1+1)]';
        
    case 2
        A(:,1) = [n(2+1), 0, -n(0+1)]';
        
    case 3
        A(:,1) = -[-n(1+1), n(0+1), 0]';
end
A(:,2) = cross(n,A(:,1));


A(:,1) = A(:,1)./norm(A(:,1));
A(:,2) = A(:,2)./norm(A(:,2));


end