function [normVec, sqrtsumvec2] = NormalizeVector(vec)
if 1
    vec2 = vec.^2;
    sumvec2 = sum(vec2,2);
    sqrtsumvec2 = sqrt(sumvec2);
    if 0
        for j = 1:size(npixels,1)
            ndescrppp(j,:) = npixels(j,:)/norm(npixels(j,:));
        end
    else
        for k = 1 : size(vec,2)
            normVec(:,k) = vec(:,k)./sqrtsumvec2;
        end
    end
else
    [a] = normc(vec')';
    b=a;
    a(a == 0) = 1;
    scal = (vec./a);
    scal
    id = ~isnan(scal);
    idd = sum(id,1) == size(scal,1);   % idd = find(sum(id) == size(scal,1));
    scall = mean(scal(:,idd),2);
    
    sqrtsumvec2 = scall;
    normVec = b;
end
end