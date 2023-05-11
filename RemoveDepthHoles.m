function depthNew = RemoveDepthHoles(disparityMap, areaRatio)
% global areaRatio
[L, LL] = bwlabel(disparityMap > 0, 8);

% areaRatio = 0.005; 0.01; 0.01;

for i = 1 : LL
    bwTemp =  L == (i);
    area(i) = sum(sum(bwTemp));
end


validLabel  = find(area > areaRatio*max(area));

validArea = false(size(disparityMap));
for i = 1 : length(validLabel)
    validTemp =  find(L == validLabel(i));
    validArea(validTemp) = true;
end
depthNew = nan(size(disparityMap));

depthNew(validArea) = disparityMap(validArea);
% validMap(~validArea) = false;
% disparityMap(~validArea) = 0;
end
