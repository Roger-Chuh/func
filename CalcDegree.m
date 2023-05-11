function theta = CalcDegree(vec1, vec2)
if 0 
    theta = min([rad2deg(acos(dot(vec1, vec2)/norm(vec1)/norm(vec2))); rad2deg(acos(dot(vec1, -vec2)/norm(vec1)/norm(vec2)))]);
else
    
    [vec1_norm, ~] = NormalizeVector(vec1);
    [vec2_norm, ~] = NormalizeVector(vec2);
    rad2deg(acos(dot(vec1_norm', vec2_norm')'));
    theta = min([rad2deg(acos(dot(vec1_norm', vec2_norm')')) rad2deg(acos(dot(vec1_norm', -vec2_norm')'))]')';
    
end
end