function [rectParamL, rectParamR] = GetRectifyParam(stereoParam, imgSize)

[rotMatLeft, rotMatRight, intrMatLeftNew, intrMatRightNew, ~, ~] = GetRectifyPlane(stereoParam, imgSize);

focLeft = stereoParam.focLeft;
focRight = stereoParam.focRight;
cenLeft = stereoParam.cenLeft;
cenRight = stereoParam.cenRight;
alphaLeft = stereoParam.alphaLeft;
alphaRight = stereoParam.alphaRight;
kcLeft = stereoParam.kcLeft;
kcRight = stereoParam.kcRight;

% Pre-compute the necessary indices and blending coefficients to enable quick rectification:
if 0
    [~,indNewL,ind1L,ind2L,ind3L,ind4L,a1L,a2L,a3L,a4L] = rect_index(zeros(imgSize(1,1:2)),rotMatLeft,focLeft,cenLeft,kcLeft,alphaLeft,intrMatLeftNew);
    [~,indNewR,ind1R,ind2R,ind3R,ind4R,a1R,a2R,a3R,a4R] = rect_index(zeros(imgSize(1,1:2)),rotMatRight,focRight,cenRight,kcRight,alphaRight,intrMatRightNew);
else
    % 20211221
    [maskL,indNewL,ind1L,ind2L,ind3L,ind4L,a1L,a2L,a3L,a4L] = rect_index(zeros(imgSize(1,1:2)),rotMatLeft,focLeft,cenLeft,kcLeft,alphaLeft,intrMatLeftNew);
    [maskR,indNewR,ind1R,ind2R,ind3R,ind4R,a1R,a2R,a3R,a4R] = rect_index(zeros(imgSize(1,1:2)),rotMatRight,focRight,cenRight,kcRight,alphaRight,intrMatRightNew);
end

rectParamL = struct('rectIntrMat', intrMatLeftNew,  'dstInd', indNewL, 'srcInd', cat(3, ind1L,ind2L,ind3L,ind4L), 'coef', cat(3, a1L,a2L,a3L,a4L));
rectParamR = struct('rectIntrMat', intrMatRightNew, 'dstInd', indNewR, 'srcInd', cat(3, ind1R,ind2R,ind3R,ind4R), 'coef', cat(3, a1R,a2R,a3R,a4R));

end

