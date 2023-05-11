function [rectParamL, rectParamR, rotMatLeft, rotMatRight, intrMatLeftNew, intrMatRightNew] = GetRectifyParam2(stereoParam, imgSize)
global cfg

if 1% isempty(cfg.KK_newL)
    [rotMatLeft, rotMatRight, intrMatLeftNew, intrMatRightNew, ~, ~] = GetRectifyPlane(stereoParam, imgSize);
else
    [rotMatLeft, rotMatRight, intrMatLeftNew, intrMatRightNew, ~, ~] = GetRectifyPlane(stereoParam, imgSize);
end


focLeft = stereoParam.focLeft;
focRight = stereoParam.focRight;
cenLeft = stereoParam.cenLeft;
cenRight = stereoParam.cenRight;
alphaLeft = stereoParam.alphaLeft;
alphaRight = stereoParam.alphaRight;
kcLeft = stereoParam.kcLeft;
kcRight = stereoParam.kcRight;

% Pre-compute the necessary indices and blending coefficients to enable quick rectification:
[~,indNewL,ind1L,ind2L,ind3L,ind4L,a1L,a2L,a3L,a4L] = rect_index(zeros(imgSize(1,1:2)),rotMatLeft,focLeft,cenLeft,kcLeft,alphaLeft,intrMatLeftNew);
[~,indNewR,ind1R,ind2R,ind3R,ind4R,a1R,a2R,a3R,a4R] = rect_index(zeros(imgSize(1,1:2)),rotMatRight,focRight,cenRight,kcRight,alphaRight,intrMatRightNew);

maskL = zeros((imgSize(1:2)));
maskL(indNewL) = 1;
maskR = zeros((imgSize(1:2)));
maskR(indNewR) = 1;
maskL(1,:) = 0;
maskR(end,:) = 0;
if kcLeft(1) ~= 0
    if 1 % abs(  length(indNewL) - imgSize(1)*imgSize(2)  ) > 100 && abs(  length(indNewR) - imgSize(1)*imgSize(2)  ) > 100
        try
            %% 20210601  break on purpose
            if 0
                asgjkb;
            end
            if 1
                crit = [1 1 0];
                minSize = [floor(imgSize(1)*.8) floor(imgSize(2)*.8)];%ny,nx
                
                [Xleft, Yleft, Wleft, Hleft] = FindLargestRectanglesPUS(maskL,crit,minSize);
                %Right image: Find borders of largest rectangle with valid image pixels in right image
                
                [Xright, Yright, Wright, Hright] = FindLargestRectanglesPUS(maskR,crit,minSize);
                %Calulate largest common borders
                YbothTop=max(Yleft,Yright);%common most top row
                YbothBottom=min(Yleft+Hleft,Yright+Hright);% common most bottom row
                Hboth=YbothBottom-YbothTop+1;%common height
                titleTxt='common height';
                %Uncomment following 8 lines if left and right images shall have same size and position
                XbothLeft=max(Xleft,Xright);%common most left column
                XbothRight=min(Xleft+Wleft-1,Xright+Wright-1);%common most right column
                Wboth=XbothRight-XbothLeft+1;%Wboth=min(Wleft,Wright);
                Xleft=XbothLeft;
                Xright=XbothLeft;
                Wleft=Wboth;
                Wright=Wboth;
                titleTxt='common rectangle and position';
                
                if 0
                    figure,
                    subplot(1,2,1);
                    imshow(maskL);
                    title(sprintf(['left image calibrated and rectified (full precision).\nred: max rectangle, green: ' titleTxt]));
                    rectangle('Position',[Xleft, Yleft, Wleft, Hleft], 'EdgeColor','r', 'LineWidth',2);%max rectangle left
                    rectangle('Position',[Xleft, YbothTop, Wleft, Hboth], 'EdgeColor','g', 'LineWidth',1);%common height
                    subplot(1,2,2);
                    imshow(maskR);
                    title(sprintf(['right image calibrated and rectified (full precision).\nred: max rectangle, green: ' titleTxt]));
                    rectangle('Position',[Xright, Yright, Wright, Hright], 'EdgeColor','r', 'LineWidth',1);%max rectangle right
                    rectangle('Position',[Xright, YbothTop, Wright, Hboth], 'EdgeColor','g', 'LineWidth',1);%common height
                    drawnow;
                end
                
                xfinal1 = XbothLeft; xfinal2 = XbothLeft +  Wboth;
                yfinal1 = YbothTop; yfinal2 = YbothTop +  Hboth;
                
                
                angX1 = atand(imgSize(2)/2/intrMatLeftNew(1));
                angX2 = atand(Wboth/2/intrMatLeftNew(1));
                fnew1 = imgSize(2)/2/tand(angX2);
                
                angY1 = atand(imgSize(1)/2/intrMatLeftNew(2,2));
                angY2 = atand(Hboth/2/intrMatLeftNew(2,2));
                fnew2 = imgSize(1)/2/tand(angY2);
            else
                
                sumHoriL = sum(maskL');
                sumVertiL = sum(maskL);
                maskSizeL = [sumVertiL(imgSize(2)/2) sumHoriL(imgSize(1)/2)];
                
                sumHoriR = sum(maskR');
                sumVertiR = sum(maskR);
                maskSizeR = [sumVertiR(imgSize(2)/2) sumHoriR(imgSize(1)/2)];
                
                
                maskLR = maskL+maskR == 2;
                
                sumHoriLR = sum(maskLR');
                sumVertiLR = sum(maskLR);
                
                range1 = sort(round(imgSize(2)/2-imgSize(2)/5:imgSize(2)/2+imgSize(2)/5));
                range2 = sort(round(imgSize(1)/2-imgSize(1)/5:imgSize(1)/2+imgSize(1)/5));
                maskSizeLR = [min(sumVertiLR(range1)) min(sumHoriLR(range2))];
                
                
                
                if 0
                    Wboth = min([maskSizeL(2) maskSizeR(2)]);
                    Hboth = min([maskSizeL(1) maskSizeR(1)]);
                else
                    Wboth = maskSizeLR(2);
                    Hboth = maskSizeLR(1);
                end
                
                angX1 = atand(imgSize(2)/2/intrMatLeftNew(1));
                angX2 = atand(Wboth/2/intrMatLeftNew(1));
                fnew1 = imgSize(2)/2/tand(angX2);
                
                angY1 = atand(imgSize(1)/2/intrMatLeftNew(2,2));
                angY2 = atand(Hboth/2/intrMatLeftNew(2,2));
                fnew2 = imgSize(1)/2/tand(angY2);
                
            end
            if 0
                fnew = max(fnew1, fnew2/cfg.fxy_ratio) + 0;  % 50;
            else
                fnew = max(fnew1, fnew2/cfg.fxy_ratio) + 0; % cfg.ext_focal;  % 50;
            end
            
            [rotMatLeft, rotMatRight, intrMatLeftNew, intrMatRightNew, ~, ~] = GetRectifyPlane(stereoParam, imgSize, fnew);
            
            focLeft = stereoParam.focLeft;
            focRight = stereoParam.focRight;
            cenLeft = stereoParam.cenLeft;
            cenRight = stereoParam.cenRight;
            alphaLeft = stereoParam.alphaLeft;
            alphaRight = stereoParam.alphaRight;
            kcLeft = stereoParam.kcLeft;
            kcRight = stereoParam.kcRight;
            
            % Pre-compute the necessary indices and blending coefficients to enable quick rectification:
            [~,indNewL,ind1L,ind2L,ind3L,ind4L,a1L,a2L,a3L,a4L] = rect_index(zeros(imgSize(1,1:2)),rotMatLeft,focLeft,cenLeft,kcLeft,alphaLeft,intrMatLeftNew);
            [~,indNewR,ind1R,ind2R,ind3R,ind4R,a1R,a2R,a3R,a4R] = rect_index(zeros(imgSize(1,1:2)),rotMatRight,focRight,cenRight,kcRight,alphaRight,intrMatRightNew);
            
            maskL2 = zeros((imgSize(1:2)));
            maskL2(indNewL) = 1;
            maskR2 = zeros((imgSize(1:2)));
            maskR2(indNewR) = 1;
            
            
            
        catch
            asdgfk = 1;
        end
    end
    
end

figure,subplot(2,1,1);imshow([maskL maskR]);
subplot(2,1,2);imshow([maskL2 maskR2]);
drawnow;

rectParamL = struct('rectIntrMat', intrMatLeftNew,  'dstInd', indNewL, 'srcInd', cat(3, ind1L,ind2L,ind3L,ind4L), 'coef', cat(3, a1L,a2L,a3L,a4L));
rectParamR = struct('rectIntrMat', intrMatRightNew, 'dstInd', indNewR, 'srcInd', cat(3, ind1R,ind2R,ind3R,ind4R), 'coef', cat(3, a1R,a2R,a3R,a4R));

end

