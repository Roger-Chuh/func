function imgRect = RectImg(xOrig2Rect, yOrig2Rect, xRect2Orig, yRect2Orig, image,KK_new, Kinit,k,R,imgOrig,varargin)



% % xOrig2Rect = xOrig2Rect(:,1:size(xRect2Orig,2)/4);
% % yOrig2Rect = yOrig2Rect(:,1:size(yRect2Orig,2)/4);


if (nargin == 10)
    paraDir = [];
    %     homoDir = 'D:\Temp\20180828\calibFishEyeHomo1';
    %     heightDir = 'D:\Temp\20181213\apriltag7';
    
elseif (nargin == 11)
    paraDir = varargin{1};
    whichOne = 'Left';
    scale = 1;
elseif (nargin == 12)
    paraDir = varargin{1};
    whichOne = varargin{2};
    scale = 1;
elseif (nargin == 13)
    paraDir = varargin{1};
    whichOne = varargin{2};
    scale = varargin{3};
else
    error('Too many input arguments');
end

if ~isempty(paraDir)
    fishEyeModel = 1;
    
    if strcmp(whichOne, 'Left')
        load(fullfile(paraDir,'calib.mat'),'oCamModelL','stereoParam')
        try
            oCamModel = oCamModelL;
        catch
            load(fullfile(paraDir,'oCamModelR.mat'));
            oCamModel = oCamModelR;
            
        end
        
        U_same = ocam_undistort_map(oCamModel,'OutputView','same');
        nim_same = ocam_undistort(imgOrig,U_same);
        try
            [rectParamL, rectParamR] = GetRectifyParam(stereoParam, size(imgOrig));
            [rectImgL] = RectifyOneImage(nim_same, rectParamL);
        catch
            asdgljkk = 1;
        end
    else
        load(fullfile(paraDir,'calib.mat'),'oCamModelR','stereoParam')
        oCamModel = oCamModelR;
        [rectParamL, rectParamR] = GetRectifyParam(stereoParam, size(imgOrig));
        U_same = ocam_undistort_map(oCamModel,'OutputView','same');
        nim_same = ocam_undistort(imgOrig,U_same);
        [rectImgR] = RectifyOneImage(nim_same, rectParamR);
    end
    
else
    fishEyeModel = 0;
end

imgRect = uint8(zeros(size(imgOrig)));



useTable = 1;0;1;0;1;0;1;1;
draw = 0;1;0;1;




useNearest = 0;



[nr,nc,~] = size(image);
nr = scale*nr;
nc = scale*nc;
numPix = nr*nc;
[xuuAllUse,yuuAllUse] = meshgrid(1:nc,1:nr);
xNum = 300;
SearchTimes = [];
SearchY = [];
Err_table = [];
PixRectStack = [];
PixRectRStack = [];
Pix2Stack = [];
Pix2StackIn = [];
step = 1; 3;
yStep = scale; 1;2;1;2;20; 3;15;3;10; 5; 10;
splitStep = 2; 3; 5;3; 5;
flag = false;
ySampled = [scale:yStep:nr]';
if ySampled(end) ~= nr
    ySampled = [ySampled;nr];
end
xAxis = [-100:nc+100];

% figure(111),imshow(zeros(nr,nc));hold on;
plotRow1 = 0;
judge = 1;
ErrDir = []; ErrWrong = [];
vldMsk = zeros(nr,nc);
pixStack = [];
figure(444),imshow(zeros(nr,nc));hold on;
figure(333),imshow(ones(nr,nc));hold on;
figure(555),imshow(ones(nr,nc));hold on;

if useNearest == 0
    figNum = 666;
else
    figNum = 777;
end
figure(figNum),hold on;
pixWrite = [];
StoreData = []; asgkdfghdfghu = 0;
Ct = {};
ct = 1;
ccccctttt = 0;
ZeroFlag = [];

Pix1 = [];
Pix2 = [];
Pix3 = [];
Pix4 = [];
Pix5 = [];

if 0
    crop1 = [1 1];
    crop2 = [nc-0 nr-0];
elseif 0
    crop1 = [146 209];
    crop1 = [400-256 300+256-16];
    crop2 = [3752 1899+256];
elseif 0
    crop1 = [11 9];
    %     crop1 = [400-256 300+256-16];
    crop2 = [1834 1048];
    
elseif 0
    crop1 = [23 15];
    %     crop1 = [400-256 300+256-16];
    crop2 = [3750 2126];
elseif 0
    crop1 = [23 15];
    crop2 = [3750-2048-32 2126-1024-128];
elseif 0 % temp
    crop1 = [23 15];
    crop2 = [3750 2126-1024-512-128];
elseif 0
    crop1 = [23 15];
    crop2 = [3750-2048 2126-1024-128];
elseif 1
    crop1 = [23 15];
    crop2 = [582 414+64];
else
    crop1 = [17 3];
    crop2 = [1136+128 642+64];
end
crop1 = [17 3];
crop2 = [1264 706];

imgR = image(:,:,1);
imgG = image(:,:,2);
imgB = image(:,:,3);

debugLine = -719;
% for ii = scale : scale:length(ySampled) - 1
for ii = scale : scale: scale*(size(image,1) - 1)
    % for ii = scale : scale: nr - scale
    i = ii/scale;
    if  i <  debugLine %0 %
        continue;
        %         draw = 1;
    end
    bufferedX = xuuAllUse(ySampled(i):ySampled(i+1),:);
    bufferedY = yuuAllUse(ySampled(i):ySampled(i+1),:);
    if useTable
        pixRect_1 = [xOrig2Rect(ySampled(i):ySampled(i)+scale-1,:)' yOrig2Rect(ySampled(i):ySampled(i)+scale-1,:)'];
        pixRect_1_tmp = zeros(size(pixRect_1,1)*scale,2);
        pixRect_1_tmp(1:scale:end,:) = [pixRect_1(:,[1 1+scale])];
        if scale == 2
            pixRect_1_tmp(2:scale:end,:) = [pixRect_1(:,[2 2+scale])];
        end
        if 0
            pixRect_1 = pixRect_1(:);
            pixRect_1 = reshape(pixRect_1,[],2);
        else
            pixRect_1 = pixRect_1_tmp;
        end
        pixRect_111 = pixRect_1;
        pixRect_1_ = Orig2Rect([xuuAllUse(ySampled(i),:); yuuAllUse(ySampled(i),:)]', Kinit, KK_new, R,k);
        if 0
            figure,imshow(zeros(nr,nc));hold on;plot(pixRect_1_(:,1),pixRect_1_(:,2),'or');plot(pixRect_1(:,1),pixRect_1(:,2),'.g');
        end
    else
        if fishEyeModel == 0
            pixRect_1 = Orig2Rect([xuuAllUse(ySampled(i),:); yuuAllUse(ySampled(i),:)]', Kinit, KK_new, R,k);
        else
            pixRect_1 = Orig2RectFishEye([xuuAllUse(ySampled(i),:); yuuAllUse(ySampled(i),:)]', Kinit, KK_new, R,k,oCamModel);
        end
    end
    if 0
        pixRect_1 = pixRect_1(scale:scale:end,:);
    end
    
    %     pixRect_1(:,1) = round(pixRect_1(:,1));
    %     pixRect_1(:,1) = floor(pixRect_1(:,1));
    pixRect_1(:,1) = ceil(pixRect_1(:,1));
    pixRect_1(:,2) = ceil(pixRect_1(:,2));
    [aa1,bb1,cc1] = unique(pixRect_1,'rows');
    cc1 = sort(cc1,'ascend');
    ccd1 = diff(cc1);
    [indexCell1, ~] = splitIndex2(find(ccd1 == 1));
    iccd1 = cell2mat(indexCell1')';
    pixOrig_1 = [xuuAllUse(ySampled(i),:); yuuAllUse(ySampled(i),:)]';
    %%    pixOrig_1 = pixOrig_1(iccd1,:);
    pixRect_1010 = pixRect_1;
    %%    pixRect_1 = pixRect_1(iccd1,:);
    pixRect_1_00 = pixRect_1;
    flag1 = pixRect_1(:,1) >=1 & pixRect_1(:,1)<=nc & pixRect_1(:,2) >=1 & pixRect_1(:,2) <= nr;
    if scale == 2
        %         flag1 = pixRect_1(:,1) >=max(crop1(1),2) & pixRect_1(:,1)<=crop2(1) & pixRect_1(:,2) >=max(2,crop1(2)) & pixRect_1(:,2) <= crop2(2);
        %         flag1 = pixRect_1(:,1) >=max(crop1(1),2) & pixRect_1(:,1)<=crop2(1) & pixRect_1(:,2) >=max(2,crop1(2)) & pixRect_1(:,2) <= crop2(2);
        flag1 = pixRect_1(:,1) >=max(crop1(1)+1) & pixRect_1(:,1)<=crop2(1) & pixRect_1(:,2) >=max(1+crop1(2)) & pixRect_1(:,2) <= crop2(2);
    else
        flag1 = pixRect_1(:,1) >=max(crop1(1),-2) & pixRect_1(:,1)<=crop2(1) & pixRect_1(:,2) >=max(-2,crop1(2)) & pixRect_1(:,2) <= crop2(2);
    end
    pixRect_1 = pixRect_1(flag1,:);
    pixOrig_1 = pixOrig_1(flag1,:);
    
    
    %     idInScale = find(mod(pixRect_1(:,1),scale) == 0 & mod(pixRect_1(:,2),scale) == 0);
    
    if ~isempty(pixRect_1)
        if mod(pixRect_1(1,1),scale) ~= 0
            %             pixRect_1(1,1) = max(2,pixRect_1(1,1) -1);
            %% pixRect_1(1,1) = max(2,pixRect_1(1,1) + 1);
            %             pixRect_1(1,1) = max(2,pixRect_1(1,1) - 1); % 20190625 for wei's convenience
            pixRect_1(1,1) = pixRect_1(1,1) - 1;
        end
        if mod(pixRect_1(1,2),scale) ~= 0
            %             pixRect_1(1,1) = max(2,pixRect_1(1,1) -1);
            %             pixRect_1(1,2) = max(2,pixRect_1(1,2) - 1);
            pixRect_1(1,2) = pixRect_1(1,2) - 1;
        end
    else
        kjhbasv =1;
    end
    
    
    %      if scale == 2
    %         flag11 = pixRect_1(:,1) >=max(crop1(1),2) & pixRect_1(:,1)<=crop2(1) & pixRect_1(:,2) >=max(2,crop1(2)) & pixRect_1(:,2) <= crop2(2);
    %     else
    %         flag11 = pixRect_1(:,1) >=max(crop1(1),-2) & pixRect_1(:,1)<=crop2(1) & pixRect_1(:,2) >=max(-2,crop1(2)) & pixRect_1(:,2) <= crop2(2);
    %     end
    %     pixRect_1 = pixRect_1(flag11,:);
    %     pixOrig_1 = pixOrig_1(flag11,:);
    
    %     imgR = image(:,:,1);
    %     imgG = image(:,:,2);
    %     imgB = image(:,:,3);
    
    
    
    if ~isempty(pixRect_1)
        startPixRect = pixRect_1(1,:);
        if draw
            figure(999),clf;imshow(zeros(nr,nc));hold on;plot(pixRect_1(:,1),pixRect_1(:,2),'.g');
        end
        
        
        pix = startPixRect;
        
        
        pixFound = []; pixHoldup1 = [];pixHolddown1 = []; % pixWrite = [];
        pixHoldup2 = [];pixHolddown2 = [];
        FlagTmp = {};
        %             for j =  startPixRect(1) : nc  % 1 : nc % endPixRect(1)
        for j =  startPixRect(1) : scale: crop2(1)
            if max(abs(diff(pixRect_1(:,2))))> 1
                askb = 1;
            end
            if j == 120
                dfoghl = 1;
            end
            pixTmp = [repmat(j,5,1) [pix(2)-2*scale pix(2)-1*scale pix(2) pix(2)+1*scale pix(2)+ 2*scale]'];
            %                 if pixTmp(3,2) <= crop1(2) || pixTmp(3,2) >= crop2(2)+scale-1
            if 0 % pixTmp(3,2) < crop1(2) || pixTmp(3,2) > crop2(2)
                continue;
            end
            flagIn = (pixTmp(:,1)>0 & pixTmp(:,1)<=nc&pixTmp(:,2)>0 & pixTmp(:,2)<=nr);
            if sum(flagIn) == 0
                sgjk = 1;
            end
            index = sub2ind([nr,nc],pixTmp(flagIn,2),pixTmp(flagIn,1));
            if useTable
                
                pixTmp_addon = [repmat(j,7,1) [pix(2)-3 pix(2)-2 pix(2)-1 pix(2) pix(2)+1 pix(2)+ 2 pix(2)+3]'];
                flagIn_addon = (pixTmp_addon(:,1)>0 & pixTmp_addon(:,1)<=nc & pixTmp_addon(:,2)>0 & pixTmp_addon(:,2)<=nr);
                index_addon = sub2ind([nr,nc],pixTmp_addon(flagIn_addon,2),pixTmp_addon(flagIn_addon,1));
                aaa = [xRect2Orig(index_addon) yRect2Orig(index_addon)];
                
                Rect2OrigTmpIn = [xRect2Orig(index) yRect2Orig(index)];
                Rect2OrigTmp = pixTmp;
                Rect2OrigTmp(~flagIn,:) = repmat([-10 -10],sum(~flagIn),1);
                if sum(~flagIn) > 1
                    asgjkf = 1;
                end
                Rect2OrigTmp(flagIn,:) = Rect2OrigTmpIn;
                Rect2OrigTmp_ = remapRect(pixTmp', KK_new,Kinit,k,R);
                
            else
                if fishEyeModel == 0
                    Rect2OrigTmp = remapRect(pixTmp', KK_new,Kinit,k,R);
                else
                    Rect2OrigTmp = remapRectFishEye(pixTmp', KK_new,Kinit,k,R,oCamModel);
                    if 1
                        %                             aaa = remapRectFishEye([pixTmp' [409;5]], KK_new,Kinit,k,R,oCamModel);
                        aaa = remapRectFishEye([[pixTmp(1,1);pixTmp(1,2)-1] pixTmp' [pixTmp(end,1);pixTmp(end,2)+1]], KK_new,Kinit,k,R,oCamModel);
                    end
                    if ismember([409 5], pixTmp,'rows')
                        asdkj = 1;
                    end
                end
            end
            zeroFlag = zeros(5,1);
            [Rect2OrigTmp__, ~] = Num2Fix_unsigned([Rect2OrigTmp]./scale, 20,9);
            flagTmp = find(floor(Rect2OrigTmp__(:,1)) >= 1& floor(Rect2OrigTmp__(:,1)) <= nc/scale-1 & Rect2OrigTmp(:,2) >= ySampled(i) & Rect2OrigTmp(:,2) < ySampled(i+1) & Rect2OrigTmp(:,2) > 0 & Rect2OrigTmp(:,2) <= nr & pixTmp(:,2) > 0 & pixTmp(:,2) <= nr);
            if nc > 3000   % && scale == 1
                flagTmp = setdiff(flagTmp,[1;5]);
            end
            if ~isempty(flagTmp)
                zeroFlag(flagTmp) = 1;
                ZeroFlag = [ZeroFlag [zeroFlag;j]];
                FlagTmp = [FlagTmp; flagTmp];
                % assume there is no way 'length(flagTmp) == 3' could happen.
                
                if ismember(2+1,flagTmp)
                    pixFound = [pixFound; pixTmp(2+1,:)];
                end
                
                if ismember(2+1,flagTmp) && ((~ismember(1+1,flagTmp) && ~isempty(pixHoldup1)) || (~ismember(0+1,flagTmp) && ~isempty(pixHoldup2)))
                    Pix1 = [Pix1; [pixHoldup2 repmat(size(pixHoldup2,1),size(pixHoldup2,1),1)]];
                    Pix2 = [Pix2; [pixHoldup1 repmat(size(pixHoldup1,1),size(pixHoldup1,1),1)]];
                    pixWrite = [pixWrite;pixHoldup1;pixHoldup2];
                    %                     StoreData = [StoreData; [mean(pixHoldup1,1) size(pixHoldup1,1)]];
                    StoreData = [StoreData; [[pixHoldup1;pixHoldup2] repmat(size([pixHoldup1;pixHoldup2],1),size([pixHoldup1;pixHoldup2],1),1)]];
                    
                    if draw
                        pixHoldupAll = [pixHoldup1; pixHoldup2];
                        figure(333),plot(pixHoldupAll(:,1),pixHoldupAll(:,2),'.','Color',rand(1,3));
                        foundInd = sub2ind([nr/1,nc/1],pixHoldupAll(:,2)./1,pixHoldupAll(:,1)./1);
                        try
                            xyOrig = round([xRect2Orig(foundInd) yRect2Orig(foundInd)]);
                            idIn = xyOrig(:,1) > 0 & xyOrig(:,1) <= nc & xyOrig(:,2) > 0 & xyOrig(:,2) <= nr;
                            [xyOrig_, ~] = Num2Fix_unsigned([xRect2Orig(foundInd) yRect2Orig(foundInd)]./scale, 20,9);
                            xyOrig_int = floor(xyOrig_);
                            xyOrig_frac = xyOrig_ - xyOrig_int;
                            
                            y_in = [ySampled(i);ySampled(i+1)];
                            
                            %                                 flag_in = xyOrig_int(:,1) >= scale & xyOrig_int(:,1) <= nc/scale-1 & xyOrig_int(:,2) >= scale & xyOrig_int(:,2) <= nr/scale-1;
                            flag_in = xyOrig_int(:,1) >= 1 & xyOrig_int(:,1) <= nc/scale-1 & xyOrig_int(:,2) >= 1 & xyOrig_int(:,2) <= nr/scale-1;
                            xyOrig_int_in = xyOrig_int(flag_in,:);
                            xyOrig_frac_in = xyOrig_frac(flag_in,:);
                            if useNearest
                                foundIndIn = foundInd(idIn);
                                origInd = sub2ind([nr,nc],xyOrig(idIn,2),xyOrig(idIn,1));
                                imgRect([foundIndIn;numPix + foundIndIn;2*numPix + foundIndIn]) = imgOrig([origInd;numPix + origInd;2*numPix + origInd]);
                                %                             figure(figNum),clf;imshow(imgRect);
                                
                            else
                                foundIndIn = foundInd(flag_in);
                                [yyy,xxx] = ind2sub([nr,nc],foundIndIn);
                                foundIndIn = sub2ind([nr,nc]/scale, yyy/scale, xxx/scale);
                                imgRect = BilinearRemap(imgRect, imgR, imgG, imgB, image, xyOrig_int_in, xyOrig_frac_in,foundIndIn);
                                
                            end
                        catch
                            dakvhb = 1;
                        end
                    end
                    
                    pixHoldup1 = []; pixHoldup2 = [];
                    
                end
                if ismember(1+1,flagTmp)|| ismember(0+1,flagTmp)
                    if ismember(1+1,flagTmp)
                        pixHoldup1 = [pixHoldup1; pixTmp(1+1,:)];
                    end
                    if ismember(0+1,flagTmp)
                        pixHoldup2 = [pixHoldup2; pixTmp(0+1,:)];
                        %                             aaa = remapRectFishEye([pixTmp' [pixTmp(end,1);pixTmp(end,2)+1]], KK_new,Kinit,k,R,oCamModel);
                        %                             aaa = remapRectFishEye([[pixTmp(1,1);pixTmp(1,2)-1] pixTmp' [pixTmp(end,1);pixTmp(end,2)+1]], KK_new,Kinit,k,R,oCamModel);
                        Ct{ct,1} = aaa;
                        Ct{ct,2} = pixTmp;
                        Ct{ct,3} = i;
                        ct = ct + 1;
                    end
                    %                         pixHoldup1 = [pixHoldup1; pixTmp(1+1,:)];
                    if ~ismember(2+1,flagTmp)
                        Pix3 = [Pix3; [pixFound repmat(size(pixFound,1),size(pixFound,1),1)]];
                        pixWrite = [pixWrite;pixFound];
                        %                         StoreData = [StoreData; [mean(pixFound,1) size(pixFound,1)]];
                        
                        if ~isempty(pixFound)
                            %                             StoreData = [StoreData; [mean(pixFound,1) size(pixFound,1)]];
                            StoreData = [StoreData; [pixFound repmat(size(pixFound,1),size(pixFound,1),1)]];
                            
                            if draw
                                figure(333),plot(pixFound(:,1),pixFound(:,2),'.','Color',rand(1,3));
                                
                                %                                     foundInd = sub2ind([nr,nc],pixFound(:,2),pixFound(:,1));
                                foundInd = sub2ind([nr/1,nc/1],pixFound(:,2)./1,pixFound(:,1)./1);
                                try
                                    xyOrig = round([xRect2Orig(foundInd) yRect2Orig(foundInd)]);
                                    xyOrig_ = [xRect2Orig(foundInd) yRect2Orig(foundInd)];
                                    [xyOrig_, ~] = Num2Fix_unsigned([xRect2Orig(foundInd) yRect2Orig(foundInd)]./scale, 20,9);
                                    
                                    xyOrig_int = floor(xyOrig_);
                                    xyOrig_frac = xyOrig_ - xyOrig_int;
                                    
                                    y_in = [ySampled(i);ySampled(i+1)];
                                    
                                    %                                         flag_in = xyOrig_int(:,1) >= scale & xyOrig_int(:,1) <= nc/scale-1 & xyOrig_int(:,2) >= scale & xyOrig_int(:,2) <= nr/scale-1;
                                    flag_in = xyOrig_int(:,1) >= 1 & xyOrig_int(:,1) <= nc/scale-1 & xyOrig_int(:,2) >= 1 & xyOrig_int(:,2) <= nr/scale-1;
                                    xyOrig_int_in = xyOrig_int(flag_in,:);
                                    xyOrig_frac_in = xyOrig_frac(flag_in,:);
                                    %
                                    %
                                    %                                     PixNew3_Bilinear_floor1 = xyOrig_int_in;
                                    %                                     PixNew3_Bilinear_floor2 = [xyOrig_int_in(:,1)+1 xyOrig_int_in(:,2)];
                                    %                                     PixNew3_Bilinear_floor3 = [xyOrig_int_in(:,1) xyOrig_int_in(:,2)+1];
                                    %                                     PixNew3_Bilinear_floor4 = [xyOrig_int_in(:,1)+1 xyOrig_int_in(:,2)+1];
                                    %
                                    %                                     indBilinear1 = sub2ind(size(image(:,:,1)),PixNew3_Bilinear_floor1(:,2),PixNew3_Bilinear_floor1(:,1));
                                    %                                     indBilinear2 = sub2ind(size(image(:,:,1)),PixNew3_Bilinear_floor2(:,2),PixNew3_Bilinear_floor2(:,1));
                                    %                                     indBilinear3 = sub2ind(size(image(:,:,1)),PixNew3_Bilinear_floor3(:,2),PixNew3_Bilinear_floor3(:,1));
                                    %                                     indBilinear4 = sub2ind(size(image(:,:,1)),PixNew3_Bilinear_floor4(:,2),PixNew3_Bilinear_floor4(:,1));
                                    %
                                    %                                     coeff1 = (1 - xyOrig_frac_in(:,2)).*(1 - xyOrig_frac_in(:,1));
                                    %                                     coeff2 = (1 - xyOrig_frac_in(:,2)).*xyOrig_frac_in(:,1);
                                    %                                     coeff3 = xyOrig_frac_in(:,2).*(1 - xyOrig_frac_in(:,1));
                                    %                                     coeff4 = xyOrig_frac_in(:,2).*xyOrig_frac_in(:,1);
                                    
                                    
                                    
                                    if ~((max(xyOrig_int_in(:,2)) == min(xyOrig_int_in(:,2))) && (min(xyOrig_int_in(:,2)) == y_in(1)))
                                        ccccctttt = ccccctttt + 1;
                                    end
                                    
                                    
                                    
                                    
                                    idIn = xyOrig(:,1) > 0 & xyOrig(:,1) <= nc & xyOrig(:,2) > 0 & xyOrig(:,2) <= nr;
                                    % % %                                     idIn = flag_in;
                                    
                                    if useNearest
                                        foundIndIn = foundInd(idIn);
                                        %                                     foundIndIn = foundInd(flag_in);
                                        origInd = sub2ind([nr,nc],xyOrig(idIn,2),xyOrig(idIn,1));
                                        imgRect([foundIndIn;numPix + foundIndIn;2*numPix + foundIndIn]) = imgOrig([origInd;numPix + origInd;2*numPix + origInd]);
                                        % % %
                                        % %
                                        % %                                     imgRect(foundIndIn) = coeff1.*double(imgR(indBilinear1)) + coeff2.*double(imgR(indBilinear2)) + coeff3.*double(imgR(indBilinear3)) + coeff4.*double(imgR(indBilinear4));
                                        % %                                     imgRect(numPix + foundIndIn) = coeff1.*double(imgG(indBilinear1)) + coeff2.*double(imgG(indBilinear2)) + coeff3.*double(imgG(indBilinear3)) + coeff4.*double(imgG(indBilinear4));
                                        % %                                     imgRect(2*numPix + foundIndIn) = coeff1.*double(imgB(indBilinear1)) + coeff2.*double(imgB(indBilinear2)) + coeff3.*double(imgB(indBilinear3)) + coeff4.*double(imgB(indBilinear4));
                                    else
                                        foundIndIn = foundInd(flag_in);
                                        [yyy,xxx] = ind2sub([nr,nc],foundIndIn);
                                        foundIndIn = sub2ind([nr,nc]/scale, yyy/scale, xxx/scale);
                                        imgRect = BilinearRemap(imgRect, imgR, imgG, imgB, image, xyOrig_int_in, xyOrig_frac_in,foundIndIn);
                                        
                                    end
                                    %                             figure(figNum),clf;imshow(imgRect);
                                catch
                                    skav = 1;
                                end
                            end
                        else
                            svig = 1;
                        end
                        pixFound = pixHoldup1;
                        pixHoldup1 = pixHoldup2;
                        pixHoldup2 = [];
                        pix = pixTmp(1+1,:);
                        
                    end
                    
                end
                
                %                     if ismember(0+1,flagTmp)
                %                         pixHoldup2 = [pixHoldup2; pixTmp(0+1,:)];
                %                         if ~ismember(2+1,flagTmp)
                %                             pixWrite = [pixWrite;pixFound];
                %                             %                         StoreData = [StoreData; [mean(pixFound,1) size(pixFound,1)]];
                %
                %                             if ~isempty(pixFound)
                %                                 %                             StoreData = [StoreData; [mean(pixFound,1) size(pixFound,1)]];
                %                                 StoreData = [StoreData; [pixFound repmat(size(pixFound,1),size(pixFound,1),1)]];
                %
                %                                 if draw
                %                                     figure(333),plot(pixFound(:,1),pixFound(:,2),'.','Color',rand(1,3));
                %
                %                                     foundInd = sub2ind([nr,nc],pixFound(:,2),pixFound(:,1));
                %                                     xyOrig = round([xRect2Orig(foundInd) yRect2Orig(foundInd)]);
                %                                     idIn = xyOrig(:,1) > 0 & xyOrig(:,1) <= nc & xyOrig(:,2) > 0 & xyOrig(:,2) <= nr;
                %                                     foundIndIn = foundInd(idIn);
                %                                     origInd = sub2ind([nr,nc],xyOrig(idIn,2),xyOrig(idIn,1));
                %                                     imgRect([foundIndIn;numPix + foundIndIn;2*numPix + foundIndIn]) = imgOrig([origInd;numPix + origInd;2*numPix + origInd]);
                %                                     %                             figure(figNum),clf;imshow(imgRect);
                %                                 end
                %                             else
                %                                 svig = 1;
                %                             end
                %                             pixFound = pixHoldup2;
                %                             pixHoldup2 = [];
                %                             pix = pixTmp(0+1,:);
                %
                %                         end
                %
                %                     end
                
                
                if ismember(2+1,flagTmp) && ((~ismember(3+1,flagTmp) && ~isempty(pixHolddown1)) || (~ismember(4+1,flagTmp) && ~isempty(pixHolddown2)))
                    
                    Pix4 = [Pix4; [pixHolddown1 repmat(size(pixHolddown1,1),size(pixHolddown1,1),1)]];
                    Pix5 = [Pix5; [pixHolddown2 repmat(size(pixHolddown2,1),size(pixHolddown2,1),1)]];
                    
                    pixWrite = [pixWrite;pixHolddown1;pixHolddown2];
                    StoreData = [StoreData; [[pixHolddown1;pixHolddown2] repmat(size([pixHolddown1;pixHolddown2],1),size([pixHolddown1;pixHolddown2],1),1)]];
                    
                    
                    if draw
                        pixHolddownAll = [pixHolddown1;pixHolddown2];
                        figure(333),plot(pixHolddownAll(:,1),pixHolddownAll(:,2),'.','Color',rand(1,3));
                        %                             foundInd = sub2ind([nr,nc],pixHolddownAll(:,2),pixHolddownAll(:,1));
                        foundInd = sub2ind([nr/1,nc/1],pixHolddownAll(:,2)./1,pixHolddownAll(:,1)./1);
                        try
                            xyOrig = round([xRect2Orig(foundInd) yRect2Orig(foundInd)]);
                            idIn = xyOrig(:,1) > 0 & xyOrig(:,1) <= nc & xyOrig(:,2) > 0 & xyOrig(:,2) <= nr;
                            
                            xyOrig_ = [xRect2Orig(foundInd) yRect2Orig(foundInd)];
                            [xyOrig_, ~] = Num2Fix_unsigned([xRect2Orig(foundInd) yRect2Orig(foundInd)]./scale, 20,9);
                            
                            xyOrig_int = floor(xyOrig_);
                            xyOrig_frac = xyOrig_ - xyOrig_int;
                            
                            y_in = [ySampled(i);ySampled(i+1)];
                            
                            flag_in = xyOrig_int(:,1) >= 1 & xyOrig_int(:,1) <= nc/scale-1 & xyOrig_int(:,2) >= 1 & xyOrig_int(:,2) <= nr/scale-1;
                            
                            xyOrig_int_in = xyOrig_int(flag_in,:);
                            xyOrig_frac_in = xyOrig_frac(flag_in,:);
                            if useNearest
                                foundIndIn = foundInd(idIn);
                                origInd = sub2ind([nr,nc],xyOrig(idIn,2),xyOrig(idIn,1));
                                imgRect([foundIndIn;numPix + foundIndIn;2*numPix + foundIndIn]) = imgOrig([origInd;numPix + origInd;2*numPix + origInd]);
                                %                             figure(figNum),clf;imshow(imgRect);
                            else
                                foundIndIn = foundInd(flag_in);
                                [yyy,xxx] = ind2sub([nr,nc],foundIndIn);
                                foundIndIn = sub2ind([nr,nc]/scale, yyy/scale, xxx/scale);
                                imgRect = BilinearRemap(imgRect, imgR, imgG, imgB, image, xyOrig_int_in, xyOrig_frac_in,foundIndIn);
                                
                            end
                        catch
                            sakvhj = 1;
                        end
                    end
                    
                    pixHolddown1 = [];
                    pixHolddown2 = [];
                end
                
                if ismember(3+1,flagTmp) || ismember(4+1,flagTmp)
                    if ismember(3+1,flagTmp)
                        pixHolddown1 = [pixHolddown1; pixTmp(3+1,:)];
                    end
                    if ismember(4+1,flagTmp)
                        pixHolddown2 = [pixHolddown2; pixTmp(4+1,:)];
                        %                             aaa = remapRectFishEye([[pixTmp(1,1);pixTmp(1,2)-1] pixTmp' [pixTmp(end,1);pixTmp(end,2)+1]], KK_new,Kinit,k,R,oCamModel);
                        Ct{ct,1} = aaa;
                        Ct{ct,2} = pixTmp;
                        Ct{ct,3} = i;
                        ct = ct + 1;
                    end
                    if ~ismember(2+1,flagTmp)
                        Pix3 = [Pix3; [pixFound repmat(size(pixFound,1),size(pixFound,1),1)]];
                        pixWrite = [pixWrite;pixFound];
                        %                         try
                        %                             StoreData = [StoreData; [mean(pixFound,1) size(pixFound,1)]];
                        %                         catch
                        %                             dsgkh = 1;
                        %                         end
                        
                        if ~isempty(pixFound)
                            StoreData = [StoreData; [pixFound repmat(size(pixFound,1),size(pixFound,1),1)]];
                            if draw
                                figure(333),plot(pixFound(:,1),pixFound(:,2),'.','Color',rand(1,3));
                                %                                     foundInd = sub2ind([nr,nc],pixFound(:,2),pixFound(:,1));
                                foundInd = sub2ind([nr/1,nc/1],pixFound(:,2)./1,pixFound(:,1)./1);
                                try
                                    xyOrig = round([xRect2Orig(foundInd) yRect2Orig(foundInd)]);
                                    idIn = xyOrig(:,1) > 0 & xyOrig(:,1) <= nc & xyOrig(:,2) > 0 & xyOrig(:,2) <= nr;
                                    
                                    xyOrig_ = [xRect2Orig(foundInd) yRect2Orig(foundInd)];
                                    [xyOrig_, ~] = Num2Fix_unsigned([xRect2Orig(foundInd) yRect2Orig(foundInd)]./scale, 20,9);
                                    
                                    xyOrig_int = floor(xyOrig_);
                                    xyOrig_frac = xyOrig_ - xyOrig_int;
                                    
                                    y_in = [ySampled(i);ySampled(i+1)];
                                    
                                    %                                         flag_in = xyOrig_int(:,1) >= scale & xyOrig_int(:,1) <= nc/scale-1 & xyOrig_int(:,2) >= scale & xyOrig_int(:,2) <= nr/scale-1;
                                    flag_in = xyOrig_int(:,1) >= 1 & xyOrig_int(:,1) <= nc/scale-1 & xyOrig_int(:,2) >= 1 & xyOrig_int(:,2) <= nr/scale-1;
                                    xyOrig_int_in = xyOrig_int(flag_in,:);
                                    xyOrig_frac_in = xyOrig_frac(flag_in,:);
                                    
                                    if useNearest
                                        foundIndIn = foundInd(idIn);
                                        origInd = sub2ind([nr,nc],xyOrig(idIn,2),xyOrig(idIn,1));
                                        imgRect([foundIndIn;numPix + foundIndIn;2*numPix + foundIndIn]) = imgOrig([origInd;numPix + origInd;2*numPix + origInd]);
                                        %                             figure(figNum),clf;imshow(imgRect);
                                    else
                                        foundIndIn = foundInd(flag_in);
                                        [yyy,xxx] = ind2sub([nr,nc],foundIndIn);
                                        foundIndIn = sub2ind([nr,nc]/scale, yyy/scale, xxx/scale);
                                        imgRect = BilinearRemap(imgRect, imgR, imgG, imgB, image, xyOrig_int_in, xyOrig_frac_in,foundIndIn);
                                        
                                    end
                                catch
                                    eqghk = 1;
                                end
                            end
                        else
                            savihu = 1;
                        end
                        pixFound = pixHolddown1;
                        pixHolddown1 = pixHolddown2;
                        pixHolddown2 = [];
                        pix = pixTmp(3+1,:);
                    end
                end
                
                
                
            else
                sigfg = 1;
            end
        end
        
        
        asdjkbh = 1;
        
        try
            pixLeft = [pixHoldup;pixHolddown];
        catch
            pixLeft = [pixHoldup1;pixHoldup2;pixHolddown1;pixHolddown2];
        end
        Pix1 = [Pix1; [pixHoldup2 repmat(size(pixHoldup2,1),size(pixHoldup2,1),1)]];
        Pix2 = [Pix2; [pixHoldup1 repmat(size(pixHoldup1,1),size(pixHoldup1,1),1)]];
        Pix3 = [Pix3; [pixFound repmat(size(pixFound,1),size(pixFound,1),1)]];
        Pix4 = [Pix4; [pixHolddown1 repmat(size(pixHolddown1,1),size(pixHolddown1,1),1)]];
        Pix5 = [Pix5; [pixHolddown2 repmat(size(pixHolddown2,1),size(pixHolddown2,1),1)]];
        pixWrite = [pixWrite;pixFound;pixLeft];
        try
            StoreData = [StoreData; [[pixFound;pixLeft] repmat(size([pixFound;pixLeft],1),size([pixFound;pixLeft],1),1)]];
        catch
            asgkdfghdfghu = asgkdfghdfghu + 1;
        end
        
        if ~isempty(pixFound)
            if draw
                figure(333),plot(pixFound(:,1),pixFound(:,2),'.','Color',rand(1,3));
                %                 foundInd = sub2ind([nr,nc],pixFound(:,2),pixFound(:,1));
                foundInd = sub2ind([nr/1,nc/1],pixFound(:,2)./1,pixFound(:,1)./1);
                try
                    xyOrig = round([xRect2Orig(foundInd) yRect2Orig(foundInd)]);
                    idIn = xyOrig(:,1) > 0 & xyOrig(:,1) <= nc & xyOrig(:,2) > 0 & xyOrig(:,2) <= nr;
                    
                    xyOrig_ = [xRect2Orig(foundInd) yRect2Orig(foundInd)];
                    [xyOrig_, ~] = Num2Fix_unsigned([xRect2Orig(foundInd) yRect2Orig(foundInd)]./scale, 20,9);
                    
                    xyOrig_int = floor(xyOrig_);
                    xyOrig_frac = xyOrig_ - xyOrig_int;
                    
                    y_in = [ySampled(i);ySampled(i+1)];
                    
                    %                     flag_in = xyOrig_int(:,1) >= scale & xyOrig_int(:,1) <= nc/scale-1 & xyOrig_int(:,2) >= scale & xyOrig_int(:,2) <= nr/scale-1;
                    flag_in = xyOrig_int(:,1) >= 1 & xyOrig_int(:,1) <= nc/scale-1 & xyOrig_int(:,2) >= 1 & xyOrig_int(:,2) <= nr/scale-1;
                    xyOrig_int_in = xyOrig_int(flag_in,:);
                    xyOrig_frac_in = xyOrig_frac(flag_in,:);
                    if useNearest
                        foundIndIn = foundInd(idIn);
                        origInd = sub2ind([nr,nc],xyOrig(idIn,2),xyOrig(idIn,1));
                        
                        imgRect([foundIndIn;numPix + foundIndIn;2*numPix + foundIndIn]) = imgOrig([origInd;numPix + origInd;2*numPix + origInd]);
                    else
                        foundIndIn = foundInd(flag_in);
                        [yyy,xxx] = ind2sub([nr,nc],foundIndIn);
                        foundIndIn = sub2ind([nr,nc]/scale, yyy/scale, xxx/scale);
                        imgRect = BilinearRemap(imgRect, imgR, imgG, imgB, image, xyOrig_int_in, xyOrig_frac_in,foundIndIn);
                        
                    end
                    
                    
                    figure(figNum);clf;imshow(imgRect);
                catch
                    advfl = 1;
                end
                drawnow;
            end
        end
        %         if ~isempty(pixHoldup1) || ~isempty(pixHolddown1)
        if ~isempty(pixLeft)
            if draw
                
                figure(333),plot(pixLeft(:,1),pixLeft(:,2),'.','Color',rand(1,3));
                %                 foundInd = sub2ind([nr,nc],pixLeft(:,2),pixLeft(:,1));
                foundInd = sub2ind([nr/1,nc/1],pixLeft(:,2)./1,pixLeft(:,1)./1);
                
                try
                    xyOrig = round([xRect2Orig(foundInd) yRect2Orig(foundInd)]);
                    idIn = xyOrig(:,1) > 0 & xyOrig(:,1) <= nc & xyOrig(:,2) > 0 & xyOrig(:,2) <= nr;
                    
                    
                    xyOrig_ = [xRect2Orig(foundInd) yRect2Orig(foundInd)];
                    [xyOrig_, ~] = Num2Fix_unsigned([xRect2Orig(foundInd) yRect2Orig(foundInd)]./scale, 20,9);
                    
                    xyOrig_int = floor(xyOrig_);
                    xyOrig_frac = xyOrig_ - xyOrig_int;
                    
                    y_in = [ySampled(i);ySampled(i+1)];
                    
                    %                     flag_in = xyOrig_int(:,1) >= scale & xyOrig_int(:,1) <= nc/scale-1 & xyOrig_int(:,2) >= scale & xyOrig_int(:,2) <= nr/scale-1;
                    flag_in = xyOrig_int(:,1) >= 1 & xyOrig_int(:,1) <= nc/scale-1 & xyOrig_int(:,2) >= 1 & xyOrig_int(:,2) <= nr/scale-1;
                    xyOrig_int_in = xyOrig_int(flag_in,:);
                    xyOrig_frac_in = xyOrig_frac(flag_in,:);
                    
                    
                    if useNearest
                        foundIndIn = foundInd(idIn);
                        origInd = sub2ind([nr,nc],xyOrig(idIn,2),xyOrig(idIn,1));
                        
                        imgRect([foundIndIn;numPix + foundIndIn;2*numPix + foundIndIn]) = imgOrig([origInd;numPix + origInd;2*numPix + origInd]);
                    else
                        foundIndIn = foundInd(flag_in);
                        [yyy,xxx] = ind2sub([nr,nc],foundIndIn);
                        foundIndIn = sub2ind([nr,nc]/scale, yyy/scale, xxx/scale);
                        imgRect = BilinearRemap(imgRect, imgR, imgG, imgB, image, xyOrig_int_in, xyOrig_frac_in,foundIndIn);
                        
                    end
                    figure(figNum);clf;imshow(imgRect);
                catch
                    dghkl = 1;
                end
                drawnow;
            end
        end
        
    end
    if mod(i,200) == 0
        i
        %         save('homo2.mat','StoreData','pixWrite');
        if scale == 1
            save(fullfile(paraDir,'homo2.mat'),'StoreData','pixWrite','ZeroFlag','Pix1','Pix2','Pix3','Pix4','Pix5');
        else
            save(fullfile(paraDir,'homo3.mat'),'StoreData','pixWrite','ZeroFlag','Pix1','Pix2','Pix3','Pix4','Pix5');
        end
        
        bvhjjkh = 1;
    end
    %     continue;
end

pixAll_1 = [xuuAllUse(:) yuuAllUse(:)];
flagAll_1 = find(pixAll_1(:,1) >=crop1(1) & pixAll_1(:,1)<=crop2(1) & pixAll_1(:,2) >=crop1(2) & pixAll_1(:,2) <= crop2(2));




return;


if 1
    ovlp = size(pixWrite,1)-size(intersect(pixWrite,pixWrite,'rows'),1);
    img = zeros(nr,nc);
    ind_ = sub2ind([nr nc],pixWrite(:,2),pixWrite(:,1));
    indind = find(pixWrite(:,1) >=crop1(1) & pixWrite(:,1)<=crop2(1) & pixWrite(:,2) >=crop1(2) & pixWrite(:,2) <= crop2(2));
    ind = ind_(indind);
    img(ind) = 1;
    figure,imshow(img);
    foundInd = ind;
    xyOrig_ = [xRect2Orig(foundInd) yRect2Orig(foundInd)];
    [xyOrig_, ~] = Num2Fix_unsigned([xRect2Orig(foundInd) yRect2Orig(foundInd)]./scale, 20,9);
    xyOrig_int = floor(xyOrig_);
    xyOrig_frac = xyOrig_ - xyOrig_int;
    flag_in = xyOrig_int(:,1) >= 1 & xyOrig_int(:,1) <= nc/scale-1 & xyOrig_int(:,2) >= 1 & xyOrig_int(:,2) <= nr/scale-1;
    sum(flag_in)
    figure,plot(xyOrig_int(:,1))
    figure,plot(xyOrig_int(:,2))
    sum(flag_in)
    xyOrig_int_in = xyOrig_int;
    xyOrig_frac_in = xyOrig_frac(flag_in,:);
    xyOrig_frac_in = xyOrig_frac;
    foundIndIn = foundInd;
    %     imgRect = BilinearRemap(imgRect, imgR, imgG, imgB, image, xyOrig_int_in, xyOrig_frac_in,foundIndIn);
    %     figure,plot(PixNew3_Bilinear_floor1)
    %     max(PixNew3_Bilinear_floor1)
    flag_in = xyOrig_int(:,1) >= 1 & xyOrig_int(:,1) <= nc/scale-1 & xyOrig_int(:,2) >= 1 & xyOrig_int(:,2) <= nr/scale-1;
    xyOrig_int_in = xyOrig_int(flag_in,:);
    xyOrig_frac_in = xyOrig_frac(flag_in,:);
    foundIndIn = foundInd(flag_in);
    [yyy,xxx] = ind2sub([nr,nc],foundIndIn);
    foundIndIn = sub2ind([nr,nc]/scale, yyy/scale, xxx/scale);
    imgRect = BilinearRemap(imgRect, imgR, imgG, imgB, image, xyOrig_int_in, xyOrig_frac_in,foundIndIn);
    figure,imshow(imgRect)
    imgL2 = imgRect;
    
    if 0
        if crop1(2)-1 ~=1
            imgL2(1:crop1(2)-1,:,1:3) = 0;
        end
        if crop2(2)+1 ~= nr
            imgL2(crop2(2)+1:end,:,1:3) = 0;
        end
        if crop2(1)+1 ~=nc
            imgL2(:,crop2(1)+1:end,1:3) = 0;
        end
    end
    % %     imgL2 = RectImg(xOrig2RectL, yOrig2RectL, xRect2OrigL, yRect2OrigL, imgLeft,KK_newL, KinitL,kL,rotMatLeft,imgLeft,paraDir,'Left');
end









ovlp = size(pixWrite,1)-size(intersect(pixWrite,pixWrite,'rows'),1);
img = zeros(nr,nc);
ind = sub2ind([nr nc],pixWrite(:,2),pixWrite(:,1));
img(ind) = 1;
figure,imshow(img);

indShow = intersect(flagAll_1, ind);

[yyy1,xxx1] = ind2sub([nr,nc],indShow);
foundIndIn1 = sub2ind([nr,nc]/scale, yyy1/scale, xxx1/scale);
foundIndOut1 = setdiff([1 : nc*nr/scale/scale]',foundIndIn1);

imgRect([foundIndOut1; foundIndOut1 + nc*nr/scale/scale; foundIndOut1 + 2*nc*nr/scale/scale]) = 0;

data = nan(nr,nc);
ind2 = sub2ind([nr nc],round(StoreData(:,2)),round(StoreData(:,1)));
data(ind2) = StoreData(:,3);
figure,imshow(data,[]);
figure,imshow(immultiply(data,data<10),[]);
figure,hist(StoreData(:,3),1000);

id1 = find(abs(diff(StoreData(:,2))) > 0); id1 = [id1; size(StoreData,1)];
id3 = find(abs(diff(StoreData(:,1))) ~= 1);id3 = [id3; size(StoreData,1)];
id = unique([id1;id3]);
id = unique([id;size(StoreData,1)]); did = diff(id); did = [id(1);did];
check = [did (ceil(did/16)*16)];
busEfficiency = size(StoreData,1)/sum((ceil(did/16)*16));
checkSum = sum(did) - size(pixWrite,1);
figure,hist(check(:,1),10000)
figure,hist(check(:,1)./check(:,2),10000)
busEfficiency2 = sum(check(:,1))/sum(check(:,2));
if scale == 1
    save(fullfile(paraDir,'homo2.mat'),'StoreData','pixWrite','ZeroFlag','Pix1','Pix2','Pix3','Pix4','Pix5');
else
    save(fullfile(paraDir,'homo3.mat'),'StoreData','pixWrite','ZeroFlag','Pix1','Pix2','Pix3','Pix4','Pix5');
end
if 0
    nr = 1080*2; nc = 1920*2;
    idd = find(StoreData(:,1) <= 2);
    figure,imshow(zeros(nr,nc));hold on;
    for i = idd(1) : 10000000;
        plot(StoreData(i,1),StoreData(i,2),'.g');
        drawnow;
    end
end


end


function pixNew = plotUndistortion(nc, nr, Kinit, kk_new,kc,R,skipNum)


% [mx,my] = meshgrid(1:nc/skipNum:(nc-0),1:nr/skipNum:(nr-0));
[mx,my] = meshgrid(0:nc/skipNum:(nc-0),0:nr/skipNum:(nr-0));
[nnx,nny]=size(mx);
px=reshape(mx',nnx*nny,1);
py=reshape(my',nnx*nny,1);
% % kk_new=[fc(1) alpha_c*fc(1) cc(1);0 fc(2) cc(2);0 0 1];
rays=inv(kk_new)*[px';py';ones(1,length(px))];
rays = R'*rays;


x=[rays(1,:)./rays(3,:);rays(2,:)./rays(3,:)];


title2=strcat('Complete Distortion Model');

% % fh1 = 2;

%if ishandle(fh1),
%    close(fh1);
%end;
% % % % % % figure; clf;
xd=apply_distortion(x,kc);
px2=Kinit(1,1)*xd(1,:)+Kinit(1,3);
py2=Kinit(2,2)*xd(2,:)+Kinit(2,3);

idIn = find(px2>=0 & px2<=nc & py2 >=0 & py2 <= nr);
idIn = [1:size(px2,2)]';

pixNew = [px2 py2];
dx=px2(idIn)'-px(idIn);
dy=py2(idIn)'-py(idIn);
% Q=quiver(px+1,py+1,dx,dy);
Q=quiver(px(idIn)+0,py(idIn)+0,dx,dy);
hold on;
% % plot(px(idIn),py(idIn),'.r');
plot(Kinit(1,3),Kinit(2,3),'o');
plot((nc-1)/2+1,(nr-1)/2+1,'x');
dr=reshape(sqrt((dx.*dx)+(dy.*dy)),nny,nnx)';
[C,h]=contour(mx,my,dr,'k');
clabel(C,h);
Mean=mean(mean(dr));
Max=max(max(dr));
title(title2);
axis equal;
end
function pixDist = remapRect(pixRect, KRect, KUndist,distCoeff, RL)

alpha = 0;

rays = inv(KRect)*pextend(pixRect);


% Rotation: (or affine transformation):

rays2 = RL'*rays;

x = [rays2(1,:)./rays2(3,:);rays2(2,:)./rays2(3,:)];


% Add distortion:
xd = apply_distortion(x,distCoeff);


% Reconvert in pixels:

px2_ = KUndist(1,1)*(xd(1,:) + alpha*xd(2,:)) + KUndist(1,3);
py2_ = KUndist(2,2)*xd(2,:) + KUndist(2,3);
pixDist = [px2_;py2_]';

end


function [v, indexx] = splineval(xPointTmp, params, xx)
coeff = params(:,4:-1:1);
% %     [Y, xx] = spline_(xPointTmp', params);
[~,index] = histc(xx',[-inf,xPointTmp(2:end-1),inf]);
xxNew = xx-xPointTmp(index)';

v = coeff(index,1);

indexx = unique(index);
for ii=2:4
    v = xxNew(:).*v + coeff(index,ii);
end
v = reshape(v,[1,length(xx)])';
end

function g1 = bicubic(f,kx,ky)
aa = [ky,0,0;0,kx,0;0.5*(1-ky),0.5*(1-kx),1];
ff = f;
[m,n]=size(f);
a=f(1,:);
c=f(m,:);             %将待插值图像矩阵前后各扩展两行两列,共扩展四行四列，当插值点为边沿点时候，确保周围有16个点
b=[f(1,1),f(1,1),f(:,1)',f(m,1),f(m,1)];
d=[f(1,n),f(1,n),f(:,n)',f(m,n),f(m,n)];
a1=[a;a;f;c;c];
b1=[b;b;a1';d;d];
f=b1';
f1=double(f);
AMat = [];
BMat = [];
for i=1:ky*m                 %利用双三次插值公式对新图象所有像素赋值
    if 1
        if 1
            temp = [i 1 1]*inv(aa);
            u = temp(1)-floor(temp(1));
            i1 = floor(temp(1)) + 2;
        else
            temp = [i 1]*inv(aa(1:2,1:2));
            u = temp(1)-floor(temp(1));
            i1 = floor(temp(1)) + 2;
        end
    else
        u=rem(i,ky)/ky;
        i1=floor(i/ky)+2;%rem取余，取余后再除以k是因为u为小数，自己代数字验证一下
    end
    %i1为B(X,Y)对应A(x,y)的横坐标整数部分，加2是因为上面扩大了四行四列，要回到原来图像的点再计算。
    A=[sw(1+u) sw(u) sw(1-u) sw(2-u)];   %四个横坐标的权重W(i)
    for j=1:kx*n
        
        
        if i == 200 && j == 200
            savdkbh = 1;
        end
        
        if 1
            if 1
                temp = [i j 1]*inv(aa);
                v = temp(2)-floor(temp(2));
                j1 = floor(temp(2)) + 2;
            else
                temp = [i j]*inv(aa(1:2,1:2));
                %             temp = [(i+0.5)*(1/kx)-0.5  (j+0.5)*(1/ky)-0.5];
                v = temp(2)-floor(temp(2));
                j1 = floor(temp(2)) + 2;
            end
        else
            v=rem(j,kx)/kx;   %纵坐标原理同上
            j1=floor(j/kx)+2;
        end
        
        C=[sw(1+v);sw(v);sw(1-v);sw(2-v)]; %转置
        B=[f1(i1-1,j1-1)     f1(i1-1,j1)     f1(i1-1,j1+1)     f1(i1-1,j1+2)    %坐标P(x+u,y+v)最近的16个点的像素值
            f1(i1    ,j1-1)     f1(i1,   j1)     f1(i1,   j1+1)     f1(i1,   j1+2)
            f1(i1+1,j1-1)     f1(i1+1,j1)    f1(i1+1,j1+1)   f1(i1+1,j1+2)
            f1(i1+2,j1-1)     f1(i1+2,j1)    f1(i1+2,j1+1)   f1(i1+2,j1+2)];
        g1(i,j)=(A*B*C);
        AMat = [AMat;A];
        BMat = [BMat;C'];
    end
end
cc = imresize(ff,[size(ff,1)*ky size(ff,2)*kx],'Method','bicubic');
figure,imshow(abs(cc(2:end-1,2:end-1)-g1(2:end-1,2:end-1)),[]);

end


function B = bicubic8(A,kx,ky)

a = [ky,0,0;0,kx,0;0.5*(1-ky),0.5*(1-kx),1];
[m,n]=size(A);
for i=1:ky*m
    for j=1:kx*n
        temp = [i j 1]*inv(a);
        di=temp(1)-floor(temp(1));
        dj=temp(2)-floor(temp(2));
        y = floor(temp(1));
        x = floor(temp(2));
        aX=max(x-3,1);
        bX=max(x-2,1);
        cX=max(x-1,1);
        dX=max(x,1);
        eX=min(x+1,n);
        fX=min(x+2,n);
        gX=min(x+3,n);
        hX=min(x+4,n);
        aY=max(y-3,1);
        bY=max(y-2,1);
        cY=max(y-1,1);
        dY=max(y,1);
        eY=min(y+1,m);
        fY=min(y+2,m);
        gY=min(y+3,m);
        hY=min(y+4,m);
        
        %  first calling for each of the 8 rows (X direction
        %  interpolation)
        temp1=bicubic8x8([A(aY,aX),A(aY,bX),A(aY,cX),A(aY,dX),A(aY,eX),A(aY,fX),A(aY,gX),A(aY,hX)],dj);
        temp2=bicubic8x8([A(bY,aX),A(bY,bX),A(bY,cX),A(bY,dX),A(bY,eX),A(bY,fX),A(bY,gX),A(bY,hX)],dj);
        temp3=bicubic8x8([A(cY,aX),A(cY,bX),A(cY,cX),A(cY,dX),A(cY,eX),A(cY,fX),A(cY,gX),A(cY,hX)],dj);
        temp4=bicubic8x8([A(dY,aX),A(dY,bX),A(dY,cX),A(dY,dX),A(dY,eX),A(dY,fX),A(dY,gX),A(dY,hX)],dj);
        temp5=bicubic8x8([A(eY,aX),A(eY,bX),A(eY,cX),A(eY,dX),A(eY,eX),A(eY,fX),A(eY,gX),A(eY,hX)],dj);
        temp6=bicubic8x8([A(fY,aX),A(fY,bX),A(fY,cX),A(fY,dX),A(fY,eX),A(fY,fX),A(fY,gX),A(fY,hX)],dj);
        temp7=bicubic8x8([A(gY,aX),A(gY,bX),A(gY,cX),A(gY,dX),A(gY,eX),A(gY,fX),A(gY,gX),A(gY,hX)],dj);
        temp8=bicubic8x8([A(hY,aX),A(hY,bX),A(hY,cX),A(hY,dX),A(hY,eX),A(hY,fX),A(hY,gX),A(hY,hX)],dj);
        
        %  Then calling bicubic for these 8 interpolated values (Y
        %  direction interpolation)
        B(i,j)=bicubic8x8([temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8],di);
    end
end
end
function [X] = bicubic8x8(v,fact)
X = 0;
if fact==0
    X = v(4);
    return
end
Y = [3 2 1 0 1 2 3 4];
Z = [1 1 1 1 -1 -1 -1 -1];
% 8x8 bicubic polynomial
% r(x) =  (67/56)|x|^3 - (123/56)|x|^2              +  1     , 0 <= |x| < 1
% r(x) = -(33/56)|x|^3 + (177/56)|x|^2 - (75/14)|x| + (39/14), 1 <= |x| < 2
% r(x) =  (9/56)|x|^3  - (75/56)|x|^2  + (51/14)|x| - (45/14), 2 <= |x| < 3
% r(x) = -(3/56)|x|^3  + (33/56)|x|^2  - (15/7)|x|  + (18/7) , 3 <= |x| < 4
% The sum of r(3+fact)*v1 + r(2+fact)*v2 + r(1+fact)*v3 + r(fact)*v4 +
% r(1-fact)*v5 + r(2-fact)*v6 + r(3-fact)*v7 + r(4-fact)*v8

A = 67/56;
B = -123/56;
C = 1;
D = -33/56;
E = 177/56;
F = -75/14;
G = 39/14;
H = 9/56;
I = -75/56;
J = 51/14;
K = -45/14;
L = -3/56;
M = 33/56;
N = -15/7;
O = 18/7;

for i = 1:8
    if v(i)~=0
        Fi = Y(i)+fact*Z(i);
        if Fi<1
            X = X + v(i)*(1 + Fi*Fi*(B + Fi*A));
        elseif Fi<2
            X = X + v(i)*(G + Fi*(F + Fi*(E + Fi*D)));
        elseif Fi<3
            X = X + v(i)*(K + Fi*(J + Fi*(I + Fi*H)));
        elseif Fi<4
            X = X + v(i)*(O + Fi*(N + Fi*(M + Fi*L)));
        end
    end
end
end



function A=sw(w1)
w=abs(w1);
a=-0.5;
if w<1&&w>=0
    A=1-(a+3)*w^2+(a+2)*w^3;
else if w>=1&&w<2
        A=a*w^3-5*a*w^2+(8*a)*w-4*a;
    else
        A=0;
    end
end
end


function B = bicubic4x4FastUse(A, kx,ky)
afact = -0.5;
[m,n] = size(A);
NewRows = size(A,1)*ky;
NewColumns = size(A,2)*kx;
a = [ky,0,0;0,kx,0;0.5*(1-ky),0.5*(1-kx),1];
[xuu,yuu] = meshgrid([1:NewRows],[1:NewColumns]);
%         Temp = [repmat(i,NewColumns,1) [1:NewColumns]' ones(NewColumns, 1)]*inv(a);
Temp = [xuu(:) yuu(:) ones(NewColumns*NewRows, 1)]*inv(a);
Di=Temp(:,1)-floor(Temp(:,1)); % di = yfract
Dj=Temp(:,2)-floor(Temp(:,2)); % dj = xfract
Y = floor(Temp(:,1));
X = floor(Temp(:,2));
AX=max(X-1,1);
BX=max(X,1);
CX=min(X+1,n);
DX=min(X+2,n);
AY=max(Y-1,1);
BY=max(Y,1);
CY=min(Y+1,m);
DY=min(Y+2,m);
ind1 = sub2ind(size(A),AY,AX);
ind2 = sub2ind(size(A),AY,BX);
ind3 = sub2ind(size(A),AY,CX);
ind4 = sub2ind(size(A),AY,DX);
ind5 = sub2ind(size(A),BY,AX);
ind6 = sub2ind(size(A),BY,BX);
ind7 = sub2ind(size(A),BY,CX);
ind8 = sub2ind(size(A),BY,DX);
ind9 = sub2ind(size(A),CY,AX);
ind10 = sub2ind(size(A),CY,BX);
ind11 = sub2ind(size(A),CY,CX);
ind12 = sub2ind(size(A),CY,DX);
ind13 = sub2ind(size(A),DY,AX);
ind14 = sub2ind(size(A),DY,BX);
ind15 = sub2ind(size(A),DY,CX);
ind16 = sub2ind(size(A),DY,DX);
Temp1=bicubic4x4Fast([A(ind1),A(ind2),A(ind3),A(ind4)],Dj,afact);
Temp2=bicubic4x4Fast([A(ind5),A(ind6),A(ind7),A(ind8)],Dj,afact);
Temp3=bicubic4x4Fast([A(ind9),A(ind10),A(ind11),A(ind12)],Dj,afact);
Temp4=bicubic4x4Fast([A(ind13),A(ind14),A(ind15),A(ind16)],Dj,afact);
%         B(i,:)=bicubic4x4Fast([Temp1,Temp2,Temp3,Temp4],Di,afact);
B=reshape(bicubic4x4Fast([Temp1,Temp2,Temp3,Temp4],Di,afact),NewColumns,NewRows)';
end
function [X] = bicubic4x4Fast(v,fact,afact)

X = zeros(size(v,1),1);
if fact==0
    X = v(2);
    return
end
Y = [1 0 1 2];
Z = [1 1 -1 -1];
% r(x) = (a + 2)|x|^3 - (a + 3)|x|^2         +  1 , 0 <= |x| < 1
% r(x) =       a|x|^3 -      5a|x|^2 + 8a|x| - 4a , 1 <= |x| < 2
% where x is {fact,fact+1,abs(fact-1),abs(fact-2)}
% then r(x) will contain the interpolation ratio which will be multiplied
% by the point value.
% The sum of r(fact+1)*v1+r(fact)*v2+r(abs(fact-1))*v3+r(abs(fact-2))*v4
% will give the interpolated value.
% The value a is taken as -0.5 as per Matlab. Play
% around with the value of a and you get some interesting results.
% ALL THE CUBIC EQUATION HAVE BEEN PICKED FROM RESEARCH PAPERS I DO NOT
% KNOW THEIR DERIVATION OR WHY THE FACTOR a EXISTS IN 4x4 but not in the
% others. The only reason I can think of is that for 6x6 & 8x8 bicubic
% equations are continous till the second differential that
% is : F(x-) = F(x+), F'(x-) = F'(x+) & is F''(x-) = F''(x+)
% while 4x4 bicubic is only : F(x-) = F(x+), F'(x-) = F'(x+)

% a = -0.5;
a = afact;
A=(2+a);
B=-(3+a);
C=1;
D=a;
E=-5*a;
F=8*a;
G=-4*a;

% for i = 1:4
%     if v(i)~=0
%         Fi = Y(i)+fact*Z(i);
%         if Fi<1
%             X = X + v(i)*(1 + Fi*Fi*(B + Fi*A));
%         elseif Fi<2
%             X = X + v(i)*(G + Fi*(F + Fi*(E + Fi*D)));
%         end
%     end
% end


% flag = v~=0;
for i = 1:4
    if 1 %v(i)~=0
        Fi = Y(i)+fact*Z(i);
        if sum(Fi<1)~=0
            X(Fi<1,1) = X(Fi<1,1) + v(Fi<1,i).*(1 + Fi(Fi<1).^2.*(B + Fi(Fi<1.)*A));
        end
        
        if sum(Fi>=1 & Fi<2)~=0
            X(Fi>=1 & Fi<2,1) = X(Fi>=1 & Fi<2,1) + v(Fi>=1 & Fi<2,i).*(G + Fi(Fi>=1 & Fi<2).*(F + Fi(Fi>=1 & Fi<2).*(E + Fi(Fi>=1 & Fi<2).*D)));
        end
        % %         if Fi<1
        % %             X = X + v(i)*(1 + Fi*Fi*(B + Fi*A));
        % %         elseif Fi<2
        % %             X = X + v(i)*(G + Fi*(F + Fi*(E + Fi*D)));
        % %         end
    end
end
end
function B = bicubic8x8FastUse(A, kx,ky)
% afact = -0.5;
[m,n] = size(A);
NewRows = size(A,1)*ky;
NewColumns = size(A,2)*kx;
a = [ky,0,0;0,kx,0;0.5*(1-ky),0.5*(1-kx),1];
[xuu,yuu] = meshgrid([1:NewRows],[1:NewColumns]);
%         Temp = [repmat(i,NewColumns,1) [1:NewColumns]' ones(NewColumns, 1)]*inv(a);
Temp = [xuu(:) yuu(:) ones(NewColumns*NewRows, 1)]*inv(a);
Di=Temp(:,1)-floor(Temp(:,1)); % di = yfract
Dj=Temp(:,2)-floor(Temp(:,2)); % dj = xfract
Y = floor(Temp(:,1));
X = floor(Temp(:,2));

AX=max(X-3,1);
BX=max(X-2,1);
CX=max(X-1,1);
DX=max(X,1);
EX=min(X+1,n);
FX=min(X+2,n);
GX=min(X+3,n);
HX=min(X+4,n);
AY=max(Y-3,1);
BY=max(Y-2,1);
CY=max(Y-1,1);
DY=max(Y,1);
EY=min(Y+1,m);
FY=min(Y+2,m);
GY=min(Y+3,m);
HY=min(Y+4,m);


ind1 = sub2ind(size(A),AY,AX);
ind2 = sub2ind(size(A),AY,BX);
ind3 = sub2ind(size(A),AY,CX);
ind4 = sub2ind(size(A),AY,DX);
ind5 = sub2ind(size(A),AY,EX);
ind6 = sub2ind(size(A),AY,FX);
ind7 = sub2ind(size(A),AY,GX);
ind8 = sub2ind(size(A),AY,HX);

ind9 = sub2ind(size(A),BY,AX);
ind10 = sub2ind(size(A),BY,BX);
ind11 = sub2ind(size(A),BY,CX);
ind12 = sub2ind(size(A),BY,DX);
ind13 = sub2ind(size(A),BY,EX);
ind14 = sub2ind(size(A),BY,FX);
ind15= sub2ind(size(A),BY,GX);
ind16 = sub2ind(size(A),BY,HX);

ind17 = sub2ind(size(A),CY,AX);
ind18 = sub2ind(size(A),CY,BX);
ind19 = sub2ind(size(A),CY,CX);
ind20 = sub2ind(size(A),CY,DX);
ind21 = sub2ind(size(A),CY,EX);
ind22 = sub2ind(size(A),CY,FX);
ind23 = sub2ind(size(A),CY,GX);
ind24 = sub2ind(size(A),CY,HX);

ind25 = sub2ind(size(A),DY,AX);
ind26 = sub2ind(size(A),DY,BX);
ind27 = sub2ind(size(A),DY,CX);
ind28 = sub2ind(size(A),DY,DX);
ind29 = sub2ind(size(A),DY,EX);
ind30 = sub2ind(size(A),DY,FX);
ind31 = sub2ind(size(A),DY,GX);
ind32 = sub2ind(size(A),DY,HX);

ind33 = sub2ind(size(A),EY,AX);
ind34 = sub2ind(size(A),EY,BX);
ind35 = sub2ind(size(A),EY,CX);
ind36 = sub2ind(size(A),EY,DX);
ind37 = sub2ind(size(A),EY,EX);
ind38 = sub2ind(size(A),EY,FX);
ind39 = sub2ind(size(A),EY,GX);
ind40 = sub2ind(size(A),EY,HX);

ind41 = sub2ind(size(A),FY,AX);
ind42 = sub2ind(size(A),FY,BX);
ind43 = sub2ind(size(A),FY,CX);
ind44 = sub2ind(size(A),FY,DX);
ind45 = sub2ind(size(A),FY,EX);
ind46 = sub2ind(size(A),FY,FX);
ind47 = sub2ind(size(A),FY,GX);
ind48 = sub2ind(size(A),FY,HX);

ind49 = sub2ind(size(A),GY,AX);
ind50 = sub2ind(size(A),GY,BX);
ind51 = sub2ind(size(A),GY,CX);
ind52 = sub2ind(size(A),GY,DX);
ind53 = sub2ind(size(A),GY,EX);
ind54 = sub2ind(size(A),GY,FX);
ind55 = sub2ind(size(A),GY,GX);
ind56 = sub2ind(size(A),GY,HX);

ind57 = sub2ind(size(A),HY,AX);
ind58 = sub2ind(size(A),HY,BX);
ind59 = sub2ind(size(A),HY,CX);
ind60 = sub2ind(size(A),HY,DX);
ind61 = sub2ind(size(A),HY,EX);
ind62 = sub2ind(size(A),HY,FX);
ind63 = sub2ind(size(A),HY,GX);
ind64 = sub2ind(size(A),HY,HX);


Temp1=bicubic8x8Fast([A(ind1),A(ind2),A(ind3),A(ind4),A(ind5),A(ind6),A(ind7),A(ind8)],Dj);
Temp2=bicubic8x8Fast([A(ind9),A(ind10),A(ind11),A(ind12),A(ind13),A(ind14),A(ind15),A(ind16)],Dj);
Temp3=bicubic8x8Fast([A(ind17),A(ind18),A(ind19),A(ind20),A(ind21),A(ind22),A(ind23),A(ind24)],Dj);
Temp4=bicubic8x8Fast([A(ind25),A(ind26),A(ind27),A(ind28),A(ind29),A(ind30),A(ind31),A(ind32)],Dj);
Temp5=bicubic8x8Fast([A(ind33),A(ind34),A(ind35),A(ind36),A(ind37),A(ind38),A(ind39),A(ind40)],Dj);
Temp6=bicubic8x8Fast([A(ind41),A(ind42),A(ind43),A(ind44),A(ind45),A(ind46),A(ind47),A(ind48)],Dj);
Temp7=bicubic8x8Fast([A(ind49),A(ind50),A(ind51),A(ind52),A(ind53),A(ind54),A(ind55),A(ind56)],Dj);
Temp8=bicubic8x8Fast([A(ind57),A(ind58),A(ind59),A(ind60),A(ind61),A(ind62),A(ind63),A(ind64)],Dj);

%  Then calling bicubic for these 8 interpolated values (Y
%  direction interpolation)

%         B(i,:)=bicubic8x8Fast([Temp1,Temp2,Temp3,Temp4],Di,afact);
B=reshape(bicubic8x8Fast([Temp1,Temp2,Temp3,Temp4,Temp5,Temp6,Temp7,Temp8],Di),NewColumns,NewRows)';
end
function [X] = bicubic8x8Fast(v,fact)
X = zeros(size(v,1),1);
if fact==0
    X = v(4);
    return
end
Y = [3 2 1 0 1 2 3 4];
Z = [1 1 1 1 -1 -1 -1 -1];
% 8x8 bicubic polynomial
% r(x) =  (67/56)|x|^3 - (123/56)|x|^2              +  1     , 0 <= |x| < 1
% r(x) = -(33/56)|x|^3 + (177/56)|x|^2 - (75/14)|x| + (39/14), 1 <= |x| < 2
% r(x) =  (9/56)|x|^3  - (75/56)|x|^2  + (51/14)|x| - (45/14), 2 <= |x| < 3
% r(x) = -(3/56)|x|^3  + (33/56)|x|^2  - (15/7)|x|  + (18/7) , 3 <= |x| < 4
% The sum of r(3+fact)*v1 + r(2+fact)*v2 + r(1+fact)*v3 + r(fact)*v4 +
% r(1-fact)*v5 + r(2-fact)*v6 + r(3-fact)*v7 + r(4-fact)*v8

A = 67/56;
B = -123/56;
C = 1;
D = -33/56;
E = 177/56;
F = -75/14;
G = 39/14;
H = 9/56;
I = -75/56;
J = 51/14;
K = -45/14;
L = -3/56;
M = 33/56;
N = -15/7;
O = 18/7;

for i = 1:8
    %     if v(i)~=0
    %         Fi = Y(i)+fact*Z(i);
    %         if Fi<1
    %             X = X + v(i)*(1 + Fi*Fi*(B + Fi*A));
    %         elseif Fi<2
    %             X = X + v(i)*(G + Fi*(F + Fi*(E + Fi*D)));
    %         elseif Fi<3
    %             X = X + v(i)*(K + Fi*(J + Fi*(I + Fi*H)));
    %         elseif Fi<4
    %             X = X + v(i)*(O + Fi*(N + Fi*(M + Fi*L)));
    %         end
    %     end
    if 1 % v(i)~=0
        Fi = Y(i)+fact*Z(i);
        
        if sum(Fi<1)~=0
            X(Fi<1,1) = X(Fi<1,1) + v(Fi<1,i).*(1 + Fi(Fi<1).^2.*(B + Fi(Fi<1.)*A));
        end
        flag2 = Fi>=1 & Fi<2;
        if sum(flag2) ~= 0
            X(flag2) = X(flag2) + v(flag2,i).*(G + Fi(flag2).*(F + Fi(flag2).*(E + Fi(flag2).*D)));
        end
        flag3 = Fi>=2 & Fi<3;
        if sum(flag3) ~= 0
            X(flag3) = X(flag3) + v(flag3,i).*(K + Fi(flag3).*(J + Fi(flag3).*(I + Fi(flag3).*H)));
        end
        flag4 = Fi>=3 & Fi<4;
        if sum(flag4) ~= 0
            X(flag4) = X(flag4) + v(flag4,i).*(O + Fi(flag4).*(N + Fi(flag4).*(M + Fi(flag4).*L)));
        end
        %
        %         if Fi<1
        %             X = X + v(i)*(1 + Fi*Fi*(B + Fi*A));
        %         elseif Fi<2
        %             X = X + v(i)*(G + Fi*(F + Fi*(E + Fi*D)));
        %         elseif Fi<3
        %             X = X + v(i)*(K + Fi*(J + Fi*(I + Fi*H)));
        %         elseif Fi<4
        %             X = X + v(i)*(O + Fi*(N + Fi*(M + Fi*L)));
        %         end
    end
end
end

function B = bicubic6x6FastUse(A,kx,ky)
[m,n] = size(A);
NewRows = size(A,1)*ky;
NewColumns = size(A,2)*kx;
a = [ky,0,0;0,kx,0;0.5*(1-ky),0.5*(1-kx),1];
[xuu,yuu] = meshgrid([1:NewRows],[1:NewColumns]);
%         Temp = [repmat(i,NewColumns,1) [1:NewColumns]' ones(NewColumns, 1)]*inv(a);
Temp = [xuu(:) yuu(:) ones(NewColumns*NewRows, 1)]*inv(a);
Di=Temp(:,1)-floor(Temp(:,1)); % di = yfract
Dj=Temp(:,2)-floor(Temp(:,2)); % dj = xfract
Y = floor(Temp(:,1));
X = floor(Temp(:,2));

AX=max(X-2,1);
BX=max(X-1,1);
CX=max(X,1);
DX=min(X+1,n);
EX=min(X+2,n);
FX=min(X+3,n);

AY=max(Y-2,1);
BY=max(Y-1,1);
CY=max(Y,1);
DY=min(Y+1,m);
EY=min(Y+2,m);
FY=min(Y+3,m);



ind1 = sub2ind(size(A),AY,AX);
ind2 = sub2ind(size(A),AY,BX);
ind3 = sub2ind(size(A),AY,CX);
ind4 = sub2ind(size(A),AY,DX);
ind5 = sub2ind(size(A),AY,EX);
ind6 = sub2ind(size(A),AY,FX);


ind7 = sub2ind(size(A),BY,AX);
ind8 = sub2ind(size(A),BY,BX);
ind9 = sub2ind(size(A),BY,CX);
ind10 = sub2ind(size(A),BY,DX);
ind11 = sub2ind(size(A),BY,EX);
ind12 = sub2ind(size(A),BY,FX);


ind13 = sub2ind(size(A),CY,AX);
ind14 = sub2ind(size(A),CY,BX);
ind15 = sub2ind(size(A),CY,CX);
ind16 = sub2ind(size(A),CY,DX);
ind17 = sub2ind(size(A),CY,EX);
ind18 = sub2ind(size(A),CY,FX);


ind19 = sub2ind(size(A),DY,AX);
ind20 = sub2ind(size(A),DY,BX);
ind21 = sub2ind(size(A),DY,CX);
ind22 = sub2ind(size(A),DY,DX);
ind23 = sub2ind(size(A),DY,EX);
ind24 = sub2ind(size(A),DY,FX);


ind25 = sub2ind(size(A),EY,AX);
ind26 = sub2ind(size(A),EY,BX);
ind27 = sub2ind(size(A),EY,CX);
ind28 = sub2ind(size(A),EY,DX);
ind29 = sub2ind(size(A),EY,EX);
ind30 = sub2ind(size(A),EY,FX);


ind31 = sub2ind(size(A),FY,AX);
ind32 = sub2ind(size(A),FY,BX);
ind33 = sub2ind(size(A),FY,CX);
ind34 = sub2ind(size(A),FY,DX);
ind35 = sub2ind(size(A),FY,EX);
ind36 = sub2ind(size(A),FY,FX);


Temp1=bicubic6x6Fast([A(ind1),A(ind2),A(ind3),A(ind4),A(ind5),A(ind6)],Dj);
Temp2=bicubic6x6Fast([A(ind7),A(ind8),A(ind9),A(ind10),A(ind11),A(ind12)],Dj);
Temp3=bicubic6x6Fast([A(ind13),A(ind14),A(ind15),A(ind16),A(ind17),A(ind18)],Dj);
Temp4=bicubic6x6Fast([A(ind19),A(ind20),A(ind21),A(ind22),A(ind23),A(ind24)],Dj);
Temp5=bicubic6x6Fast([A(ind25),A(ind26),A(ind27),A(ind28),A(ind29),A(ind30)],Dj);
Temp6=bicubic6x6Fast([A(ind31),A(ind32),A(ind33),A(ind34),A(ind35),A(ind36)],Dj);

B=reshape(bicubic6x6Fast([Temp1,Temp2,Temp3,Temp4,Temp5,Temp6],Di),NewColumns,NewRows)';
end
function [X] = bicubic6x6Fast(v,fact)
X = zeros(size(v,1),1);
if fact==0
    X = v(3);
    return
end
Y = [2 1 0 1 2 3];
Z = [1 1 1 -1 -1 -1];
% 6x6 bicubic polynomial
% r(x) =  (6/5)|x|^3 - (11/5)|x|^2              +  1    , 0 <= |x| < 1
% r(x) = -(3/5)|x|^3 + (16/5)|x|^2 -(27/5)|x|   +  14/5 , 1 <= |x| < 2
% r(x) =  (1/5)|x|^3 - (8/5)|x|^2  +(21/5)|x|   -  18/5 , 2 <= |x| < 3
%
% The sum of r(2+fact)*v1 + r(1+fact)*v2 + r(fact)*v3 + r(1-fact)*v4 +
% r(2-fact)*v5 + r(3-fact)*v6 represents the interpolated value

A =  4/3; 6/5;
B = -7/3; -11/5;
C =  1;
D = -7/12; -3/5;
E =  3; 16/5;
F = -59/12; -27/5;
G =  15/6; 14/5;
H =  1/12; 1/5;
I = -2/3; -8/5;
J =  21/12; 21/5;
K = -3/2; -18/5;

for i = 1:6
    %     if v(i)~=0
    %         Fi = Y(i)+fact*Z(i);
    %         if Fi<1
    %             X = X + v(i)*(1 + Fi*Fi*(B + Fi*A));
    %         elseif Fi<2
    %             X = X + v(i)*(G + Fi*(F + Fi*(E + Fi*D)));
    %         elseif Fi<3
    %             X = X + v(i)*(K + Fi*(J + Fi*(I + Fi*H)));
    %         end
    %     end
    
    if 1 %v(i)~=0
        Fi = Y(i)+fact*Z(i);
        flag1 = Fi < 1;
        if sum(flag1) ~= 0
            X(flag1) = X(flag1) + v(flag1,i).*(1 + Fi(flag1).*Fi(flag1).*(B + Fi(flag1).*A));
        end
        flag2 = Fi >= 1& Fi < 2;
        if sum(flag2) ~= 0
            X(flag2) = X(flag2) + v(flag2,i).*(G + Fi(flag2).*(F + Fi(flag2).*(E + Fi(flag2).*D)));
        end
        flag3 = Fi >= 2& Fi < 3;
        if sum(flag3) ~= 0
            X(flag3) = X(flag3) + v(flag3,i).*(K + Fi(flag3).*(J + Fi(flag3).*(I + Fi(flag3).*H)));
        end
        
        %
        %         if Fi<1
        %             X = X + v(i)*(1 + Fi*Fi*(B + Fi*A));
        %         elseif Fi<2
        %             X = X + v(i)*(G + Fi*(F + Fi*(E + Fi*D)));
        %         elseif Fi<3
        %             X = X + v(i)*(K + Fi*(J + Fi*(I + Fi*H)));
        %         end
    end
end


end

function pixRect = Orig2Rect(pix, KOrig, KRect, R,kc)

[pixUndist] = normalize_pixel(pix',[KOrig(1,1);KOrig(2,2)],[KOrig(1,3);KOrig(2,3)],kc,0);
pixUndistR = R*pextend(pixUndist);
pixRect = pflat(KRect*pixUndistR);
pixRect = pixRect(1:2,:)';



end

function pixDist = remapRectFishEye(pixRect, KRect, KUndist,distCoeff, RL,oCamModel)

alpha = 0;

rays = inv(KRect)*pextend(pixRect);


% Rotation: (or affine transformation):

rays2 = RL'*rays;

x = [rays2(1,:)./rays2(3,:);rays2(2,:)./rays2(3,:)];
if 0
    pixxDistort = pixDistort([2 1],:);
    % M2 = cam2world_fast(pixxDistort, calib_data.ocam_model);
    M2 = cam2world(pixxDistort, oCamModel);
    M2 = [M2(1,:)./M2(3,:);M2(2,:)./M2(3,:);ones(1,size(pixxDistort,2))];
    pixxUndist = intrMat_same*[-M2(2,:); -M2(1,:); M2(3,:)];
    pixUndist3D = [-M2(2,:); -M2(1,:); M2(3,:)];
    pixUndist = [pixxUndist(1,:)./pixxUndist(3,:); pixxUndist(2,:)./pixxUndist(3,:)];   %;ones(1,size(pixDistort,2))];
else
    % %     MM = inv(intrMat_same)*[pixUnDistort; ones(1,size(pixUnDistort,2))];
    MM = rays2;
    % % pixxDistort = pixDistort([2 1],:);
    % % % M2 = cam2world_fast(pixxDistort, calib_data.ocam_model);
    % % M2 = cam2world(pixxDistort, oCamModel);
    % % M2 = [M2(1,:)./M2(3,:);M2(2,:)./M2(3,:);ones(1,size(pixxDistort,2))];
    % % pixxUndist = intrMat_same*[-M2(2,:); -M2(1,:); M2(3,:)];
    % % pixUndist3D = [-M2(2,:); -M2(1,:); M2(3,:)];
    % % pixUndist = [pixxUndist(1,:)./pixxUndist(3,:); pixxUndist(2,:)./pixxUndist(3,:)];   %;ones(1,size(pixDistort,2))];
    % %
    % % MM = [-M2(2,:); -M2(1,:); M2(3,:)];
    MM_ = [-MM(2,:);-MM(1,:);MM(3,:)];
    [M4,~] = NormalizeVector(MM_');
    try
        % %         a
        M33 = world2cam_fast(-M4',oCamModel);
    catch
        M33 = world2cam(-M4',oCamModel);
    end
    M33_ = M33([2 1],:);
    pixDistorted = M33_';
    pixDist = pixDistorted;
end

% Add distortion:
% xd = apply_distortion(x,distCoeff);
%
%
% % Reconvert in pixels:
%
% px2_ = KUndist(1,1)*(xd(1,:) + alpha*xd(2,:)) + KUndist(1,3);
% py2_ = KUndist(2,2)*xd(2,:) + KUndist(2,3);
% pixDist = [px2_;py2_]';

end
function pixRect = Orig2RectFishEye(pix, KOrig, KRect, R,kc,oCamModel)

% % % [pixUndist] = normalize_pixel(pix',[KOrig(1,1);KOrig(2,2)],[KOrig(1,3);KOrig(2,3)],kc,0);

pixxDistort = pix(:,[2 1])';
% M2 = cam2world_fast(pixxDistort, calib_data.ocam_model);
M2 = cam2world(pixxDistort, oCamModel);
M2 = [M2(1,:)./M2(3,:);M2(2,:)./M2(3,:);ones(1,size(pixxDistort,2))];
% pixxUndist = KOrig*[-M2(2,:); -M2(1,:); M2(3,:)];
pixUndist3D = [-M2(2,:); -M2(1,:); M2(3,:)];
% pixUndist = [pixxUndist(1,:)./pixxUndist(3,:); pixxUndist(2,:)./pixxUndist(3,:)];

pixUndistR = R*(pixUndist3D);
pixRect = pflat(KRect*pixUndistR);
pixRect = pixRect(1:2,:)';



end

function pixNew = plotUndistortionFishEye(nc, nr, Kinit, kk_new,kc,R,skipNum,oCamModel)


% [mx,my] = meshgrid(1:nc/skipNum:(nc-0),1:nr/skipNum:(nr-0));
[mx,my] = meshgrid(0:nc/skipNum:(nc-0),0:nr/skipNum:(nr-0));
[nnx,nny]=size(mx);
px=reshape(mx',nnx*nny,1);
py=reshape(my',nnx*nny,1);
% % kk_new=[fc(1) alpha_c*fc(1) cc(1);0 fc(2) cc(2);0 0 1];
rays=inv(kk_new)*[px';py';ones(1,length(px))];
rays = R'*rays;


x=[rays(1,:)./rays(3,:);rays(2,:)./rays(3,:)];


title2=strcat('Complete Distortion Model');

% % fh1 = 2;

%if ishandle(fh1),
%    close(fh1);
%end;
% % % % % % figure; clf;
if 0
    xd=apply_distortion(x,kc);
    px2=Kinit(1,1)*xd(1,:)+Kinit(1,3);
    py2=Kinit(2,2)*xd(2,:)+Kinit(2,3);
else
    MM = rays;
    % % pixxDistort = pixDistort([2 1],:);
    % % % M2 = cam2world_fast(pixxDistort, calib_data.ocam_model);
    % % M2 = cam2world(pixxDistort, oCamModel);
    % % M2 = [M2(1,:)./M2(3,:);M2(2,:)./M2(3,:);ones(1,size(pixxDistort,2))];
    % % pixxUndist = intrMat_same*[-M2(2,:); -M2(1,:); M2(3,:)];
    % % pixUndist3D = [-M2(2,:); -M2(1,:); M2(3,:)];
    % % pixUndist = [pixxUndist(1,:)./pixxUndist(3,:); pixxUndist(2,:)./pixxUndist(3,:)];   %;ones(1,size(pixDistort,2))];
    % %
    % % MM = [-M2(2,:); -M2(1,:); M2(3,:)];
    MM_ = [-MM(2,:);-MM(1,:);MM(3,:)];
    [M4,~] = NormalizeVector(MM_');
    try
        M33 = world2cam_fast(-M4',oCamModel);
    catch
        M33 = world2cam(-M4',oCamModel);
    end
    M33_ = M33([2 1],:);
    pixDistorted = M33_';
    pixDist = pixDistorted;
    px2 = pixDist(:,1)';
    py2 = pixDist(:,2)';
end



idIn = find(px2>=0 & px2<=nc & py2 >=0 & py2 <= nr);
idIn = [1:size(px2,2)]';

pixNew = [px2 py2];
dx=px2(idIn)'-px(idIn);
dy=py2(idIn)'-py(idIn);
% Q=quiver(px+1,py+1,dx,dy);
Q=quiver(px(idIn)+0,py(idIn)+0,dx,dy);
hold on;
% % plot(px(idIn),py(idIn),'.r');
plot(Kinit(1,3),Kinit(2,3),'o');
plot((nc-1)/2+1,(nr-1)/2+1,'x');
dr=reshape(sqrt((dx.*dx)+(dy.*dy)),nny,nnx)';
[C,h]=contour(mx,my,dr,'k');
clabel(C,h);
Mean=mean(mean(dr));
Max=max(max(dr));
title(title2);
axis equal;
end

function imgRect = BilinearRemap(imgRect, imgR, imgG, imgB, image, xyOrig_int_in, xyOrig_frac_in,foundIndIn)

numPix = size(imgR,1)*size(imgR,2);

PixNew3_Bilinear_floor1 = xyOrig_int_in;
PixNew3_Bilinear_floor2 = [xyOrig_int_in(:,1)+1 xyOrig_int_in(:,2)];
PixNew3_Bilinear_floor3 = [xyOrig_int_in(:,1) xyOrig_int_in(:,2)+1];
PixNew3_Bilinear_floor4 = [xyOrig_int_in(:,1)+1 xyOrig_int_in(:,2)+1];

indBilinear1 = sub2ind(size(image(:,:,1)),PixNew3_Bilinear_floor1(:,2),PixNew3_Bilinear_floor1(:,1));
indBilinear2 = sub2ind(size(image(:,:,1)),PixNew3_Bilinear_floor2(:,2),PixNew3_Bilinear_floor2(:,1));
indBilinear3 = sub2ind(size(image(:,:,1)),PixNew3_Bilinear_floor3(:,2),PixNew3_Bilinear_floor3(:,1));
indBilinear4 = sub2ind(size(image(:,:,1)),PixNew3_Bilinear_floor4(:,2),PixNew3_Bilinear_floor4(:,1));

coeff1 = (1 - xyOrig_frac_in(:,2)).*(1 - xyOrig_frac_in(:,1));
coeff2 = (1 - xyOrig_frac_in(:,2)).*xyOrig_frac_in(:,1);
coeff3 = xyOrig_frac_in(:,2).*(1 - xyOrig_frac_in(:,1));
coeff4 = xyOrig_frac_in(:,2).*xyOrig_frac_in(:,1);

[coeff1, ~] = Num2Fix_unsigned(coeff1, 1,9);
[coeff2, ~] = Num2Fix_unsigned(coeff2, 0,9);
[coeff3, ~] = Num2Fix_unsigned(coeff3, 0,9);
[coeff4, ~] = Num2Fix_unsigned(coeff4, 0,9);




% if ~((max(xyOrig_int_in(:,2)) == min(xyOrig_int_in(:,2))) && (min(xyOrig_int_in(:,2)) == y_in(1)))
%     ccccctttt = ccccctttt + 1;
% end

% % idIn = xyOrig(:,1) > 0 & xyOrig(:,1) <= nc & xyOrig(:,2) > 0 & xyOrig(:,2) <= nr;
% % %                                     idIn = flag_in;
% % foundIndIn = foundInd(idIn);
% % foundIndIn = foundInd(flag_in);
%                                     origInd = sub2ind([nr,nc],xyOrig(idIn,2),xyOrig(idIn,1));
%                                     imgRect([foundIndIn;numPix + foundIndIn;2*numPix + foundIndIn]) = imgOrig([origInd;numPix + origInd;2*numPix + origInd]);
%
cof1 = Num2Fix_unsigned(coeff1.*double(imgR(indBilinear1)), 8,9);
cof2 = Num2Fix_unsigned(coeff2.*double(imgR(indBilinear2)), 8,9);
cof3 = Num2Fix_unsigned(coeff3.*double(imgR(indBilinear3)), 8,9);
cof4 = Num2Fix_unsigned(coeff4.*double(imgR(indBilinear4)), 8,9);

cof5 = Num2Fix_unsigned(coeff1.*double(imgG(indBilinear1)), 8,9);
cof6 = Num2Fix_unsigned(coeff2.*double(imgG(indBilinear2)), 8,9);
cof7 = Num2Fix_unsigned(coeff3.*double(imgG(indBilinear3)), 8,9);
cof8 = Num2Fix_unsigned(coeff4.*double(imgG(indBilinear4)), 8,9);

cof9 = Num2Fix_unsigned(coeff1.*double(imgB(indBilinear1)), 8,9);
cof10 = Num2Fix_unsigned(coeff2.*double(imgB(indBilinear2)), 8,9);
cof11 = Num2Fix_unsigned(coeff3.*double(imgB(indBilinear3)), 8,9);
cof12 = Num2Fix_unsigned(coeff4.*double(imgR(indBilinear4)), 8,9);

imgRect(foundIndIn) =  floor(cof1 + cof2 + cof3 + cof4);
imgRect(numPix + foundIndIn) = floor(cof5 + cof6 + cof7 + cof8);
imgRect(2*numPix + foundIndIn) = floor(cof9 + cof10 + cof11 + cof12);


% imgRect(numPix + foundIndIn) = coeff1.*double(imgG(indBilinear1)) + coeff2.*double(imgG(indBilinear2)) + coeff3.*double(imgG(indBilinear3)) + coeff4.*double(imgG(indBilinear4));
% imgRect(2*numPix + foundIndIn) = coeff1.*double(imgB(indBilinear1)) + coeff2.*double(imgB(indBilinear2)) + coeff3.*double(imgB(indBilinear3)) + coeff4.*double(imgB(indBilinear4));

end


function [fixPt, aa] = Num2Fix_unsigned(num, intLen,fracLen)
flag = num>=0;
% a = fi(num, 0, intLen+fracLen, fracLen,'RoundingMethod','Floor');
aa = fi(num(flag), 0, intLen+fracLen, fracLen,'RoundingMethod','Floor');
bb = fi(num(~flag), 0, intLen+fracLen, fracLen,'RoundingMethod','Floor');
% % a = fi(num, 0, intLen+fracLen, fracLen);
% fixPt = a.data;
fixPt = zeros(size(num));
fixPt(flag) = aa.data;
fixPt(~flag) = bb.data;
end