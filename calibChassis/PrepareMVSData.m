function PrepareMVSData(inputDir0, probPath, kfDir)

doHist = 0;
%  PrepareMVSData('C:\Users\rongjiezhu\Pictures\Camera Roll\zedRect\office010\image_1', 'C:\Users\rongjiezhu\Pictures\Camera Roll\zedRect\office010\depthOpt')

% PrepareMVSData('C:\Users\rongjiezhu\Pictures\Camera Roll\zedRect\office011\data', 'C:\Users\rongjiezhu\Pictures\Camera Roll\zedRect\office011\depthOpt','C:\Users\rongjiezhu\Pictures\Camera Roll\zedRect\office011\kf')



if 0
    rotMat = rodrigues([0.02 0.03 -0.05]);
    q = rotMat2quatern(rotMat);
    quatern2rotMat(q) - rotMat
end






try
    dirInfoKey = dir(fullfile(kfDir,'key*.mat'));
    keyframeInfoStack = {};
    for u = 1:length(dirInfoKey)
        load(fullfile(kfDir, dirInfoKey(u).name));
        keyframeInfoStack = [keyframeInfoStack; keyframeInfo(1:1:end,:)];
    end
catch
    sdgjh = 1;
end

clear 'intrMatL';

scale =  0.5;

inputDir = inputDir0;
dirInfo = dir(fullfile(inputDir,'depth_*.mat'));

load(fullfile(inputDir, dirInfo(1).name));
intrMat = scale.*intrMatL0;
intrMat(3,3) = 1;
 

% baseline = 0.11964611;
dispRng = 80;


cnt = 1;


tumData = 1;

if ~tumData
    fid1 = fopen(fullfile(probPath,'cameras.txt'),'w');%????
    fprintf(fid1,sprintf('%d\n\n',length(dirInfo)));
else
    scale = 1;
    baseTime = 1341845688;
    Hz = 30;
    
    fid1 = fopen(fullfile(probPath,'rgb.txt'),'w');
    fprintf(fid1,sprintf('# rgb maps\n# a\n# b\n'));
    
    fid2 = fopen(fullfile(probPath,'depth.txt'),'w');
    fprintf(fid2,sprintf('# depth maps\n# a\n# b\n'));
    
    fid3 = fopen(fullfile(probPath,'groundtruth.txt'),'w');
    fprintf(fid3,sprintf('# gt\n# a\n# b\n'));
    
    MakeDirIfMissing(fullfile(probPath,'rgb'));
    MakeDirIfMissing(fullfile(probPath,'depth'));
end


for i = 1 : length(dirInfo) %size(keyframeInfoStack,1) %
    if 0
        imgId = keyframeInfoStack{i,1}(1);
    else
        imgId = i;
    end
    if 0
        load(fullfile(inputDir, dirInfo(imgId).name));
    else
        load(fullfile(inputDir, sprintf('depth_%06d.mat',imgId-1)));
    end
    
    if doHist
        rgb = histeq(rgb);
        rgbR = histeq(rgbR);
    end
    
    
    curTime = baseTime + 1/Hz*cnt;
    colEnd = floor(size(rgb,2)/2)*2;
    rowEnd = floor(size(rgb,1)/2)*2;
    rgb = imresize(rgb(1:rowEnd, 1:colEnd), scale);
    rgbR = imresize(rgbR(1:rowEnd, 1:colEnd), scale);
    img = cat(3, rgb,rgb,rgb);
    
    
    
    if tumData
        
        depth = 5000.*depthInit(1:rowEnd, 1:colEnd);
        depth(isnan(depth)) = 0;
        depth = uint16(depth);
        
        qq = rotMat2quatern(pose(1:3,1:3));
        q = qq([2 3 4 1]);
        t = pose(1:3,4)';
        
        imwrite(img, fullfile(probPath, 'rgb',sprintf('%0.6f.png',curTime)));
        imwrite(depth,fullfile(probPath, 'depth',sprintf('%0.6f.png',curTime)),'png','bitdepth',16);
        fprintf(fid1,sprintf('%0.6f rgb/%0.6f.png\n',curTime, curTime));
        fprintf(fid2,sprintf('%0.6f depth/%0.6f.png\n',curTime, curTime));
        fprintf(fid3,sprintf('%0.4f %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f\n',curTime, t, q));
        
    else
        if 1
            try
                dispMat = 1.*imresize(uint16(round(intrMat(1,1)*baseline./depthInit(1:rowEnd, 1:colEnd))), scale);
            catch
                sgkhj = 1;
            end
        else
            dispMat = disparity(rgb, rgbR,'DisparityRange', [0 dispRng],'BlockSize', 15);
            dispMat(dispMat < 10) = 0;
            dispMat(dispMat > dispRng-5) = 0;
            depthMat = intrMat(1,1).*baseline./dispMat;
            depthMat(isinf(depthMat)) = 0;
            dispMat = uint16(dispMat);
        end
        
        dispMat(isnan(dispMat)) = 0;
        r = pose(1:3,1:3);
        
        
        t = 100*pose(1:3,4)';  % m->cm
        rt = [r;t];
        krt = [intrMat;rt];
        
        if 1
            imwrite(img, fullfile(probPath, sprintf('test%04d.jpg',cnt-1)));
        else
            imwrite(img, fullfile(probPath, sprintf('test%04d.jpg',imgId-1)));
        end
        
        if 0
            dispMatT = dispMat';
            calibRomFid = fopen(fullfile(probPath, sprintf('depth_%04d.bin', cnt-1)),'w');
            fwrite(calibRomFid, dispMatT(:),'uint16');
            fclose(calibRomFid);
        end
        
        
        if 0
            imageT =  (rgb)';
            calibRomFid = fopen(fullfile(probPath, sprintf('img_%04d.bin', cnt)),'w');
            fwrite(calibRomFid, im2double(imageT(:)),'double');
            fclose(calibRomFid);
            
            depthMatT = double(1.*double(1.*depthMat'));
            calibRomFid = fopen(fullfile(probPath, sprintf('dep_%04d.bin', cnt)),'w');
            fwrite(calibRomFid, (depthMatT(:)),'double');
            fclose(calibRomFid);
            
            poseT = inv(pose)';
            poseT(4,1:3) = poseT(4,1:3)./1;
            calibRomFid = fopen(fullfile(probPath, sprintf('pose_%04d.bin', cnt)),'w');
            fwrite(calibRomFid, poseT(:),'single');
            fclose(calibRomFid);
            
            intrMatT = intrMat';
            calibRomFid = fopen(fullfile(probPath, sprintf('intrMat_%04d.bin', 1)),'w');
            fwrite(calibRomFid, intrMatT(:),'double');
            fclose(calibRomFid);
            %         cnt = cnt+1;
        end
        
        
        for ji = 1 : size(krt,1)
            fprintf(fid1,sprintf('%0.7f    %0.7f    %0.7f \n',krt(ji,1), krt(ji,2), krt(ji,3)));
            %         fprintf(fid1,'\n');
        end
        fprintf(fid1,'\n\n');
    end
    cnt = cnt + 1;
end
if ~tumData
    fprintf(fid1,'\n');
end
fclose(fid1);
% fprintf(fid2,'\n');
fclose(fid2);
% fprintf(fid3,'\n');
fclose(fid3);

end