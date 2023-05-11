function LoopCalibration()
inputDir = 'C:\Users\Roger\Desktop\Depth_Test_Tool\demo-stereocapture\loop';
dirInfo = dir(fullfile(inputDir,'0*'));

loopNum = 100; 10; 100;
eachSampleNum = 4; 3;
TempTrialDirStack = {};
validId = {};
for i = 1 : length(dirInfo)
    
   tempDir = fullfile(inputDir, dirInfo(i).name); 
   tempImgDir = fullfile(tempDir,sprintf('%02d',i));
   MakeDirIfMissing(tempImgDir);
   dirInfoL = dir(fullfile(tempDir,'img_l','*.png'));
   dirInfoR = dir(fullfile(tempDir,'img_r','*.png'));
   imgNum = length(dirInfoL);
   TempTrialDir = {};
   TempTrialDirStack{i} = {};
   validId{i} = [];
%    TempTrialDirStack = [TempTrialDirStack; TempTrialDir];
   for j = 1 : loopNum
       randInd = randperm(imgNum);
       randId = randInd(1:eachSampleNum);
       tempTrialDir = fullfile(tempImgDir, sprintf('%03d',j));
%        TempTrialDirStack{i} = [TempTrialDirStack{i};tempTrialDir];
       
       MakeDirIfMissing(tempTrialDir);
       delete(fullfile(tempTrialDir,'*'));
       imgDirL = fullfile(tempTrialDir,'img_l');
       imgDirR = fullfile(tempTrialDir,'img_r');
       MakeDirIfMissing(imgDirL);
       MakeDirIfMissing(imgDirR);
       delete(fullfile(imgDirL,'*'));
       delete(fullfile(imgDirR,'*'));
       for k = 1 : eachSampleNum
          copyfile(fullfile(tempDir,'img_l',dirInfoL(randId(k)).name), imgDirL);
          copyfile(fullfile(tempDir,'img_r',dirInfoL(randId(k)).name), imgDirR);
       end
       cmdList = ['C:\vs\test\ConsoleApplication4\x64\Release\calib_tools.exe 640 400 1 1 ',tempTrialDir,' C:/vs/test/ConsoleApplication4/param'];
       system(cmdList);
       if length(dir(fullfile(tempTrialDir,'calib_param_320x200_01.eeprom'))) > 0
           validId{i} = [validId{i}; j];
           TempTrialDirStack{i} = [TempTrialDirStack{i};tempTrialDir];
       end
   end
    
end

binCnt = 30; 
for l = 1: length(TempTrialDirStack)
    intrVecBuf{l,1} = zeros(length(TempTrialDirStack{l}),8);
    distVecBuf{l,1} = zeros(length(TempTrialDirStack{l}),10);
    intrVecOldBuf{l,1} = zeros(length(TempTrialDirStack{l}),8);
    distVecOldBuf{l,1} = zeros(length(TempTrialDirStack{l}),10);
    rtVecBuf{l,1} = zeros(length(TempTrialDirStack{l}),6);
    rtVecOldBuf{l,1} = zeros(length(TempTrialDirStack{l}),6);
    validList = validId{l};
    for ll = 1 : length(TempTrialDirStack{l})%length(validList) % length(TempTrialDirStack{l})
%         [intrVecBuf{l,1}(ll,:), distVecBuf{l,1}(ll,:), intrVecOldBuf{l,1}(ll,:), distVecOldBuf{l,1}(ll,:),  rtVecBuf{l,1}(ll,:), rtVecOldBuf{l,1}(ll,:)] = parseParam(TempTrialDirStack{l}{validList(ll)});
        [intrVecBuf{l,1}(ll,:), distVecBuf{l,1}(ll,:), intrVecOldBuf{l,1}(ll,:), distVecOldBuf{l,1}(ll,:),  rtVecBuf{l,1}(ll,:), rtVecOldBuf{l,1}(ll,:)] = parseParam(TempTrialDirStack{l}{(ll)});
    end
    
    %    figure,subplot(2,2,1),hist([intrVecBuf{l}(:,1) intrVecOldBuf{l}(:,1)]);title('fx');
    %    subplot(2,2,1),hist([intrVecBuf{l}(:,1) intrVecOldBuf{l}(:,2)]);title('fx');
    
    figure, subplot(2,2,1);hist(intrVecBuf{l,1}(:,[1 2 5 6]), binCnt);title('fx fy l r');
    subplot(2,2,2);hist(intrVecBuf{l,1}(:,[3 4 7 8]), binCnt);title('cx cy l r');
    subplot(2,2,3);hist(rtVecBuf{l,1}(:,[1:3]), binCnt);title('rot vec');
    subplot(2,2,4);hist(rtVecBuf{l,1}(:,[5:6]), binCnt);title('trans vec');
    
    figure, subplot(2,2,1);hist(intrVecOldBuf{l,1}(:,[1 2 5 6]), binCnt);title('fx fy l r');
    subplot(2,2,2);hist(intrVecOldBuf{l,1}(:,[3 4 7 8]), binCnt);title('cx cy l r');
    subplot(2,2,3);hist(rtVecOldBuf{l,1}(:,[1:3]), binCnt);title('rot vec');
    subplot(2,2,4);hist(rtVecOldBuf{l,1}(:,[5:6]), binCnt);title('trans vec');
end

end
function [intrVec, distVec, intrVecOld, distVecOld, rtVec, rtVecOld] = parseParam(tempTrialDir)
intrMat = ReadYaml(fullfile(tempTrialDir,'intrinsics.yml'));
intrMatOld = ReadYaml(fullfile(tempTrialDir,'intrinsics_old.yml'));

extrMat = ReadYaml(fullfile(tempTrialDir,'extrinsics.yml'));
extrMatOld = ReadYaml(fullfile(tempTrialDir,'extrinsics_old.yml'));

intrVec = [intrMat.M1.data(1,1) intrMat.M1.data(2,2) intrMat.M1.data(1,3) intrMat.M1.data(2,3) ...
    intrMat.M2.data(1,1) intrMat.M2.data(2,2) intrMat.M2.data(1,3) intrMat.M2.data(2,3)];
distVec = [intrMat.D1.data  intrMat.D2.data];

intrVecOld = [intrMatOld.M1.data(1,1) intrMatOld.M1.data(2,2) intrMatOld.M1.data(1,3) intrMatOld.M1.data(2,3) ...
    intrMatOld.M2.data(1,1) intrMatOld.M2.data(2,2) intrMatOld.M2.data(1,3) intrMatOld.M2.data(2,3)];
distVecOld = [intrMatOld.D1.data intrMatOld.D2.data];

rtVec = [rodrigues(extrMat.R.data)' extrMat.T.data];
rtVecOld = [rodrigues(extrMatOld.R.data)' extrMatOld.T.data];
end