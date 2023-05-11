% inputDir = 'D:\temp5\20171010\1\1';
function changeExtensionByEvalBin(inputDir, varargin)

global cfg
if (nargin == 1)
    nameStr = 'imgR';
elseif (nargin == 2)
    nameStr = varargin{1};
else
    error('Too many input arguments');
end

isBin = 0;

if ~isBin
    dirInfo = dir(fullfile(inputDir,'*.png'));
else
    dirInfo = dir(fullfile(inputDir,'*.bin'));
end
if length(dirInfo) == 0
    dirInfo = dir(fullfile(inputDir,'*.bmp'));
end
if length(dirInfo) == 0
    dirInfo = dir(fullfile(inputDir,'*.jpg'));
end
if 	0
    dircell=struct2cell(dirInfo);
    if 0
        aaAa = dircell(2,:)';
        sss = datetime(aaAa,'Format','dd-MMM-y HH:mm:ss');
        [B,I] = sort(sss);
    else
        [~, ind1] = sort([dirInfo(1:length(dirInfo)).datenum], 'ascend');
        I = ind1;
    end
    tempInfo = cell(length(I),1);
    for fg = 1:length(I)
        tempInfo{fg,1} = dirInfo(I(fg)).name;
    end
    for fgg = 1:length(I)
        dirInfo(fgg).name = tempInfo{fgg};
    end
end


frameNum = length(dirInfo); 1/2;  %/2;
workDir = pwd;
cd (inputDir);
ct = 5725;
for i = 1:    frameNum  %length(dirInfo)
    
    if 0
        img = imread(fullfile(inputDir,dirInfo(i).name));
        %    name = str2double(dirInfo(i).name(5:end-4));
        %    name = name.*100; name = round(name);
        %    name = name*10000000 + 1000000000000000000;
        if 0
            imwrite(img,strcat(inputDir,'\',sprintf('%19.0f',name),'.png'));
        else
%             imwrite(img, sprintf(strcat(nameStr, '_%04d.png'), i));
            imwrite(img, sprintf('a%06d.png', i-1));
            
        end
    elseif 1
        if isBin
            eval(['!rename' 32 dirInfo(i).name 32 strcat(dirInfo(i).name(end-4),'_',dirInfo(i).name(1:end-6),'.bin')]);
        else
            if 0
                eval(['!rename' 32 strcat(dirInfo(i).name(end-4),'_',dirInfo(i).name(1:end-6),'.jpg') 32 sprintf(strcat(nameStr, '_%04d.jpg'), i)]);
            else
                eval(['!rename' 32 dirInfo(i).name 32 sprintf(strcat(nameStr, '_%04d.jpg'), i)]);
            end
        end
    elseif 0
        eval(['!rename' 32 dirInfo(i).name 32 sprintf('a%06d.png', i-1)]);
    else
        eval(['!rename' 32 dirInfo(i).name 32 sprintf('%06d.png', ct)]);
    end
    ct = ct+1;
    %    eval(['!rename' 32 dirInfo(i).name 32 sprintf('t1_%04d.png', i)]);
    %     eval(['!rename' 32 dirInfo(i).name 32 sprintf(strcat(nameStr, '_%04d.png'), i)]);
    
    % % % % % %     eval(['!rename' 32 dirInfo(i).name 32 sprintf( '%06d.png', i-1)])0;
    
    %    eval(['!rename' 32 dirInfo(i + frameNum).name 32 sprintf('imgR_%04d.png', i)]);
    
    %     eval(['!rename' 32 dirInfo(i).name 32 sprintf('rectR_%04d.png', i)]);
    
    
end
cd (workDir)
end