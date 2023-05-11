function read_board_depth(srcdirname, w0, h0)
global cfg


if cfg.isRotate
    w = min([h0 w0]);
    h = max([h0 w0]);
else
    h = min([h0 w0]);
    w = max([h0 w0]);
end
dstdirname = srcdirname;
% srcdirname = 'D:\Temp\20201112\record';
% dstdirname = 'D:\Temp\20201112\record\1';
% 

% srcdirname = 'D:\Temp\20201114\raw\raw';
% dstdirname = 'D:\Temp\20201114\raw\raw\1';

if ~exist(dstdirname, 'dir')
    mkdir(dstdirname);
end

subdirs = dir(fullfile(srcdirname,'depth_src*.raw'));
if length(subdirs) == 0
    subdirs = dir(fullfile(srcdirname,'src*.raw'));
end
dispInfo = dir(fullfile(srcdirname,'disp*.raw'));
[~, ind11] = sort([dispInfo(1:end).datenum], 'ascend');
tempInfo1 = cell(length(ind11),1);
for fg = 1:length(ind11)
    tempInfo1{fg,1} = dispInfo(ind11(fg)).name;
end
for fgg = 1:length(ind11)
    dispInfo(fgg).name = tempInfo1{fgg};
end
cfg.dispInfo = dispInfo;
cfg.inputDir = srcdirname;


dirInfo = subdirs;
[~, ind1] = sort([dirInfo(1:length(dirInfo)/2).datenum], 'ascend');
[~, ind2] = sort([dirInfo(length(dirInfo)/2+1:end).datenum], 'ascend');
ind2 = ind2 + length(dirInfo)/2;
I = [ind1 ind2];
tempInfo = cell(length(I),1);
for fg = 1:length(I)
    tempInfo{fg,1} = dirInfo(I(fg)).name;
end
for fgg = 1:length(I)
    dirInfo(fgg).name = tempInfo{fgg};
end
cfg.dirInfo = dirInfo;
subdirs = dirInfo;
for i = 1 : length(subdirs)/2
    if strcmp(subdirs(i).name, '.') || strcmp(subdirs(i).name, '..')
        continue
    end   
    subnameL = fullfile(srcdirname, subdirs(i).name);
    numstrL = subdirs(i).name(end-7:end-4);
    depth_imgL = getImg(subnameL, w, h);
    subnameR = fullfile(srcdirname, subdirs(i + length(subdirs)/2).name);
    numstrR = subdirs(i + length(subdirs)/2).name(end-7:end-4);
    depth_imgR = getImg(subnameR, w, h);
    
%     filename= fullfile(dstdirname, ['depth_',numstr,'.png']);
    imwrite(uint8(depth_imgL), fullfile(dstdirname,sprintf('imgL_%05d.png', i)));
    imwrite(uint8(depth_imgR), fullfile(dstdirname,sprintf('imgR_%05d.png', i)));
end
end
function depth_img = getImg(subname, w, h)
disp(subname);
fp=fopen(subname,'rb');
%     depth= fread(fp, 'uint16=>uint16');
depth= fread(fp, 'uint8');
if 0
    depth_img=reshape(depth(:),[w, h])';
else
%     depth_img=reshape(depth(w*h+1:end),[w, h])';
    depth_img=reshape(depth(1:w*h),[w, h])';
end

end