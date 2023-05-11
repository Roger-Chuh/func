function showTrackingFrame()

inputDir = 'G:\matlab\data\dso\1';

inputDir1 = strcat(inputDir,'\imgs_different_radius');
inputDir2 = strcat(inputDir,'\imgs_same_radius');

dirInfo1 = dir(fullfile(inputDir1, '*.png'));
dirInfo2 = dir(fullfile(inputDir2, '*.png'));


timestamp1 = [];
timestamp2 = [];

timestamp1 = extractTimeInfo(dirInfo1);

timestamp2 = extractTimeInfo(dirInfo2);

[comm, idx1, idx2] = intersect(timestamp1, timestamp2);

for i = 1 : length(comm)
   img1 = imread(fullfile(inputDir1, dirInfo1(idx1(i)).name)); 
   img2 = imread(fullfile(inputDir2, dirInfo2(idx2(i)).name)); 
   figure(1),imshow([img1 img2]);
   drawnow;
end

end
function timestamp1 = extractTimeInfo(dirInfo1)
timestamp1 = [];
for i = 1 : length(dirInfo1)
    name = dirInfo1(i).name;
    timestamp1 = [timestamp1; str2double(name(7:end-4))];
    
end



end