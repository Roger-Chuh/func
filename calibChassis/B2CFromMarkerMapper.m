function B2CFromMarkerMapper

close all;

inputDir = 'D:\Work\marker_mapper1.0.15\20201228_11633\video\L';

log = load(fullfile(inputDir, 'test_001.log'));

intrMat = [  5.19278503e+02, 0., 6.64875183e+02; 0., 5.19322693e+02,       3.93181030e+02; 0., 0., 1.];

dist = [-4.05976782e-03, 3.78094008e-03, -4.94747190e-04,       3.95255513e-04, -1.58973224e-03]';

dirInfo = dir(fullfile(inputDir, '*.png'));


for i = 1 : size(log,1)
    xList(i,:) = log(i,8:2:end);
    yList(i,:) = log(i,9:2:end);
end
for i = 2 : length(dirInfo)
    img = imread(fullfile(inputDir, dirInfo(i).name));
   
    corner = cornerfinder([xList(i-1,:); yList(i-1,:)],double(img), 2, 2);
    cornerX(:,i-1) = corner(1,:)';
    cornerY(:,i-1) = corner(2,:)';

    figure,imshow(img);hold on;plot(xList(i-1,:), yList(i-1,:),'.r');
    plot(corner(1,:), corner(2,:),'.g');
%     plot(xList(i-1,1), yList(i-1,1),'*b');
    drawnow;
    
end


end