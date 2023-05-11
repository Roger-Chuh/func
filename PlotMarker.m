function PlotMarker()
inputDir = 'C:\vs\test\ConsoleApplication2\ConsoleApplication2\3';

imgInfo = dir(fullfile(inputDir, '*.jpg'));
cornerInfo = dir(fullfile(inputDir, '*.txt'));

scale = 0.25;

for i = 1:length(imgInfo)
   corner = load(fullfile(inputDir, cornerInfo(i).name));
   id1 = find(corner(:,1)<10000);
   id2 = find(corner(:,1)>10000);
   corner = corner(:,2:end);
   corner = corner + 1;
   
   cnr1 = [corner(id1,1:2); corner(id1,3:4); corner(id1,5:6); corner(id1,7:8)];
   cnr2 = corner(id2,1:2);
   img = imresize(imread(fullfile(inputDir, imgInfo(i).name)),scale);
   figure,imshow(img);hold on;plot(cnr1(:,1), cnr1(:,2),'.r') ;plot(cnr2(:,1), cnr2(:,2),'.g') ;
end

end