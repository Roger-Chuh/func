function landmark_use = showDepthFilterLandmark(varargin)

if nargin < 1
   filename =  'landmark_statistics.txt';
elseif nargin == 1
      filename = varargin{1};
else
end
inputDir = 'G:\matlab\data\direct\sim';

pc = load(fullfile(inputDir,filename));
search_fail = pc(:, end-1:end);
pc = pc(1:end,:);
idepth = pc(:,2);
landmark = pc(:,3:5);
sigma = pc(:,6);
goodness = pc(:,7);
[~,depths_] = NormalizeVector(landmark);
%  figure(4000),pcshow(landmark(goodness == 1 & depths_<2 & sigma < 1,:),'MarkerSize', 100);

timing = load(fullfile(inputDir, 'time_statistics.txt'));
landmark_use = landmark(goodness == 1 & depths_<15.5 & sigma <  0.5 & idepth > -0.8,:);
figure,
% subplot(1,2,1),
pcshow(landmark_use,'MarkerSize', 100);title(sprintf(strcat(filename, ', %d frames, %d landmarks'), length(timing), length(sigma)));
% subplot(1,2,2);hist(timing(timing(:,2) < 300,2),200);grid on;title(sprintf('%d frames, %d landmarks', length(timing), length(sigma)));
% subplot(1,2,2);hist(sigma,200);
fprintf(sprintf('sigma num: %d\n', length(sigma)));

figure,subplot(1,2,1);hist(search_fail, 100);
subplot(1,2,2),hist(search_fail(:,2)./search_fail(:,1), 100);title('fail / search ratio');
end