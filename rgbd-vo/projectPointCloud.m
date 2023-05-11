function [img_coords, valid_mask, ind] = projectPointCloud(pointclouds, K, height, width)
% project the pointcloud into image coorindates
%
% INPUT:
%   pointclouds: a matrix of size [num_points, 3]
%   K: intrinsic camera parameters
%   height: image height
%   width: image width
%
% OUTPUT:
%   img_coords: a matrix of size [num_pixels, 2]
%   valid_mask: binary mask for valid pixels

margin = 5; 2; 5; 10; 5; 2;

img_coords = bsxfun(@rdivide, pointclouds, pointclouds(:, 3));
img_coords = img_coords * K';
img_coords = img_coords(:, 1:2);

% valid_mask = img_coords(:, 1) > 0 & img_coords(:, 1) <= width & img_coords(:, 2) > 0 & img_coords(:, 2) <= height;
valid_mask = img_coords(:, 1) > margin & img_coords(:, 1) <= width - margin & img_coords(:, 2) > margin & img_coords(:, 2) <= height - margin;
img_coords = img_coords(valid_mask, :);

ind = sub2ind([height, width], ceil(img_coords(:,2)), ceil(img_coords(:,1)));
end