function [y, x] = spline_(X, params)

	% Initialize x and y with 0's.
	x = 0; y = 0;

	% Number of nodes.
% 	n = size(X);
    n = length(X);

	% Calculate each subspline.
	for i = 1:n-1
		xi = [X(i):0.1:X(i + 1)]';
		hi = xi - X(i);
		Hi = [ones(size(xi)), hi, hi .^ 2, hi .^ 3];
		x = [x; xi];
		y = [y; Hi * params(i, :)'];
    end

	% Remove 0's from the beginning.
	x = x(2:end); y = y(2:end);

end