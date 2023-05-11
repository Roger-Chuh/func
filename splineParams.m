function params = splineParams(X, y)
	
	% Quantity of nodes.
	n = length(X);
	
	l = zeros(n, 1);
	m = zeros(n, 1);
	z = zeros(n, 1);

	% First parameter vector.
	a = y(1:end-1);

	% Second parameter incomplete vector.	
	b = zeros(n - 1, 1);

	% Third parameter incomplete vector.
	c = zeros(n, 1);

	% Fourth parameter incomplete vector.
	d = zeros(n - 1, 1);

	h = X(2:end) - X(1:end-1);

	a1 = y(3:end) - y(2:end-1);
	a2 = y(2:end-1) - y(1:end-2);
	h1 = h(2:end);
	h2 = h(1:end-1);

	alpha = 3 .* (a1./h1 - a2./h2);

	for i = 2:(n - 1)
		l(i) = 2 * (X(i + 1) - X(i - 1)) - h(i - 1) * m(i - 1);
		m(i) = h(i) / l(i);
		z(i) = (alpha(i - 1) - h(i - 1) * z(i - 1))/l(i);
    end

	l(n) = 1; 

	for i = (n-1):-1:1
		c(i) = z(i) - m(i) * c(i + 1);
		b(i) = (y(i + 1)-a(i))/h(i) - h(i)*(c(i + 1) + 2 * c(i))/3;
		d(i) = (c(i + 1) - c(i)) / (3 * h(i));
    end

	c = c(1:end-1);

	params = [a, b, c, d];

end