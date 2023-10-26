function corners = Detect(img0)


[X, dir, scores, middle] = extract_edgelets(img0, 3, 15);

% subplot(1,2,1);
% %figure,imshow(img0);hold on;plot(X(1,:), X(2,:),'.r');
% imshow(img0);hold on; plot(X(1,:), X(2,:),'.b');

P = size(X,2);
Q = 200;

di_mat = ones(P, 1) * 999999;
order_mat = zeros(Q, 1);
[a, b] = max( scores );
order_mat(1) = b;

for j = 2 : Q
   for i = 1:P
      diff = dir(:, i) - dir( :, order_mat( j -1) );
      diff_2 = diff' * diff;
      if di_mat(i) > diff_2
          di_mat(i) = diff_2;
      end 
   end
    sd_mat = scores.*di_mat;
    [c, d] = max(sd_mat);
    order_mat(j) = d;
end

filterX = X(:, order_mat);
corners = filterX';
% subplot(1,2,2);
% %figure,imshow(img0);hold on;plot(filterX(1,:), filterX(2,:),'.r');
% imshow(img0);hold on;plot(filterX(1,:), filterX(2,:),'.b');
end