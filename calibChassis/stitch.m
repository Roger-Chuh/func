function [img_bb, pad_up] = stitch(img_i, img_b, rt)
    hei_i = size(img_i,1);
    wid_i = size(img_i,2);
    hei_b = size(img_b,1);
    wid_b = size(img_b,2);
    if 0
        cimg_b = corner_detector(img_b);
        [x_b, y_b, rmax_b] = anms(cimg_b, 1500);
        [descs_b] = feat_desc_geoblur(img_b, x_b, y_b);
        
        cimg_i = corner_detector(img_i);
        [x_i, y_i, rmax_i] = anms(cimg_i, 1500);
        [descs_i] = feat_desc_geoblur(img_i, x_i, y_i);
        [match_i] = feat_match(descs_i, descs_b);
        
        % extract matched points
        xx_i = x_i(match_i ~= -1);
        yy_i = y_i(match_i ~= -1);
        xx_b = x_b(match_i(match_i ~= -1));
        yy_b = y_b(match_i(match_i ~= -1));
        
        % reject outliers
        [H, inlier_ind] = ransac_est_homography(xx_i, yy_i, xx_b, yy_b, 10);
    end
    % compute d2b
    img_i_dist = dist2border(img_i);
    img_b_dist = dist2border(img_b);
    
    % mapped coordinates of upperleft, upperright, bottomleft, bottomright
    % for img_i
    if 0
        ul = H*[1 1 1]'; ul = ul/ul(end);
        ur = H*[wid_i, 1, 1]'; ur = ur/ur(end);
        bl = H*[1, hei_i, 1]'; bl = bl/bl(end);
        br = H*[wid_i, hei_i, 1]'; br = br/br(end);
    else
        ul = (rt(1:3,1:3)*[1 1 1]' + repmat(rt(1:3,4),1,1));ul = ul/ul(end);
        ur = (rt(1:3,1:3)*[wid_i, 1, 1]'+ repmat(rt(1:3,4),1,1)); ur = ur/ur(end);
        bl = (rt(1:3,1:3)*[1, hei_i, 1]'+ repmat(rt(1:3,4),1,1)); bl = bl/bl(end);
        br = (rt(1:3,1:3)*[wid_i, hei_i, 1]'+ repmat(rt(1:3,4),1,1)); br = br/br(end);
    end
%     tic
    
    % find out how much padding we need to make
    pad_up = 0; % update pad_up and pad_left inorder to do an offset mapping 
    pad_left = 0;
    pad_down = 0;
    pad_right = 0;
    if max(br(1),ur(1)) > wid_b
        pad_right = round(max(br(1),ur(1))-wid_b+30);
        img_b = padarray(img_b, [0, pad_right], 'post');
    end
    if max(br(2), bl(2)) > hei_b
        pad_down = round(max(br(2), bl(2))-hei_b+30);
        img_b = padarray(img_b, [pad_down, 0], 'post');
    end
    if min(ul(1), bl(1)) <= 0 
        pad_left = round(-min(ul(1), bl(1))+30);
        img_b = padarray(img_b, [0, pad_left], 'pre');
    end
    if min(ul(2), ur(2)) <= 0
        pad_up = round(-min(ul(2), ur(2)) + 30); 0;
        img_b = padarray(img_b, [pad_up, 0], 'pre');
    end
    if 0
        H_inv = inv(H);
    else
        rt_inv = inv(rt);
    end
    
    % - Mapping coordinates from img_b to img_i to retrieve pixels
    % - Implemented distance to border blending
    % - Vecterized
    
    [y_b, x_b] = meshgrid(round(pad_up+min(ul(2), ur(2))):round(pad_up+max(bl(2), br(2))), ...
        round(pad_left+min(ul(1), bl(1))):round(pad_left+max(br(1),ur(1))));
    y_b = y_b(:); x_b = x_b(:);
    if 0
        xy = H_inv*[x_b - pad_left, y_b - pad_up, ones(size(x_b,1),1)]';
    else
         xy =  (rt_inv(1:3,1:3)*[x_b - pad_left, y_b - pad_up, ones(size(x_b,1),1)]'+ repmat(rt_inv(1:3,4),1,length(x_b)));
    end
    x_i = int64(xy(1,:)'./xy(3,:)'); y_i = int64(xy(2,:)'./xy(3,:)');
    
    % Blend img_b and img_i pixels according to dist to boundary
    % Only the coordinates that are possible to map to img_i are considerred.
    indices = x_i > 0 & x_i <= size(img_i, 2) & y_i > 0 & y_i <= size(img_i, 1) & y_b-pad_up > 0 & y_b - pad_up <= size(img_b_dist,1) & x_b - pad_left > 0 & x_b - pad_left <= size(img_b_dist,2);
    idx_i = (x_i(indices)-1)*size(img_i_dist,1) + y_i(indices);
    idx_b = (x_b(indices)-pad_left-1)*size(img_b_dist,1) + y_b(indices)-pad_up;
    p = img_i_dist(idx_i)./(img_i_dist(idx_i) + img_b_dist(idx_b));
    p(isnan(p)) = 0;
    % map rgb channels
    img_b((x_b(indices)-1)*size(img_b,1)+y_b(indices)) = p.*img_i((x_i(indices)-1)*size(img_i,1) + y_i(indices)) + ... 
        (1-p).*img_b((x_b(indices)-1)*size(img_b,1) + y_b(indices));
    img_b((x_b(indices)-1)*size(img_b,1)+y_b(indices) + size(img_b,1)*size(img_b,2) ) = p.*img_i((x_i(indices)-1)*size(img_i,1) + y_i(indices) + size(img_i,1)*size(img_i,2)) + ... 
        (1-p).*img_b((x_b(indices)-1)*size(img_b,1) + y_b(indices) + size(img_b,1)*size(img_b,2));
    img_b((x_b(indices)-1)*size(img_b,1)+y_b(indices) + size(img_b,1)*size(img_b,2)*2 ) = p.*img_i((x_i(indices)-1)*size(img_i,1) + y_i(indices) + size(img_i,1)*size(img_i,2)*2) + ... 
        (1-p).*img_b((x_b(indices)-1)*size(img_b,1) + y_b(indices) + size(img_b,1)*size(img_b,2)*2);
    
    % if it goes beyond border of img_b, just copy back the pixels from img_i
    indices = x_i > 0 & x_i <= size(img_i, 2) & y_i > 0 & y_i <= size(img_i, 1) & (y_b-pad_up <= 0 | y_b - pad_up > size(img_b_dist,1) | x_b - pad_left <= 0 | x_b - pad_left > size(img_b_dist,2));
    % rgb channels
    try
        img_b((x_b(indices)-1)*size(img_b,1)+y_b(indices)) = img_i((x_i(indices)-1)*size(img_i,1) + y_i(indices));
    catch
        sgdkhj = 1;
    end
    img_b((x_b(indices)-1)*size(img_b,1)+y_b(indices) + size(img_b,1)*size(img_b,2) ) = img_i((x_i(indices)-1)*size(img_i,1) + y_i(indices) + size(img_i,1)*size(img_i,2));
    img_b((x_b(indices)-1)*size(img_b,1)+y_b(indices) + size(img_b,1)*size(img_b,2)*2 ) = img_i((x_i(indices)-1)*size(img_i,1) + y_i(indices) + size(img_i,1)*size(img_i,2)*2);
    
    img_bb = uint8(img_b);
%     toc
end
function [img_dist] = dist2border(img)
    if size(img,3) > 1
        img = rgb2gray(img);
    end
    img_dist = (img == 0);
    img_dist(1:size(img,1), 1) = 1;
    img_dist(1:size(img,1), size(img,2)) = 1;
    img_dist(1, 1:size(img,2)) = 1;
    img_dist(size(img,1), 1:size(img,2)) = 1;
    img_dist = bwdist(img_dist, 'chessboard');
end