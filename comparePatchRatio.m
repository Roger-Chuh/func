function comparePatchRatio()

% close all;

host_uvs = [375 379.781 377.391 364.746 369.873 376.368 377.736 379.104 373.632 376.023 372.264 368.505 378.759 371.241 380.127 372.609
    184 173.746 178.873 179.219 181.609 187.759 191.518 195.276 180.241 175.114 176.482 177.851 182.632 185.368 186.391 189.127];
target_uvs = [588.87 585.155 587.056 596.792 592.822 587.771 586.649 585.502 589.943 588.077 590.991 593.897 586.009 591.721 584.937 590.596
  122.8 130.327 126.592 126.147 124.471 120.058 117.325 114.601 125.549 129.332 128.305 127.243 123.859 121.707 121.133  118.95];
figure,imshow(ones(480, 640)); hold on;plot(host_uvs(1,:), host_uvs(2,:), '.r'); hold on;plot(target_uvs(1,:), target_uvs(2,:), '.b');



len1_host = norm(host_uvs(:,2) - host_uvs(:,4));
len2_host = norm((host_uvs(:,2) + host_uvs(:,4)) / 2 - host_uvs(:,8));

len1_target = norm(target_uvs(:,2) - target_uvs(:,4));

  if (abs(target_uvs(1, 2) - target_uvs(1, 4)) < 0.0001) 
    len2_target = abs(target_uvs(1, 8) - target_uvs(1, 2));
   else 
    slope = (target_uvs(2, 2) - target_uvs(2, 4)) / (target_uvs(1, 2) - target_uvs(1, 4));
    line_para = [slope, -1, -slope * target_uvs(1, 2) + target_uvs(2, 2)];
    line_para = line_para./norm(line_para(1:2));
   
    err2 = dot(line_para, [target_uvs(:, 2); 1]');
    err4 = dot(line_para, [target_uvs(:, 4); 1]');
    
    len2_target = abs(dot(line_para,[target_uvs(1, 8), target_uvs(2, 8), 1]));
    
    slope2 = (target_uvs(2, 1) - target_uvs(2, 8)) / (target_uvs(1, 1) - target_uvs(1, 8));
    line_para2 = [slope2, -1, -slope2 * target_uvs(1, 1) + target_uvs(2, 1)];
    line_para2 = line_para2./norm(line_para2(1:2));
    
    err1 = dot(line_para2, [target_uvs(:, 1); 1]');
    err8 = dot(line_para2, [target_uvs(:, 8); 1]');
    
    
    intersection = cross(line_para, line_para2);
    intersection = intersection ./ intersection(3);
    
    len2_check = norm(target_uvs(:,8) - intersection(1:2)');
    
    len2_err = len2_check - len2_target;
    
  end
plot(intersection(1), intersection(2),'og');

depth_scale_tmp = ((len1_target * len2_target) / (len1_host * len2_host))


end