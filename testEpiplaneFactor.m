function testEpiplaneFactor()

n = rand(3,1);
n = n./norm(n);
pt = 0.5*n;
n1 = rand(3,1);
n1 = n1./norm(n1);
n2 = cross(n,n1);
dir = n1;

pt=dir-dot(dir,n) * n;
a = (dir-pt)./norm(dir-pt);
[a-n]
dot(pt, n)



K = [300 0 320; 0 300 240; 0 0 1];

xyz = [10 20 1000;20 -20 1002;-10 50 1003;-10 -50 990];

pose{1,1} = eye(4);
for i = 1 : size(xyz,1)
    if i == 1
        pose{i,1} = eye(4);
    else
        pose{i,1} = [rodrigues(0.1*(rand(3,1)-0.5)) 1*(rand(3,1)-0.5); 0 0 0 1];
    end
    
    bearing = pose{i,1}(1:3,1:3)*xyz' + repmat(pose{i,1}(1:3,4), 1, size(xyz,1));
    [bearing_] = NormalizeVector(bearing');
    bearings{i,1} = bearing_;
    
end


opt_fids = 2:size(xyz,1);

pose_Twc = pose;
pose_err = {};
for i = opt_fids
   pose_err{i,1} = [rodrigues(0.2*(rand(3,1)-0.5)) 0.1*(rand(3,1)-0.5); 0 0 0 1];
   pose_Twc{i,1} = pose_err{i,1} * inv(pose{i,1}); 
end


host_fid = 1;



Twc_host = pose_Twc{host_fid,1};
Twc_host_gt = inv(pose{host_fid,1});
for pid = 1 : size(xyz,1)
   pt = xyz(pid,:);
   host_bearing = bearings{host_fid}(pid,:);
   assert(norm(pt./norm(pt) - host_bearing) < 0.000001);
   for target_fid = opt_fids
      target_bearing = bearings{target_fid}(pid,:);
      reproj =  pose{target_fid,1}(1:3,1:3)*pt' + pose{target_fid,1}(1:3,4);
      assert(norm(target_bearing' - reproj./norm(reproj)) < 0.000001);
      Twc_target = pose_Twc{target_fid,1};
      Twc_target_gt = inv(pose{target_fid,1});
      [err, d_err_d_host_pose, d_err_d_target_pose] = computeEpiplaneFacor(Twc_host, Twc_target, host_bearing, target_bearing);
      delta_target_pose = pose_err{target_fid,1};
      delta_target_pose_vec = [rodrigues(delta_target_pose(1:3,1:3));delta_target_pose(1:3,4)];
      err_comp = d_err_d_target_pose * delta_target_pose_vec;
      err - err_comp
%       d_err_d_host_pose * 
   end
       
       
  
    
    
end




end
function [err, d_err_d_host_pose, d_err_d_target_pose] = computeEpiplaneFacor(Twc_host, Twc_target, host_bearing, target_bearing)
T_th = inv(Twc_target) * Twc_host;
plane = SkewSymMat(T_th(1:3,4)) * T_th(1:3,1:3) * host_bearing';
plane_norm = plane./norm(plane);

pt_in_plane = target_bearing' - (target_bearing * plane_norm) .*plane_norm;
plane_norm_check = (pt_in_plane - target_bearing')./norm(pt_in_plane - target_bearing');

if 0
    min([norm(plane_norm_check + plane_norm) norm(plane_norm_check - plane_norm)])
    dot(pt_in_plane,plane_norm)
end

err = pt_in_plane./norm(pt_in_plane) - target_bearing';

d_err_d_pt = compute_d_bearing_d_pt_jac(pt_in_plane);

d_pt_d_plane_norm = compute_d_pt_d_plane_norm_jac(target_bearing', plane_norm);

d_plane_norm_d_plane = compute_d_bearing_d_pt_jac(plane);


[d_plane_d_R_host, d_plane_d_t_host, d_plane_d_R_target, d_plane_d_t_target] = compute_d_plane_d_pose(Twc_host, Twc_target, host_bearing');


d_err_d_plane = d_err_d_pt * d_pt_d_plane_norm * d_plane_norm_d_plane;

d_err_d_host_pose = d_err_d_plane * [d_plane_d_R_host d_plane_d_t_host];

d_err_d_target_pose = d_err_d_plane * [d_plane_d_R_target d_plane_d_t_target];

end
function d_bearing_d_pt = compute_d_bearing_d_pt_jac(pt)

d_bearing_d_pt = (norm(pt).*eye(3) - pt * pt'./norm(pt))./(norm(pt)^2);

end
function d_bearing_d_pt = compute_d_pt_d_plane_norm_jac(target_bearing, plane_norm)
d_bearing_d_pt = zeros(3,3);

d_bearing_d_pt(1,1) = -2*target_bearing(1)*plane_norm(1) - target_bearing(2)*plane_norm(2) - target_bearing(3)*plane_norm(3);
d_bearing_d_pt(1,2) = -target_bearing(2)*plane_norm(1);
d_bearing_d_pt(1,3) = -target_bearing(3)*plane_norm(1);

d_bearing_d_pt(2,1) = -target_bearing(1)*plane_norm(2);
d_bearing_d_pt(2,2) = -target_bearing(1)*plane_norm(1) - 2*target_bearing(2)*plane_norm(2) - target_bearing(3)*plane_norm(3);
d_bearing_d_pt(2,3) = -target_bearing(3)*plane_norm(2);

d_bearing_d_pt(3,1) = -target_bearing(1)*plane_norm(3);
d_bearing_d_pt(3,2) = -target_bearing(2)*plane_norm(3);
d_bearing_d_pt(3,3) = -target_bearing(1)*plane_norm(1) - target_bearing(2)*plane_norm(2) - 2*target_bearing(3)*plane_norm(3);

end
function [d_plane_d_R_host, d_plane_d_t_host, d_plane_d_R_target, d_plane_d_t_target] = compute_d_plane_d_pose(Twc_host, Twc_target, host_bearing)

Rwc_j = Twc_target(1:3,1:3);
Rwc_i = Twc_host(1:3,1:3);
twc_j = Twc_target(1:3,4);
twc_i = Twc_host(1:3,4);
R_cj_ci = Twc_target(1:3,1:3)' * Twc_host(1:3,1:3)';

d_plane_d_R_host = SkewSymMat(-R_cj_ci * host_bearing) * Rwc_j' * SkewSymMat(twc_i) - SkewSymMat(Rwc_j' * (twc_i - twc_j)) * Rwc_j' * SkewSymMat(Rwc_i * host_bearing);

d_plane_d_t_host = -SkewSymMat(R_cj_ci * host_bearing) * Rwc_j';

d_plane_d_R_target = -d_plane_d_R_host;

d_plane_d_t_target = -d_plane_d_t_host;



end