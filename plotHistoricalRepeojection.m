function plotHistoricalRepeojection()
%% | 1   |      2      |   3        |    4      |     5         |      6          |     7      |      8     |      9       |      10      |  
%% | pid |  host_fid   | target_fid | host_cid  |  target_cid   |   reproj_error  |  rho_raw   |   rho_filt |   rho_final  | sigma_bound  |

a = load('G:\matlab\data\direct\sim\reproject_statistics.txt');

raw_final_bound = abs((a(:,7) - a(:,9))./a(:,10));
filt_final_bound = abs((a(:,8) - a(:,9))./a(:,10));
figure,subplot(3,1,1);plot(a(:,6)); grid on; title('reprojection error');
subplot(3,1,2);plot(raw_final_bound); grid on; title('raw - bound ratio');
subplot(3,1,3);plot(filt_final_bound); grid on; title('filt - bound ratio');


invalid_rep_id = find(a(:,6) > 2);


invalid_raw_id = find(raw_final_bound > 20);
invalid_filt_id = find(filt_final_bound > 20);

aa = intersect(invalid_rep_id, invalid_raw_id);
bb = intersect(invalid_rep_id, invalid_filt_id);



end